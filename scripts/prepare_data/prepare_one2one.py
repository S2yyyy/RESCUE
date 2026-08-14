#!/usr/bin/env python3
"""Build the one-to-one benchmark query set (data/v3/<species>/).

For each target species:
  1. Load human<->target ortholog_one2one pairs (Compara TSV)
  2. Keep pairs with the human protein in the canonical proteome and the
     target gene present in the CDS BED
  3. Sort by protein identity and select up to 2000 pairs by evenly spaced
     sampling over the identity-sorted list (deterministic, no RNG)
  4. Ground truth = CDS exons of the target gene's canonical transcript
     (longest total CDS), 0-based half-open
  5. Stratification labels: identity tier, micro_exon, long_intron

Outputs per species: proteins.fa, ground_truth.json, strata.json, pairs.tsv,
meta.json.

Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import GENE_STRUCTURES, GENOMES, ORTHOLOGS, PROTEOMES, V3_DATA  # noqa: E402

SPECIES = {
    "mouse": "mus_musculus",
    "dog": "canis_lupus_familiaris",
    "chicken": "gallus_gallus",
    "zebrafish": "danio_rerio",
}

MAX_ONE2ONE = 2000


def load_fasta(fa_path, strip_version=True):
    """Load FASTA -> {id: seq} and {gene_id: prot_id}."""
    seqs, gene_to_prot = {}, {}
    cur_id, cur_gene, cur_seq = "", "", []
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(cur_seq)
                    if cur_gene:
                        gene_to_prot[cur_gene] = cur_id
                header = line[1:].strip()
                raw_id = header.split()[0]
                cur_id = raw_id.rsplit(".", 1)[0] if strip_version else raw_id
                cur_seq, cur_gene = [], ""
                for field in header.split():
                    if field.startswith("gene:"):
                        cur_gene = field.split(":")[1]
                        if strip_version:
                            cur_gene = cur_gene.rsplit(".", 1)[0]
            else:
                cur_seq.append(line.strip())
    if cur_id:
        seqs[cur_id] = "".join(cur_seq)
        if cur_gene:
            gene_to_prot[cur_gene] = cur_id
    return seqs, gene_to_prot


def load_cds_bed(bed_path):
    """Load CDS BED -> {gene_id: {transcript_id: [exon dicts]}}."""
    cds = defaultdict(lambda: defaultdict(list))
    with open(bed_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            chrom = parts[0]
            start, end = int(parts[1]), int(parts[2])
            gene_id = parts[3]
            tx_id = parts[4]
            strand = parts[5]
            phase = int(parts[6]) if len(parts) > 6 else 0
            cds[gene_id][tx_id].append({
                "chrom": chrom, "start": start, "end": end,
                "strand": strand, "phase": phase,
            })
    return cds


def pick_canonical_transcript(tx_dict):
    """Pick the transcript with the longest total CDS."""
    best_tx, best_len = None, 0
    for tx_id, exons in tx_dict.items():
        total = sum(e["end"] - e["start"] for e in exons)
        if total > best_len:
            best_tx, best_len = tx_id, total
    return best_tx


def compute_strata(exons, identity):
    """Compute stratification labels for a protein."""
    labels = []
    if identity >= 80:
        labels.append("1:1_easy")
    elif identity < 60:
        labels.append("1:1_distant")
    else:
        labels.append("1:1_moderate")

    for e in exons:
        if e["end"] - e["start"] < 30:
            labels.append("micro_exon")
            break

    sorted_exons = sorted(exons, key=lambda e: e["start"])
    max_intron = 0
    for i in range(1, len(sorted_exons)):
        intron = sorted_exons[i]["start"] - sorted_exons[i - 1]["end"]
        max_intron = max(max_intron, intron)
    if max_intron > 100_000:
        labels.append("long_intron")

    return labels


def load_ortholog_pairs(species_ensembl):
    """Load 1:1 ortholog pairs from the Compara TSV."""
    orth_file = ORTHOLOGS / f"human_vs_{species_ensembl}_one2one.tsv"
    pairs = []
    with open(orth_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["species"] == "homo_sapiens":
                human_gene = row["gene_stable_id"]
                human_prot = row["protein_stable_id"]
                target_gene = row["homology_gene_stable_id"]
                identity = float(row["homology_identity"]) if row.get("homology_identity") else 0
            else:
                target_gene = row["gene_stable_id"]
                human_gene = row["homology_gene_stable_id"]
                human_prot = row["homology_protein_stable_id"]
                identity = float(row["identity"]) if row.get("identity") else 0
            human_gene = human_gene.rsplit(".", 1)[0] if "." in human_gene else human_gene
            human_prot = human_prot.rsplit(".", 1)[0] if "." in human_prot else human_prot
            target_gene = target_gene.rsplit(".", 1)[0] if "." in target_gene else target_gene
            pairs.append({
                "human_gene": human_gene,
                "human_protein": human_prot,
                "target_gene": target_gene,
                "identity": identity,
            })
    return pairs


def prepare_species(species_short, species_ensembl, human_prots, human_gene_to_prot, max_pairs):
    print(f"\n{'=' * 60}\nPreparing {species_short} ({species_ensembl})\n{'=' * 60}")

    out_dir = V3_DATA / species_short
    out_dir.mkdir(parents=True, exist_ok=True)

    target_cds = load_cds_bed(GENE_STRUCTURES / f"{species_short}_CDS.bed")
    print(f"  Target CDS: {len(target_cds)} genes")

    pairs = load_ortholog_pairs(species_ensembl)
    print(f"  Ortholog pairs: {len(pairs)}")

    valid_pairs = []
    for p in pairs:
        hp = p["human_protein"]
        tg = p["target_gene"]
        if hp not in human_prots and p["human_gene"] in human_gene_to_prot:
            hp = human_gene_to_prot[p["human_gene"]]
            p["human_protein"] = hp
        if hp in human_prots and tg in target_cds:
            valid_pairs.append(p)
    print(f"  Valid pairs (protein + CDS): {len(valid_pairs)}")

    # Deterministic stratified sampling: evenly spaced over identity order.
    valid_pairs.sort(key=lambda p: p["identity"])
    if len(valid_pairs) > max_pairs:
        step = len(valid_pairs) / max_pairs
        selected = [valid_pairs[int(i * step)] for i in range(max_pairs)]
    else:
        selected = valid_pairs
    print(f"  Selected: {len(selected)} pairs")

    ground_truth, strata, protein_ids = {}, {}, []
    for p in selected:
        hp = p["human_protein"]
        tg = p["target_gene"]
        tx_id = pick_canonical_transcript(target_cds[tg])
        exons = target_cds[tg][tx_id]
        sorted_exons = sorted(exons, key=lambda e: e["start"])

        ground_truth[hp] = {
            "target_gene": tg,
            "transcript": tx_id,
            "chrom": sorted_exons[0]["chrom"],
            "strand": sorted_exons[0]["strand"],
            "exons": [{"chrom": e["chrom"], "start": e["start"], "end": e["end"]}
                      for e in sorted_exons],
            "identity": p["identity"],
        }
        strata[hp] = compute_strata(sorted_exons, p["identity"])
        protein_ids.append(hp)

    strata_counts = defaultdict(int)
    for labels in strata.values():
        for label in labels:
            strata_counts[label] += 1
    print("  Strata distribution:")
    for s, c in sorted(strata_counts.items()):
        print(f"    {s}: {c}")

    with open(out_dir / "proteins.fa", "w") as fh:
        for pid in protein_ids:
            fh.write(f">{pid}\n{human_prots[pid]}\n")
    with open(out_dir / "ground_truth.json", "w") as fh:
        json.dump(ground_truth, fh, indent=2)
    with open(out_dir / "strata.json", "w") as fh:
        json.dump(strata, fh, indent=2)
    with open(out_dir / "pairs.tsv", "w") as fh:
        writer = csv.DictWriter(fh, delimiter="\t",
                                fieldnames=["human_protein", "human_gene",
                                            "target_gene", "identity"])
        writer.writeheader()
        for p in selected:
            writer.writerow(p)

    meta = {
        "species_short": species_short,
        "species_ensembl": species_ensembl,
        "genome": str(GENOMES / f"{species_short}.fa"),
        "genome_index": str(GENOMES / f"{species_short}.mpi"),
        "n_proteins": len(protein_ids),
        "n_1to1": len(selected),
    }
    with open(out_dir / "meta.json", "w") as fh:
        json.dump(meta, fh, indent=2)
    return len(selected)


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare one-to-one benchmark data")
    parser.add_argument("--species", nargs="*", default=list(SPECIES.keys()))
    parser.add_argument("--max-pairs", type=int, default=MAX_ONE2ONE)
    args = parser.parse_args()

    print("Loading human proteins...")
    human_prots, human_gene_to_prot = load_fasta(PROTEOMES / "human_canonical.fa")
    print(f"  Loaded {len(human_prots)} human proteins")

    total = 0
    for sp in args.species:
        if sp not in SPECIES:
            print(f"Unknown species: {sp}")
            continue
        total += prepare_species(sp, SPECIES[sp], human_prots,
                                 human_gene_to_prot, args.max_pairs)

    print(f"\nTotal: {total} proteins across {len(args.species)} species")
    print(f"Data written to {V3_DATA}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
