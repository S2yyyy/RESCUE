#!/usr/bin/env python3
"""Build the paralog benchmark query set (data/v3/<species>_paralog/).

For each target species:
  1. Load human<->target ortholog_one2many / ortholog_many2many pairs
  2. Keep human genes with n_target_loci >= 2
  3. Primary ground truth = highest-identity target ortholog; all other
     target loci are recorded as acceptable alternatives (paralog_info.json)
  4. Sort candidates by paralog count (descending, ties by identity) and
     take the top 500 per species (deterministic, no RNG)

Outputs per species: proteins.fa, ground_truth.json, paralog_info.json,
pairs.tsv, meta.json.

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

MAX_PARALOG = 500


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


def load_paralog_pairs(species_ensembl):
    """Load paralog-rich pairs, grouped by human gene."""
    orth_file = ORTHOLOGS / f"human_vs_{species_ensembl}_paralog_rich.tsv"
    gene_groups = defaultdict(list)
    with open(orth_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            human_gene = row["gene_stable_id"]
            human_prot = row["protein_stable_id"]
            target_gene = row["homology_gene_stable_id"]
            identity = float(row["identity"]) if row.get("identity") else 0
            n_target_loci = int(row["n_target_loci"]) if row.get("n_target_loci") else 0

            human_gene = human_gene.rsplit(".", 1)[0] if "." in human_gene else human_gene
            human_prot = human_prot.rsplit(".", 1)[0] if "." in human_prot else human_prot
            target_gene = target_gene.rsplit(".", 1)[0] if "." in target_gene else target_gene

            gene_groups[human_gene].append({
                "human_gene": human_gene,
                "human_protein": human_prot,
                "target_gene": target_gene,
                "identity": identity,
                "n_target_loci": n_target_loci,
            })
    return gene_groups


def prepare_species(species_short, species_ensembl, human_prots, human_gene_to_prot, max_pairs):
    print(f"\n{'=' * 60}\nPreparing {species_short} ({species_ensembl}) — paralog\n{'=' * 60}")

    out_dir = V3_DATA / f"{species_short}_paralog"
    out_dir.mkdir(parents=True, exist_ok=True)

    target_cds = load_cds_bed(GENE_STRUCTURES / f"{species_short}_CDS.bed")
    print(f"  Target CDS: {len(target_cds)} genes")

    gene_groups = load_paralog_pairs(species_ensembl)
    print(f"  Human genes with paralog-rich orthologs: {len(gene_groups)}")

    filtered_groups = {
        hg: orths for hg, orths in gene_groups.items()
        if any(o["n_target_loci"] >= 2 for o in orths)
    }
    print(f"  Genes with n_target_loci >= 2: {len(filtered_groups)}")

    candidates = []
    for hg, orthologs in filtered_groups.items():
        hp = orthologs[0]["human_protein"]
        if hp not in human_prots and hg in human_gene_to_prot:
            hp = human_gene_to_prot[hg]
        if hp not in human_prots:
            continue

        valid_orthologs = [o for o in orthologs if o["target_gene"] in target_cds]
        if not valid_orthologs:
            continue

        best = max(valid_orthologs, key=lambda o: o["identity"])

        alt_loci = []
        for o in valid_orthologs:
            if o["target_gene"] != best["target_gene"]:
                tx_id = pick_canonical_transcript(target_cds[o["target_gene"]])
                if tx_id is None:
                    continue
                exons = target_cds[o["target_gene"]][tx_id]
                sorted_exons = sorted(exons, key=lambda e: e["start"])
                alt_loci.append({
                    "target_gene": o["target_gene"],
                    "identity": o["identity"],
                    "chrom": sorted_exons[0]["chrom"],
                    "strand": sorted_exons[0]["strand"],
                })

        if not alt_loci:
            continue

        candidates.append({
            "human_gene": hg,
            "human_protein": hp,
            "best": best,
            "alt_loci": alt_loci,
            "n_paralogs": len(valid_orthologs),
            "n_target_loci": best["n_target_loci"],
        })
    print(f"  Valid candidates (protein + CDS + alt loci): {len(candidates)}")

    # Deterministic selection: most paralogs first, ties by identity.
    candidates.sort(key=lambda c: (-c["n_paralogs"], c["best"]["identity"]))
    selected = candidates[:max_pairs] if len(candidates) > max_pairs else candidates
    print(f"  Selected: {len(selected)} proteins")

    ground_truth, paralog_info, pairs_out, protein_ids = {}, {}, [], []
    for c in selected:
        hp = c["human_protein"]
        best = c["best"]
        tg = best["target_gene"]

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
            "identity": best["identity"],
            "paralog_loci": c["alt_loci"],
        }
        paralog_info[hp] = {
            "human_gene": c["human_gene"],
            "n_target_loci": c["n_target_loci"],
            "n_paralogs_with_cds": c["n_paralogs"],
            "best_target_gene": tg,
            "best_identity": best["identity"],
            "all_targets": [
                {"target_gene": tg, "identity": best["identity"],
                 "chrom": sorted_exons[0]["chrom"], "is_best": True}
            ] + [
                {"target_gene": al["target_gene"], "identity": al["identity"],
                 "chrom": al["chrom"], "is_best": False}
                for al in c["alt_loci"]
            ],
        }
        pairs_out.append({
            "human_protein": hp,
            "human_gene": c["human_gene"],
            "target_gene": tg,
            "identity": best["identity"],
        })
        protein_ids.append(hp)

    with open(out_dir / "proteins.fa", "w") as fh:
        for pid in protein_ids:
            fh.write(f">{pid}\n{human_prots[pid]}\n")
    with open(out_dir / "ground_truth.json", "w") as fh:
        json.dump(ground_truth, fh, indent=2)
    with open(out_dir / "paralog_info.json", "w") as fh:
        json.dump(paralog_info, fh, indent=2)
    with open(out_dir / "pairs.tsv", "w") as fh:
        writer = csv.DictWriter(fh, delimiter="\t",
                                fieldnames=["human_protein", "human_gene",
                                            "target_gene", "identity"])
        writer.writeheader()
        for p in pairs_out:
            writer.writerow(p)

    meta = {
        "species_short": species_short,
        "species_ensembl": species_ensembl,
        "genome": str(GENOMES / f"{species_short}.fa"),
        "genome_index": str(GENOMES / f"{species_short}.mpi"),
        "n_proteins": len(protein_ids),
        "benchmark_type": "paralog_rich",
    }
    with open(out_dir / "meta.json", "w") as fh:
        json.dump(meta, fh, indent=2)
    return len(selected)


def main() -> int:
    parser = argparse.ArgumentParser(description="Prepare paralog benchmark data")
    parser.add_argument("--species", nargs="*", default=list(SPECIES.keys()))
    parser.add_argument("--max-pairs", type=int, default=MAX_PARALOG)
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

    print(f"\nTotal: {total} paralog proteins across {len(args.species)} species")
    print(f"Data written to {V3_DATA}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
