#!/usr/bin/env python3
"""Build the small RESCUE test dataset (mouse, paralog regime).

Derives a self-contained mini data tree from the full benchmark data:
40 mouse paralog proteins, a mini genome containing each ground-truth /
alternative locus with flanking sequence, and all dependent files with
coordinates remapped to the mini genome. Used to produce the Zenodo test
dataset referenced by the README quick start.

Selection is deterministic: proteins with 2-6 target loci, all locus spans
<= 100 kb, evenly spaced over the primary-identity order.

Requires: pysam (for indexed FASTA access).

Usage:
  RESCUE_DATA_ROOT=/path/to/biodata RESCUE_V3_DATA=/path/to/data/v3 \
      python scripts/prepare_data/make_test_dataset.py --out test_dataset
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

SPECIES = "mouse"
ENSEMBL = "mus_musculus"
N_PROTEINS = 40
MAX_LOCI = 6
MAX_SPAN = 100_000
FLANK = 250_000


def region_name(chrom: str, start: int, end: int) -> str:
    return f"{chrom}_{start}_{end}"


def load_cds_by_gene(bed_path: Path):
    """gene -> list of (chrom, start, end, tx, strand, phase)."""
    genes: dict[str, list[tuple]] = defaultdict(list)
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            genes[f[3]].append((f[0], int(f[1]), int(f[2]), f[4], f[5],
                                f[6] if len(f) > 6 else "0"))
    return genes


def canonical_span(blocks: list[tuple]) -> tuple[str, int, int]:
    """(chrom, min_start, max_end) of the longest-total-CDS transcript."""
    by_tx: dict[str, list[tuple]] = defaultdict(list)
    for b in blocks:
        by_tx[b[3]].append(b)
    best_tx = max(sorted(by_tx), key=lambda t: sum(b[2] - b[1] for b in by_tx[t]))
    tx = by_tx[best_tx]
    return tx[0][0], min(b[1] for b in tx), max(b[2] for b in tx)


def main() -> int:
    ap = argparse.ArgumentParser(description="Build RESCUE test dataset")
    ap.add_argument("--out", required=True, help="output directory (data tree root)")
    ap.add_argument("--n", type=int, default=N_PROTEINS)
    ap.add_argument("--flank", type=int, default=FLANK)
    args = ap.parse_args()

    import pysam  # deferred: only this script needs it

    v3p = V3_DATA / f"{SPECIES}_paralog"
    gt = json.loads((v3p / "ground_truth.json").read_text())
    pinfo = json.loads((v3p / "paralog_info.json").read_text())
    cds = load_cds_by_gene(GENE_STRUCTURES / f"{SPECIES}_CDS.bed")

    # --- select proteins -------------------------------------------------
    candidates = []
    for pid, info in pinfo.items():
        targets = info.get("all_targets", [])
        if not (2 <= len(targets) <= MAX_LOCI):
            continue
        if any(t["target_gene"] not in cds for t in targets):
            continue
        spans = []
        ok = True
        for t in targets:
            chrom, s, e = canonical_span(cds[t["target_gene"]])
            if e - s > MAX_SPAN:
                ok = False
                break
            spans.append((chrom, s, e))
        if ok:
            candidates.append((pid, info["best_identity"], spans))
    candidates.sort(key=lambda c: c[1])
    step = len(candidates) / args.n
    selected = [candidates[int(i * step)] for i in range(args.n)]
    print(f"selected {len(selected)} of {len(candidates)} eligible proteins")
    selected_ids = [c[0] for c in selected]

    # --- build + merge regions -------------------------------------------
    fai = {
        line.split("\t")[0]: int(line.split("\t")[1])
        for line in open(str(GENOMES / f"{SPECIES}.fa") + ".fai")
    }
    intervals: dict[str, list[list[int]]] = defaultdict(list)
    for _, _, spans in selected:
        for chrom, s, e in spans:
            intervals[chrom].append([max(0, s - args.flank),
                                     min(fai[chrom], e + args.flank)])
    regions = []  # (name, chrom, start, end)
    for chrom, ivs in sorted(intervals.items()):
        ivs.sort()
        merged = []
        for s, e in ivs:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        for s, e in merged:
            regions.append((region_name(chrom, s, e), chrom, s, e))
    print(f"{len(regions)} merged regions, "
          f"{sum(e - s for _, _, s, e in regions) / 1e6:.1f} Mb total")

    def remap(chrom: str, pos_s: int, pos_e: int):
        """Map an interval to (region_name, rel_start, rel_end) or None."""
        for name, rchrom, rs, re_ in regions:
            if rchrom == chrom and rs <= pos_s and pos_e <= re_:
                return name, pos_s - rs, pos_e - rs
        return None

    # --- write mini genome ------------------------------------------------
    out = Path(args.out)
    gdir = out / "data/processed/genomes"
    gdir.mkdir(parents=True, exist_ok=True)
    fa = pysam.FastaFile(str(GENOMES / f"{SPECIES}.fa"))
    with open(gdir / f"{SPECIES}.fa", "w") as fh:
        for name, chrom, s, e in regions:
            seq = fa.fetch(chrom, s, e)
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    # --- remap v3 ground truth / paralog info -----------------------------
    v3out = out / f"data/v3/{SPECIES}_paralog"
    v3out.mkdir(parents=True, exist_ok=True)
    gt_out = {}
    for pid in selected_ids:
        entry = gt[pid]
        new_exons = []
        new_chrom = None
        for ex in entry["exons"]:
            m = remap(ex["chrom"], ex["start"], ex["end"])
            if m is None:
                raise SystemExit(f"GT exon outside regions: {pid} {ex}")
            new_chrom = m[0]
            new_exons.append({"chrom": m[0], "start": m[1], "end": m[2]})
        gt_out[pid] = {**entry, "chrom": new_chrom, "exons": new_exons}
        for locus in gt_out[pid].get("paralog_loci", []):
            gene_blocks = cds[locus["target_gene"]]
            chrom, s, e = canonical_span(gene_blocks)
            m = remap(chrom, s, e)
            if m is None:
                raise SystemExit(f"alt locus outside regions: {pid} {locus}")
            locus["chrom"] = m[0]
    (v3out / "ground_truth.json").write_text(json.dumps(gt_out, indent=2))

    pinfo_out = {}
    for pid in selected_ids:
        info = pinfo[pid]
        for t in info["all_targets"]:
            chrom, s, e = canonical_span(cds[t["target_gene"]])
            t["chrom"] = remap(chrom, s, e)[0]
        pinfo_out[pid] = info
    (v3out / "paralog_info.json").write_text(json.dumps(pinfo_out, indent=2))

    # proteins.fa and pairs.tsv subset (sequences/IDs unchanged)
    keep = set(selected_ids)
    with open(v3p / "proteins.fa") as fin, open(v3out / "proteins.fa", "w") as fout:
        write = False
        for line in fin:
            if line.startswith(">"):
                write = line[1:].split()[0] in keep
            if write:
                fout.write(line)
    with open(v3p / "pairs.tsv") as fin, open(v3out / "pairs.tsv", "w") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames,
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            if row["human_protein"] in keep:
                writer.writerow(row)
    (v3out / "meta.json").write_text(json.dumps({
        "species_short": SPECIES,
        "species_ensembl": ENSEMBL,
        "genome": "data/processed/genomes/mouse.fa",
        "n_proteins": len(selected_ids),
        "benchmark_type": "paralog_rich",
        "test_dataset": True,
    }, indent=2))

    # --- CDS BED subset (genes fully contained in a region) ---------------
    gsdir = out / "data/processed/gene_structures"
    gsdir.mkdir(parents=True, exist_ok=True)
    kept_genes = set()
    with open(gsdir / f"{SPECIES}_CDS.bed", "w") as fh:
        for gene, blocks in sorted(cds.items()):
            m_first = remap(blocks[0][0], blocks[0][1], blocks[0][2])
            if m_first is None:
                continue
            name = m_first[0]
            rs = None
            for rname, rchrom, rstart, rend in regions:
                if rname == name:
                    rs = rstart
                    break
            if all(remap(b[0], b[1], b[2]) is not None and
                   remap(b[0], b[1], b[2])[0] == name for b in blocks):
                for b in blocks:
                    m = remap(b[0], b[1], b[2])
                    fh.write(f"{m[0]}\t{m[1]}\t{m[2]}\t{gene}\t{b[3]}\t{b[4]}\t{b[5]}\n")
                kept_genes.add(gene)

    # --- proteome subset --------------------------------------------------
    pdir = out / "data/processed/proteins"
    pdir.mkdir(parents=True, exist_ok=True)
    with open(PROTEOMES / f"{SPECIES}_canonical.fa") as fin, \
            open(pdir / f"{SPECIES}_canonical.fa", "w") as fout:
        write = False
        for line in fin:
            if line.startswith(">"):
                gene = ""
                for field in line[1:].split():
                    if field.startswith("gene:"):
                        gene = field.split(":", 1)[1].rsplit(".", 1)[0]
                        break
                write = gene in kept_genes
            if write:
                fout.write(line)

    # --- ortholog table subset ---------------------------------------------
    odir = out / "data/processed/orthologs"
    odir.mkdir(parents=True, exist_ok=True)
    human_genes = {pinfo[pid]["human_gene"] for pid in selected_ids}
    with open(ORTHOLOGS / f"human_vs_{ENSEMBL}_one2one.tsv") as fin, \
            open(odir / f"human_vs_{ENSEMBL}_one2one.tsv", "w") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames,
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            g1 = row["gene_stable_id"].rsplit(".", 1)[0]
            g2 = row["homology_gene_stable_id"].rsplit(".", 1)[0]
            if g1 in human_genes or g2 in human_genes or \
                    g1 in kept_genes or g2 in kept_genes:
                writer.writerow(row)

    print(f"test dataset written to {out}")
    print(f"kept {len(kept_genes)} genes in CDS/proteome subsets")
    return 0


if __name__ == "__main__":
    sys.exit(main())
