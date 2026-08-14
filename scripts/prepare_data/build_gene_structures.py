#!/usr/bin/env python3
"""Build gene-structure BED files from Ensembl release 113 GTFs.

Reads $RESCUE_DATA_ROOT/raw/release-113/<Species>.<assembly>.113.gtf.gz and
writes (under $RESCUE_GENE_STRUCTURE_DIR, default data/processed/gene_structures/):

  <species>_CDS.bed    - CDS features of protein_coding genes
                         chrom start end gene_id transcript_id strand phase
  <species>_exons.bed  - exon features of protein_coding genes
                         chrom start end gene_id transcript_id strand

Coordinates are converted from GTF 1-based inclusive to 0-based half-open.
Ensembl version suffixes are stripped from gene/transcript IDs.

Standard library only.
"""
from __future__ import annotations

import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import DATA_ROOT, GENE_STRUCTURES  # noqa: E402

SPECIES = {
    "mouse": ("Mus_musculus", "GRCm39"),
    "dog": ("Canis_lupus_familiaris", "ROS_Cfam_1.0"),
    "chicken": ("Gallus_gallus", "bGalGal1.mat.broiler.GRCg7b"),
    "zebrafish": ("Danio_rerio", "GRCz11"),
    "human": ("Homo_sapiens", "GRCh38"),
}


def strip_version(identifier: str) -> str:
    return identifier.rsplit(".", 1)[0] if "." in identifier else identifier


def parse_attrs(attr_str: str) -> dict[str, str]:
    """Parse Ensembl GTF attributes: `key1 "v1"; key2 "v2"; ...`."""
    out: dict[str, str] = {}
    for part in attr_str.split(";"):
        part = part.strip()
        if not part or " " not in part:
            continue
        key, _, val = part.partition(" ")
        out[key.strip()] = val.strip().strip('"')
    return out


def build_bed(gtf_path: Path, feature: str, with_phase: bool, out_path: Path) -> int:
    n = 0
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    with gzip.open(gtf_path, "rt") as fh, open(tmp, "w") as out:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != feature:
                continue
            attrs = parse_attrs(f[8])
            if attrs.get("gene_biotype") != "protein_coding":
                continue
            gene_id = strip_version(attrs.get("gene_id", ""))
            tx_id = strip_version(attrs.get("transcript_id", ""))
            if not gene_id or not tx_id:
                continue
            start = int(f[3]) - 1  # 1-based inclusive -> 0-based half-open
            end = int(f[4])
            row = [f[0], str(start), str(end), gene_id, tx_id, f[6]]
            if with_phase:
                try:
                    phase = str(int(f[7]))
                except ValueError:
                    phase = "0"
                row.append(phase)
            out.write("\t".join(row) + "\n")
            n += 1
    tmp.rename(out_path)
    return n


def main() -> int:
    raw = DATA_ROOT / "raw" / "release-113"
    GENE_STRUCTURES.mkdir(parents=True, exist_ok=True)
    for short, (spfile, asm) in SPECIES.items():
        gtf = raw / f"{spfile}.{asm}.113.gtf.gz"
        if not gtf.exists():
            print(f"SKIP {short}: missing {gtf} (run download_ensembl.sh first)")
            continue
        n_cds = build_bed(gtf, "CDS", True, GENE_STRUCTURES / f"{short}_CDS.bed")
        n_exon = build_bed(gtf, "exon", False, GENE_STRUCTURES / f"{short}_exons.bed")
        print(f"{short}: CDS rows={n_cds}, exon rows={n_exon}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
