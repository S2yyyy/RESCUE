#!/usr/bin/env python3
"""Build per-species canonical proteomes from Ensembl release 113.

For each protein-coding gene, the canonical protein is the translation of
the transcript with the longest total CDS length (ties broken by the
lexicographically smallest transcript ID, for determinism).

Inputs:
  $RESCUE_DATA_ROOT/raw/release-113/<Species>.<assembly>.pep.all.fa.gz
  $RESCUE_GENE_STRUCTURE_DIR/<species>_CDS.bed   (run build_gene_structures.py first)

Output (under $RESCUE_PROTEOME_DIR, default data/processed/proteins/):
  <species>_canonical.fa   - full original pep.all header lines and sequences

Standard library only.
"""
from __future__ import annotations

import gzip
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import DATA_ROOT, GENE_STRUCTURES, PROTEOMES  # noqa: E402

SPECIES = {
    "mouse": ("Mus_musculus", "GRCm39"),
    "dog": ("Canis_lupus_familiaris", "ROS_Cfam_1.0"),
    "chicken": ("Gallus_gallus", "bGalGal1.mat.broiler.GRCg7b"),
    "zebrafish": ("Danio_rerio", "GRCz11"),
    "human": ("Homo_sapiens", "GRCh38"),
}


def strip_version(identifier: str) -> str:
    return identifier.rsplit(".", 1)[0] if "." in identifier else identifier


def canonical_transcripts(cds_bed: Path) -> set[str]:
    """{gene: {tx: total_cds_len}} -> set of canonical (longest-CDS) tx IDs."""
    gene_tx_len: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    with open(cds_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            gene_tx_len[f[3]][f[4]] += int(f[2]) - int(f[1])
    chosen = set()
    for gene, tx_lens in gene_tx_len.items():
        best = max(sorted(tx_lens), key=lambda t: tx_lens[t])
        chosen.add(best)
    return chosen


def filter_proteome(pep_gz: Path, keep_txs: set[str], out_path: Path) -> int:
    n = 0
    tmp = out_path.with_suffix(out_path.suffix + ".tmp")
    with gzip.open(pep_gz, "rt") as fh, open(tmp, "w") as out:
        keep_record = False
        for line in fh:
            if line.startswith(">"):
                tx = ""
                for field in line[1:].split():
                    if field.startswith("transcript:"):
                        tx = strip_version(field.split(":", 1)[1])
                        break
                keep_record = tx in keep_txs
                if keep_record:
                    n += 1
            if keep_record:
                out.write(line)
    tmp.rename(out_path)
    return n


def main() -> int:
    raw = DATA_ROOT / "raw" / "release-113"
    PROTEOMES.mkdir(parents=True, exist_ok=True)
    for short, (spfile, asm) in SPECIES.items():
        pep = raw / f"{spfile}.{asm}.pep.all.fa.gz"
        cds = GENE_STRUCTURES / f"{short}_CDS.bed"
        if not pep.exists() or not cds.exists():
            print(f"SKIP {short}: need {pep.name} and {cds.name}")
            continue
        keep = canonical_transcripts(cds)
        n = filter_proteome(pep, keep, PROTEOMES / f"{short}_canonical.fa")
        print(f"{short}: canonical transcripts={len(keep)}, proteins written={n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
