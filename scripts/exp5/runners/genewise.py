"""GeneWise (wise2) forced-locus runner."""
from __future__ import annotations

import os
import sys
from pathlib import Path

from ..common import run_subproc
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from rescue_paths import binary  # noqa: E402

GENEWISE = binary("GENEWISE_BIN", "genewise")
WISECONFIGDIR = os.environ.get("WISECONFIGDIR", "")


def run_one(protein_id, protein_seq, region_fa, region_chrom,
            region_offset, region_strand, scratch, options) -> dict | None:
    prot_fa = scratch / f"{protein_id}.prot.fa"
    prot_fa.write_text(f">{protein_id}\n{protein_seq}\n")

    if region_strand == "+":
        strand_flag = "-tfor"
    elif region_strand == "-":
        strand_flag = "-trev"
    else:
        strand_flag = "-both"

    env = os.environ.copy()
    if WISECONFIGDIR:
        env["WISECONFIGDIR"] = WISECONFIGDIR

    cmd = [GENEWISE, "-gff", strand_flag, str(prot_fa), str(region_fa)]
    rc, stdout, _ = run_subproc(cmd, timeout=options.get("timeout", 300), env=env)
    prot_fa.unlink(missing_ok=True)
    if rc != 0 or not stdout.strip():
        return None
    return _parse(stdout, region_chrom, region_offset)


def _parse(text: str, region_chrom: str, region_offset: int) -> dict | None:
    hits = []
    cur_cds: list[dict] = []
    cur_score = 0.0
    cur_strand = "+"

    def _flush():
        nonlocal cur_cds, cur_score, cur_strand
        if cur_cds:
            hits.append({
                "rank": 0, "chrom": region_chrom,
                "strand": cur_strand, "score": cur_score,
                "identity": 0,
                "query_start": 1, "query_end": 0,
                "cds": sorted(cur_cds, key=lambda c: c["start"]),
            })
        cur_cds = []
        cur_score = 0.0
        cur_strand = "+"

    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line == "//":
            _flush()
            continue
        parts = line.split("\t")
        if len(parts) < 8:
            continue
        feat = parts[2]
        strand = parts[6]
        if feat == "match":
            try:
                cur_score = float(parts[5])
                cur_strand = strand
            except ValueError:
                pass
        elif feat in ("cds", "CDS"):
            try:
                s = int(parts[3]); e = int(parts[4])
                cur_cds.append({
                    "chrom": region_chrom,
                    "start": min(s, e) + region_offset,
                    "end":   max(s, e) + region_offset,
                })
                cur_strand = strand
            except ValueError:
                pass
    _flush()
    if not hits:
        return None
    hits.sort(key=lambda h: -h["score"])
    return hits[0]
