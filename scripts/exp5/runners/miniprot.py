"""miniprot forced-locus runner.

Runs `miniprot --gff -N 1` on a single-record region FASTA with one query.
Returns a single hit dict (best by score) or None on failure.
"""
from __future__ import annotations

import re
from pathlib import Path

from ..common import run_subproc

MINIPROT_BIN = "miniprot"


def run_one(protein_id: str, protein_seq: str, region_fa: Path,
            region_chrom: str, region_offset: int, region_strand: str,
            scratch: Path, options: dict) -> dict | None:
    """Run miniprot on a single protein × region.

    Returns a hit dict in the v3_evaluate-compatible schema, or None on failure.
    Coordinates in the returned `cds` are GENOME-ABSOLUTE 1-based inclusive.
    """
    prot_fa = scratch / f"{protein_id}.prot.fa"
    prot_fa.write_text(f">{protein_id}\n{protein_seq}\n")

    extra = options.get("miniprot_extra", [])
    cmd = [MINIPROT_BIN, "--gff", "-t", "1", "-N", "1", *extra,
           str(region_fa), str(prot_fa)]
    rc, stdout, stderr = run_subproc(cmd, timeout=options.get("timeout", 60))
    prot_fa.unlink(missing_ok=True)
    if rc != 0 or not stdout.strip():
        return None

    return _parse(stdout, region_chrom, region_offset)


def _parse(gff_text: str, region_chrom: str, region_offset: int) -> dict | None:
    """Parse miniprot GFF (single-region, single-query). Return best hit or None.

    miniprot emits the region name (e.g. "12:1000-2000") in the chrom column;
    we rewrite to the actual genome chromosome and add the offset.
    """
    cur = None
    cur_cds: list[dict] = []
    best: dict | None = None
    best_score = -1e18

    def _flush():
        nonlocal best, best_score
        if cur is None or not cur_cds:
            return
        sc = cur.get("score", 0)
        if sc > best_score:
            best_score = sc
            cur["cds"] = [
                {"chrom": region_chrom,
                 "start": int(c["start"]) + region_offset,
                 "end": int(c["end"]) + region_offset}
                for c in sorted(cur_cds, key=lambda c: int(c["start"]))
            ]
            cur["chrom"] = region_chrom
            best = cur.copy()

    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            if line.startswith("##PAF"):
                _flush()
                cur = None
                cur_cds = []
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        feat = f[2]
        if feat == "mRNA":
            _flush()
            cur_cds = []
            attr = dict(a.split("=", 1) for a in f[8].split(";") if "=" in a)
            tgt = attr.get("Target", "").split()
            qs = int(tgt[1]) if len(tgt) > 1 else 1
            qe = int(tgt[2]) if len(tgt) > 2 else 0
            cur = {
                "rank": int(attr.get("Rank", 1)),
                "strand": f[6],
                "score": float(f[5]) if f[5] != "." else 0.0,
                "identity": float(attr.get("Identity", 0)),
                "query_start": qs,
                "query_end": qe,
            }
        elif feat == "CDS" and cur is not None:
            cur_cds.append({"start": f[3], "end": f[4]})
    _flush()

    if best is None:
        return None
    return best
