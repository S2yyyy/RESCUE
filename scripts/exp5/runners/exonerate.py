"""Exonerate forced-locus runner."""
from __future__ import annotations

from pathlib import Path

from ..common import run_subproc

EXONERATE = "exonerate"


def run_one(protein_id, protein_seq, region_fa, region_chrom,
            region_offset, region_strand, scratch, options) -> dict | None:
    prot_fa = scratch / f"{protein_id}.prot.fa"
    prot_fa.write_text(f">{protein_id}\n{protein_seq}\n")
    cmd = [EXONERATE,
           "--model", "protein2genome",
           "--bestn", "1",
           "--showtargetgff", "yes",
           "--showvulgar", "yes",
           "--verbose", "0",
           "--showalignment", "no",
           str(prot_fa), str(region_fa)]
    rc, stdout, _ = run_subproc(cmd, timeout=options.get("timeout", 300))
    prot_fa.unlink(missing_ok=True)
    if rc != 0 or not stdout.strip():
        return None
    return _parse(stdout, region_chrom, region_offset)


def _parse(text: str, region_chrom: str, region_offset: int) -> dict | None:
    hits = []
    cur_hit = None
    cur_cds: list[dict] = []

    def _flush():
        nonlocal cur_hit, cur_cds
        if cur_hit and cur_cds:
            cur_hit["cds"] = sorted(cur_cds, key=lambda c: c["start"])
            hits.append(cur_hit)
        cur_hit = None
        cur_cds = []

    for line in text.splitlines():
        line = line.strip()
        if line.startswith("vulgar:"):
            _flush()
            parts = line.split()
            if len(parts) >= 10:
                cur_hit = {
                    "rank": 0, "chrom": region_chrom,
                    "strand": "+" if parts[8] == "+" else "-",
                    "score": float(parts[9]),
                    "identity": 0,
                    "query_start": int(parts[2]) + 1,
                    "query_end": int(parts[3]),
                }
        elif "\texonerate:" in line and cur_hit:
            parts = line.split("\t")
            if len(parts) >= 5 and parts[2] == "cds":
                cur_cds.append({
                    "chrom": region_chrom,
                    "start": int(parts[3]) + region_offset,
                    "end":   int(parts[4]) + region_offset,
                })
    _flush()
    if not hits:
        return None
    hits.sort(key=lambda h: -h["score"])
    return hits[0]
