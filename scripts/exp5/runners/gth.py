"""GTH (GenomeThreader) forced-locus runner."""
from __future__ import annotations

import re
from pathlib import Path

from ..common import run_subproc

GTH = "gth"
SPECIES_MODEL = {"mouse": "mouse", "dog": "human",
                 "chicken": "chicken", "zebrafish": "human"}


def run_one(protein_id, protein_seq, region_fa, region_chrom,
            region_offset, region_strand, scratch, options) -> dict | None:
    species = options.get("species", "mouse")
    sp_model = SPECIES_MODEL.get(species, "human")
    prot_fa = scratch / f"{protein_id}.prot.fa"
    prot_fa.write_text(f">{protein_id}\n{protein_seq}\n")
    cmd = [GTH, "-genomic", str(region_fa), "-protein", str(prot_fa),
           "-species", sp_model, "-prseedlength", "4", "-prminmatchlen", "8"]
    rc, stdout, _ = run_subproc(cmd, timeout=options.get("timeout", 180))
    prot_fa.unlink(missing_ok=True)
    if rc != 0 or not stdout.strip():
        return None
    return _parse(stdout, region_chrom, region_offset)


def _parse(text: str, region_chrom: str, region_offset: int) -> dict | None:
    hits = []
    current_exons: list[dict] = []
    current_score = 0.0
    current_strand = "+"

    def _flush():
        nonlocal current_exons, current_score, current_strand
        if current_exons:
            hits.append({
                "rank": 0, "chrom": region_chrom,
                "strand": current_strand, "score": current_score,
                "identity": 0,
                "query_start": 1, "query_end": 0,
                "cds": sorted(current_exons, key=lambda c: c["start"]),
            })
        current_exons = []
        current_score = 0.0
        current_strand = "+"

    for line in text.splitlines():
        line = line.strip()
        if line.startswith("MATCH"):
            _flush()
            parts = line.split("\t")
            if len(parts) >= 4:
                try:
                    current_score = float(parts[3])
                except ValueError:
                    pass
            continue
        m = re.match(r"PGS_\S+([+-])_\S+", line)
        if m:
            current_strand = m.group(1)
            continue
        m = re.match(r"Exon\s+\d+\s+(\d+)\s+(\d+)\s+\(\s*\d+\s+n\)", line)
        if m:
            s = int(m.group(1)); e = int(m.group(2))
            current_exons.append({
                "chrom": region_chrom,
                "start": min(s, e) + region_offset,
                "end":   max(s, e) + region_offset,
            })
    _flush()
    if not hits:
        return None
    hits.sort(key=lambda h: -h["score"])
    return hits[0]
