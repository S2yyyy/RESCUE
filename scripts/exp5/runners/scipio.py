"""Scipio forced-locus runner.

Runs scipio.1.4.1.pl in the `scipio` conda env, then yaml2gff.1.4.pl to
extract per-exon coords. Coordinates in the GFF are 1-based inclusive
on the region; we offset to genome coordinates.
"""
from __future__ import annotations

import re
import os
import sys
from pathlib import Path

from ..common import run_subproc
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from rescue_paths import PERL, TOOLS_ROOT  # noqa: E402

SCIPIO_ROOT = Path(os.environ.get("SCIPIO_ROOT", TOOLS_ROOT / "scipio/scipio-1.4"))
SCIPIO_PL = str(Path(os.environ.get("SCIPIO_PL", SCIPIO_ROOT / "scipio.1.4.1.pl")))
YAML2GFF = str(Path(os.environ.get("YAML2GFF_PL", SCIPIO_ROOT / "yaml2gff.1.4.pl")))
SCIPIO_PERL = os.environ.get("SCIPIO_PERL", PERL)


def run_one(protein_id, protein_seq, region_fa, region_chrom,
            region_offset, region_strand, scratch, options) -> dict | None:
    prot_fa = scratch / f"{protein_id}.prot.fa"
    prot_fa.write_text(f">{protein_id}\n{protein_seq}\n")
    yaml_path = scratch / f"{protein_id}.scipio.yaml"
    cmd = [SCIPIO_PERL, SCIPIO_PL,
           "--min_score=0.3", "--min_identity=50",
           "--min_coverage=0.3",
           str(region_fa), str(prot_fa)]
    rc, stdout, _ = run_subproc(cmd, timeout=options.get("timeout", 180))
    prot_fa.unlink(missing_ok=True)
    if rc != 0 or not stdout.strip():
        return None
    yaml_path.write_text(stdout)
    rc2, gff_stdout, _ = run_subproc(
        [SCIPIO_PERL, YAML2GFF, str(yaml_path)],
        timeout=options.get("timeout", 60))
    yaml_path.unlink(missing_ok=True)
    if rc2 != 0 or not gff_stdout.strip():
        return None
    return _parse(gff_stdout, region_chrom, region_offset)


_CHROM_UNESC = re.compile(r"%([0-9A-Fa-f]{2})")


def _unescape(s: str) -> str:
    return _CHROM_UNESC.sub(lambda m: chr(int(m.group(1), 16)), s)


def _parse(gff_text: str, region_chrom: str, region_offset: int) -> dict | None:
    cds: list[dict] = []
    score = 0.0
    strand = "+"
    qs, qe = 1, 0
    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        feat = f[2]
        if feat != "protein_match":
            continue
        try:
            s = int(f[3]); e = int(f[4])
            sc = float(f[5])
        except ValueError:
            continue
        strand = f[6]
        score = max(score, sc)
        cds.append({
            "chrom": region_chrom,
            "start": min(s, e) + region_offset,
            "end":   max(s, e) + region_offset,
        })
        # Query coords from attributes: Query=PROT_ID Q_S Q_E
        attrs = f[8]
        m = re.search(r"Query=\S+\s+(\d+)\s+(\d+)", attrs)
        if m:
            cqs, cqe = int(m.group(1)), int(m.group(2))
            qs = min(qs, cqs) if qe > 0 else cqs
            qe = max(qe, cqe)
    if not cds:
        return None
    return {
        "rank": 1, "chrom": region_chrom, "strand": strand,
        "score": score, "identity": 0,
        "query_start": qs, "query_end": qe,
        "cds": sorted(cds, key=lambda c: c["start"]),
    }
