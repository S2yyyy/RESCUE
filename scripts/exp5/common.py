"""Shared utilities for Exp 5 forced-locus runners.

Provides:
  - GT region extraction from data/v3/<sp>/ground_truth.json
  - Protein FASTA / genome FASTA loading
  - Per-tool subprocess invocation helpers
  - Output schema (per_protein TSV columns, hit dict schema)

The "hit" dict produced by every runner must be compatible with
v3_evaluate.evaluate_hit, i.e., contain:
  {chrom, strand, score, query_start, query_end, query_len,
   cds: [{chrom, start, end}, ...]}
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import GENOMES, PROJECT_ROOT, V3_DATA  # noqa: E402

REPO = PROJECT_ROOT
V3 = V3_DATA
EVAL_HIT_FN_PATH = REPO / "scripts/v3_evaluate.py"

PER_PROTEIN_COLS = [
    "protein_id", "tool", "species", "regime", "padding",
    "failed", "wallclock_s",
    "exon_f1", "exon_precision", "exon_recall",
    "junction_f1", "junction_precision", "junction_recall",
    "cds_recovery", "query_coverage", "full_length",
    "n_pred_exons", "n_gt_exons", "exon_count_diff",
    "correct_chrom", "correct_locus",
]


_FASTA_CACHE: dict[str, pysam.FastaFile] = {}
_FASTA_CACHE_PID: int = 0


def get_fasta(species: str) -> pysam.FastaFile:
    """Return a per-process pysam.FastaFile handle.

    The cache is keyed by species but invalidated on PID change so
    forked workers get their own handle (pysam handles are not safe
    to share across processes).
    """
    global _FASTA_CACHE_PID
    pid = os.getpid()
    if _FASTA_CACHE_PID != pid:
        _FASTA_CACHE.clear()
        _FASTA_CACHE_PID = pid
    if species not in _FASTA_CACHE:
        meta = json.load(open(V3 / species / "meta.json"))
        _FASTA_CACHE[species] = pysam.FastaFile(meta["genome"])
    return _FASTA_CACHE[species]


def load_proteins(species: str, regime: str = "one2one") -> dict[str, str]:
    """Load protein sequences. regime ∈ {one2one, paralog}."""
    if regime == "one2one":
        fa = V3 / species / "proteins.fa"
    else:
        fa = V3 / f"{species}_paralog" / "proteins.fa"
    out: dict[str, str] = {}
    cur_id, cur_seq = "", []
    with fa.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id:
                    out[cur_id] = "".join(cur_seq)
                cur_id = line[1:].strip().split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_id:
            out[cur_id] = "".join(cur_seq)
    return out


def load_ground_truth(species: str, regime: str = "one2one") -> dict:
    if regime == "one2one":
        path = V3 / species / "ground_truth.json"
    else:
        path = V3 / f"{species}_paralog" / "ground_truth.json"
    return json.load(path.open())


def gt_region(gt_entry: dict, padding: int, chrom_lengths: dict[str, int]
              ) -> tuple[str, int, int, str, int]:
    """Return (chrom, start, end, strand, region_offset).
    Coordinates 0-based half-open. region_offset is the start.
    """
    exons = gt_entry["exons"]
    chrom = gt_entry["chrom"]
    strand = gt_entry.get("strand", "+")
    gt_start = min(e["start"] for e in exons)
    gt_end = max(e["end"] for e in exons)
    chrom_len = chrom_lengths.get(chrom, 1 << 31)
    rstart = max(0, gt_start - padding)
    rend = min(chrom_len, gt_end + padding)
    return chrom, rstart, rend, strand, rstart


def extract_region_fa(species: str, chrom: str, start: int, end: int,
                       out_path: Path) -> int:
    """Extract genomic region to a single-record FASTA. Returns region length.

    The header uses the original chrom name only — runners report this same
    name in their output and we add `start` (region_offset) to coordinates
    when computing genome-absolute positions.
    """
    fa = get_fasta(species)
    seq = fa.fetch(chrom, start, end)
    with out_path.open("w") as fh:
        fh.write(f">{chrom}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i+60] + "\n")
    return len(seq)


def chrom_lengths(species: str) -> dict[str, int]:
    fa = get_fasta(species)
    return {ref: fa.get_reference_length(ref) for ref in fa.references}


def write_protein_fa(protein_id: str, seq: str, out_path: Path) -> None:
    with out_path.open("w") as fh:
        fh.write(f">{protein_id}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i+60] + "\n")


def run_subproc(cmd: list[str], timeout: int = 300, cwd: str | None = None,
                env: dict | None = None) -> tuple[int, str, str]:
    """Run a subprocess with a wallclock cap. Returns (returncode, stdout, stderr)."""
    try:
        r = subprocess.run(cmd, capture_output=True, text=True,
                           timeout=timeout, cwd=cwd, env=env)
        return r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return -1, "", "TIMEOUT"
    except Exception as e:
        return -2, "", f"EXCEPTION: {e}"


def _import_evaluate_hit():
    """Lazy import to avoid circulating sys.path."""
    sys.path.insert(0, str(REPO / "scripts"))
    from v3_evaluate import evaluate_hit
    return evaluate_hit


def evaluate_hit_default(hit: dict | None, gt_entry: dict,
                          protein_len: int) -> dict:
    """Wrapper that returns a dict with `failed=True` if hit is None."""
    fn = _import_evaluate_hit()
    if hit is None:
        return {"failed": True, **{k: None for k in (
            "exon_f1", "exon_precision", "exon_recall",
            "junction_f1", "junction_precision", "junction_recall",
            "cds_recovery", "query_coverage", "full_length",
            "n_pred_exons", "n_gt_exons", "exon_count_diff",
            "correct_chrom", "correct_locus")}}
    metrics = fn(hit, gt_entry, protein_len)
    return {"failed": False, **metrics}


def offset_cds(cds: list[dict], chrom: str, region_offset: int) -> list[dict]:
    """Convert region-relative coords (1-based inclusive) to genome-absolute (0-based)."""
    out = []
    for c in cds:
        s = int(c["start"]) - 1 + region_offset
        e = int(c["end"]) + region_offset
        out.append({"chrom": chrom, "start": s, "end": e})
    return out
