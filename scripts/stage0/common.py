"""Stage 0 shared types and geometry.

Design contract (see EXPERIMENT_PLAN_V4_PER_TOOL_ABC.md §Stage 0):
  - A Candidate is an 8-tuple locus-level record emitted by any tool after
    Stage 0 normalization. Fields are tool-neutral; native scores are kept
    per-tool and are never rescaled across tools.
  - A HSP is a single local alignment used as raw material for HSP-only
    tools before aggregation.
  - Overlap math is fixed here; callers must not redefine it.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Iterable, Iterator


# ---------------------------------------------------------------------------
# Record types
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class HSP:
    """One raw local hit, pre-aggregation."""
    chrom: str
    start: int          # 0-based, half-open
    end: int
    strand: str         # '+' / '-'
    bitscore: float
    q_start: int        # 0-based on the query protein
    q_end: int

    def __post_init__(self):
        assert self.strand in ("+", "-"), f"bad strand {self.strand!r}"
        assert self.start < self.end, f"bad interval {self.start}-{self.end}"
        assert self.q_start < self.q_end, \
            f"bad query interval {self.q_start}-{self.q_end}"


@dataclass(frozen=True)
class Candidate:
    """Locus-level record after Stage 0 normalization. Top-K operates on this."""
    chrom: str
    start: int
    end: int
    strand: str
    score: float        # tool-native, higher-is-better convention (see normalize.py)
    q_start: int
    q_end: int
    q_cov: float        # (q_end - q_start) / query_length, in [0, 1]

    def to_bed_row(self, name: str = ".") -> str:
        """BED12-ish 8-col row: chrom start end name score strand q_start q_end q_cov."""
        return "\t".join([
            self.chrom,
            str(self.start),
            str(self.end),
            name,
            f"{self.score:.6g}",
            self.strand,
            str(self.q_start),
            str(self.q_end),
            f"{self.q_cov:.4f}",
        ])

    @classmethod
    def from_bed_row(cls, line: str) -> "Candidate":
        f = line.rstrip("\n").split("\t")
        return cls(
            chrom=f[0],
            start=int(f[1]),
            end=int(f[2]),
            strand=f[5],
            score=float(f[4]),
            q_start=int(f[6]),
            q_end=int(f[7]),
            q_cov=float(f[8]),
        )


BED_HEADER = "#chrom\tstart\tend\tname\tscore\tstrand\tq_start\tq_end\tq_cov"


# ---------------------------------------------------------------------------
# Interval geometry (single source of truth)
# ---------------------------------------------------------------------------

def interval_overlap(a_start: int, a_end: int,
                     b_start: int, b_end: int) -> int:
    """Length of intersection, 0 if disjoint."""
    lo = max(a_start, b_start)
    hi = min(a_end, b_end)
    return max(0, hi - lo)


def reciprocal_overlap(a_start: int, a_end: int,
                       b_start: int, b_end: int,
                       threshold: float = 0.5) -> bool:
    """Stage 0 correct-locus criterion: |A∩B|/|A| ≥ t AND |A∩B|/|B| ≥ t."""
    inter = interval_overlap(a_start, a_end, b_start, b_end)
    if inter == 0:
        return False
    a_len = a_end - a_start
    b_len = b_end - b_start
    if a_len <= 0 or b_len <= 0:
        return False
    return (inter / a_len >= threshold) and (inter / b_len >= threshold)


def reciprocal_overlap_frac(a_start: int, a_end: int,
                            b_start: int, b_end: int) -> float:
    """min(|A∩B|/|A|, |A∩B|/|B|). Useful for continuous sensitivity sweeps."""
    inter = interval_overlap(a_start, a_end, b_start, b_end)
    if inter == 0:
        return 0.0
    a_len = a_end - a_start
    b_len = b_end - b_start
    if a_len <= 0 or b_len <= 0:
        return 0.0
    return min(inter / a_len, inter / b_len)


# ---------------------------------------------------------------------------
# Query-length lookup
# ---------------------------------------------------------------------------

def read_query_lengths(path: str) -> dict[str, int]:
    """Read protein FASTA (or 2-col TSV `id<TAB>length`) → {id: length}.

    FASTA: id is the first whitespace-separated token after '>'.
    """
    lengths: dict[str, int] = {}
    with open(path) as fh:
        first = fh.readline()
        fh.seek(0)
        if first.startswith(">"):
            cur_id = None
            cur_len = 0
            for line in fh:
                line = line.rstrip()
                if line.startswith(">"):
                    if cur_id is not None:
                        lengths[cur_id] = cur_len
                    cur_id = line[1:].split()[0]
                    cur_len = 0
                else:
                    cur_len += len(line)
            if cur_id is not None:
                lengths[cur_id] = cur_len
        else:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                qid, qlen = line.rstrip().split("\t")[:2]
                lengths[qid] = int(qlen)
    return lengths


# ---------------------------------------------------------------------------
# Query-coverage computation for a set of HSPs on the same query
# ---------------------------------------------------------------------------

def query_coverage(hsps: Iterable[HSP], query_len: int) -> tuple[int, int, float]:
    """Return (q_start_min, q_end_max, covered_fraction).

    q_cov uses the UNION of HSP query ranges divided by query_len, so
    overlapping HSPs don't double-count. q_start / q_end record the envelope.
    """
    intervals = sorted((h.q_start, h.q_end) for h in hsps)
    if not intervals:
        return 0, 0, 0.0
    merged: list[tuple[int, int]] = []
    for s, e in intervals:
        if merged and s <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    covered = sum(e - s for s, e in merged)
    qs = merged[0][0]
    qe = merged[-1][1]
    frac = covered / query_len if query_len > 0 else 0.0
    return qs, qe, frac
