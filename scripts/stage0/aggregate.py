"""Stage 0 HSP aggregation, dedup, and top-K selection.

These are post-processing operations that ONLY act on normalized records.
No tool's native output is modified; existing parsers feed this module.

Rules (EXPERIMENT_PLAN_V4_PER_TOOL_ABC.md §Stage 0):
  - HSP → region: same-strand, same-chrom, gap ≤ 200 kb are merged.
    region score = sum of member bit scores (BLAST-family convention).
    region coords = (min start, max end); q_start/q_end = union envelope;
    q_cov = union / query_length.
  - Dedup: among candidates from the SAME tool, if two reciprocally overlap
    by ≥ 50% on the genome, keep the higher-scoring one. No cross-tool dedup.
  - Top-K: after dedup, sort by (-score, chrom_lex, start, end), keep K.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable

from .common import HSP, Candidate, interval_overlap, query_coverage, reciprocal_overlap


# ---------------------------------------------------------------------------
# HSP → Candidate aggregation (HSP-only tools)
# ---------------------------------------------------------------------------

def aggregate_hsps(hsps: Iterable[HSP],
                   query_len: int,
                   gap: int = 200_000,
                   combiner: str = "sum") -> list[Candidate]:
    """Merge HSPs on the same chrom+strand with adjacent gaps ≤ `gap`.

    `gap` is the distance between consecutive HSP endpoints (next.start − prev.end).
    Overlapping HSPs (gap < 0) are always merged.

    `combiner` controls the per-cluster score:
      - "sum" (default): score = Σ member bit scores (BLAST-family convention).
      - "max": score = max member bit score (used by Exp 4 B4 sensitivity).
    """
    by_key: dict[tuple[str, str], list[HSP]] = defaultdict(list)
    for h in hsps:
        by_key[(h.chrom, h.strand)].append(h)

    out: list[Candidate] = []
    for (chrom, strand), group in by_key.items():
        group.sort(key=lambda h: (h.start, h.end))
        cluster: list[HSP] = []
        cluster_end = -1
        for h in group:
            if cluster and h.start - cluster_end > gap:
                out.append(_flush_cluster(cluster, chrom, strand, query_len, combiner))
                cluster = []
                cluster_end = -1
            cluster.append(h)
            cluster_end = max(cluster_end, h.end)
        if cluster:
            out.append(_flush_cluster(cluster, chrom, strand, query_len, combiner))
    return out


def _flush_cluster(cluster: list[HSP], chrom: str, strand: str,
                   query_len: int, combiner: str = "sum") -> Candidate:
    start = min(h.start for h in cluster)
    end = max(h.end for h in cluster)
    if combiner == "max":
        score = max(h.bitscore for h in cluster)
    else:
        score = sum(h.bitscore for h in cluster)
    qs, qe, qcov = query_coverage(cluster, query_len)
    return Candidate(chrom=chrom, start=start, end=end, strand=strand,
                     score=score, q_start=qs, q_end=qe, q_cov=qcov)


# ---------------------------------------------------------------------------
# Dedup: within-tool, genome-side reciprocal 50% overlap
# ---------------------------------------------------------------------------

def dedup_candidates(cands: Iterable[Candidate],
                     threshold: float = 0.5) -> list[Candidate]:
    """Remove near-duplicate candidates from the same tool.

    Greedy by score descending: a candidate is dropped if it reciprocally
    overlaps a previously kept candidate on the same chrom+strand.
    """
    sorted_cands = sorted(cands, key=lambda c: (-c.score, c.chrom, c.start, c.end))
    kept: list[Candidate] = []
    kept_by_key: dict[tuple[str, str], list[Candidate]] = defaultdict(list)
    for c in sorted_cands:
        key = (c.chrom, c.strand)
        dup = False
        for k in kept_by_key[key]:
            if reciprocal_overlap(c.start, c.end, k.start, k.end, threshold):
                dup = True
                break
        if not dup:
            kept.append(c)
            kept_by_key[key].append(c)
    return kept


# ---------------------------------------------------------------------------
# Top-K (global across the whole assembly, after dedup)
# ---------------------------------------------------------------------------

def top_k(cands: Iterable[Candidate], k: int) -> list[Candidate]:
    """Stable sort by (-score, chrom, start, end); truncate at k."""
    ordered = sorted(cands, key=lambda c: (-c.score, c.chrom, c.start, c.end))
    return ordered[:k] if k > 0 else ordered


# ---------------------------------------------------------------------------
# Pipeline convenience
# ---------------------------------------------------------------------------

def finalize_candidates(cands: Iterable[Candidate],
                        k: int = 50,
                        dedup_threshold: float = 0.5) -> list[Candidate]:
    """Dedup then top-K. Used by both HSP and native-output paths."""
    return top_k(dedup_candidates(cands, dedup_threshold), k)
