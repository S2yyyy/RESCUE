"""Stage 0 locus-level evaluator: Search-OK / Select-OK per query.

Pure post-processing. Takes a tool's finalized Candidate BED + a GT table
of protein → GT locus (chrom, start, end, strand) and reports per-query
A/B status for the search and selection stages.

GT = gene body span = (min(exon_start), max(exon_end)) on the annotated
strand. Alignment-level scoring (Align-OK / C class) is not handled here;
that's the job of a separate exon-level F1 evaluator on downstream CDS
predictions.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable

from .common import Candidate, reciprocal_overlap


# ---------------------------------------------------------------------------
# GT loading
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GTLocus:
    chrom: str
    start: int      # 0-based gene body start (min exon_start)
    end: int        # gene body end
    strand: str
    gene_id: str    # for broad/strict 1:many expansion


def load_gt_table(path: str) -> dict[str, list[GTLocus]]:
    """Load GT table. TSV with header:
        query_id  chrom  start  end  strand  gene_id

    One query may map to multiple GT loci (for 1:many broad/strict sets);
    the caller decides which subset to evaluate against.
    """
    gts: dict[str, list[GTLocus]] = defaultdict(list)
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"query_id", "chrom", "start", "end", "strand", "gene_id"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"GT table missing columns: {sorted(missing)}")
        for row in reader:
            gts[row["query_id"]].append(GTLocus(
                chrom=row["chrom"],
                start=int(row["start"]),
                end=int(row["end"]),
                strand=row["strand"],
                gene_id=row["gene_id"],
            ))
    return gts


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------

def _locus_match(cand: Candidate, gt: GTLocus, threshold: float,
                 enforce_strand: bool) -> bool:
    if cand.chrom != gt.chrom:
        return False
    if enforce_strand and cand.strand != gt.strand:
        return False
    return reciprocal_overlap(cand.start, cand.end, gt.start, gt.end, threshold)


def any_match(cands: Iterable[Candidate], gts: Iterable[GTLocus],
              threshold: float = 0.5, enforce_strand: bool = True) -> bool:
    """True if at least one candidate reciprocally overlaps at least one GT locus."""
    gts = list(gts)
    for c in cands:
        for g in gts:
            if _locus_match(c, g, threshold, enforce_strand):
                return True
    return False


# ---------------------------------------------------------------------------
# Per-query classification
# ---------------------------------------------------------------------------

@dataclass
class LocusVerdict:
    query_id: str
    search_ok: bool
    select_ok: bool
    top1_chrom: str | None
    top1_start: int | None
    top1_end: int | None
    matched_gene_ids: tuple[str, ...]  # among candidates


def classify_query(query_id: str,
                   cands: list[Candidate],
                   gts: list[GTLocus],
                   threshold: float = 0.5,
                   enforce_strand: bool = True) -> LocusVerdict:
    """cands must already be finalized (dedup + top-K) and score-sorted."""
    if not cands:
        return LocusVerdict(query_id, False, False, None, None, None, ())

    matched = []
    for c in cands:
        for g in gts:
            if _locus_match(c, g, threshold, enforce_strand):
                matched.append(g.gene_id)

    search_ok = len(matched) > 0
    top1 = cands[0]
    select_ok = search_ok and any(
        _locus_match(top1, g, threshold, enforce_strand) for g in gts
    )
    return LocusVerdict(
        query_id=query_id,
        search_ok=search_ok,
        select_ok=select_ok,
        top1_chrom=top1.chrom,
        top1_start=top1.start,
        top1_end=top1.end,
        matched_gene_ids=tuple(sorted(set(matched))),
    )


# ---------------------------------------------------------------------------
# Batch evaluation over a normalized BED
# ---------------------------------------------------------------------------

def load_candidates_by_query(bed_path: str) -> dict[str, list[Candidate]]:
    """Read a Stage-0 BED where column 4 (name) = query_id.

    Preserves input order within each query, so the file producer is
    responsible for sorting by score desc before writing.
    """
    per_query: dict[str, list[Candidate]] = defaultdict(list)
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            qid = f[3]
            per_query[qid].append(Candidate.from_bed_row(line))
    return per_query


def evaluate_bed(bed_path: str, gt_path: str,
                 threshold: float = 0.5,
                 enforce_strand: bool = True,
                 out_tsv: str | None = None) -> list[LocusVerdict]:
    cands_by_query = load_candidates_by_query(bed_path)
    gts_by_query = load_gt_table(gt_path)

    verdicts: list[LocusVerdict] = []
    all_queries = sorted(set(cands_by_query) | set(gts_by_query))
    for qid in all_queries:
        v = classify_query(
            qid,
            cands_by_query.get(qid, []),
            gts_by_query.get(qid, []),
            threshold=threshold,
            enforce_strand=enforce_strand,
        )
        verdicts.append(v)

    if out_tsv:
        with open(out_tsv, "w") as out:
            out.write("query_id\tsearch_ok\tselect_ok\ttop1_chrom\ttop1_start\ttop1_end\tmatched_gene_ids\n")
            for v in verdicts:
                out.write("\t".join([
                    v.query_id,
                    str(int(v.search_ok)),
                    str(int(v.select_ok)),
                    v.top1_chrom or "",
                    "" if v.top1_start is None else str(v.top1_start),
                    "" if v.top1_end is None else str(v.top1_end),
                    ",".join(v.matched_gene_ids),
                ]) + "\n")
    return verdicts


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _main() -> int:
    p = argparse.ArgumentParser(description="Stage 0 locus evaluator")
    p.add_argument("--bed", required=True, help="Finalized Stage-0 BED")
    p.add_argument("--gt", required=True, help="GT table TSV")
    p.add_argument("--threshold", type=float, default=0.5,
                   help="Reciprocal overlap threshold (default 0.5)")
    p.add_argument("--ignore-strand", action="store_true",
                   help="Disable strand enforcement for locus match")
    p.add_argument("-o", "--out", help="Output TSV")
    args = p.parse_args()

    verdicts = evaluate_bed(
        args.bed, args.gt,
        threshold=args.threshold,
        enforce_strand=not args.ignore_strand,
        out_tsv=args.out,
    )
    n = len(verdicts)
    s = sum(v.search_ok for v in verdicts)
    sel = sum(v.select_ok for v in verdicts)
    sys.stderr.write(f"queries={n}  search_ok={s} ({s/n:.2%})  select_ok={sel} ({sel/n:.2%})\n")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
