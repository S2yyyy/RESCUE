"""Stage 0 driver.

CLI that takes a tool's native output file and writes a finalized
Stage-0 BED (dedup + top-K applied). No tool binaries are invoked here;
this operates only on files already produced by a separate run step.

Usage
-----
    python -m scripts.stage0.normalize \
        --tool tblastn \
        --input path/to/hits.tsv \
        --query-lengths path/to/proteins.fa \
        --top-k 50 \
        -o out.bed
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from typing import Iterator

from .common import BED_HEADER, Candidate, HSP, read_query_lengths
from .aggregate import aggregate_hsps, finalize_candidates
from . import parsers


# ---------------------------------------------------------------------------
# Tool registry
# ---------------------------------------------------------------------------

# HSP-only tools: raw HSPs feed aggregate_hsps per query before finalize.
HSP_TOOLS = {
    "tblastn":  parsers.iter_hsps_blast_tab,
    "tblastn-report": parsers.iter_hsps_blast_report,   # genBlastA `.blast.report`
    "diamond":  parsers.iter_hsps_blast_tab_swapped,    # DIAMOND blastx (genome=query, proteins=db)
    "mmseqs":   parsers.iter_hsps_blast_tab,            # MMseqs2 easy-search (protein=query)
    "last":     parsers.iter_hsps_last_tab,
    "last-swap": parsers.iter_hsps_blast_tab_swapped,   # LAST DNA-as-query (matches DIAMOND orientation)
    "blat":     parsers.iter_hsps_blat_psl,
}

# Native locus-level tools: Candidates are emitted directly; only dedup + top-K apply.
NATIVE_TOOLS = {
    "miniprot":             parsers.iter_candidates_miniprot_paf,
    "spaln":                parsers.iter_candidates_spaln_gff,
    "gth":                  parsers.iter_candidates_gth_gff,
    "scipio":               parsers.iter_candidates_scipio_yaml,
    "procompart":           parsers.iter_candidates_procompart,
    "genblasta":            parsers.iter_candidates_genblasta,
    "genblastg":            parsers.iter_candidates_genblastg,
    # miniprot variants reuse the miniprot PAF parser
    "minisplice-miniprot":  parsers.iter_candidates_miniprot_paf,
    "syncmer-miniprot":     parsers.iter_candidates_miniprot_paf,
}

ALL_TOOLS = set(HSP_TOOLS) | set(NATIVE_TOOLS)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def normalize(tool: str, input_path: str, query_lengths: dict[str, int],
              top_k: int = 50, agg_gap: int = 200_000,
              dedup_threshold: float = 0.5,
              hsp_combiner: str = "sum",
              ) -> Iterator[tuple[str, Candidate]]:
    """Yield (query_id, Candidate) after Stage 0 normalization."""
    tool = tool.lower()
    if tool in HSP_TOOLS:
        per_query: dict[str, list[HSP]] = defaultdict(list)
        for qid, hsp in HSP_TOOLS[tool](input_path):
            per_query[qid].append(hsp)
        for qid, hsps in per_query.items():
            qlen = query_lengths.get(qid, 0)
            cands = aggregate_hsps(hsps, query_len=qlen, gap=agg_gap,
                                    combiner=hsp_combiner)
            for c in finalize_candidates(cands, k=top_k,
                                         dedup_threshold=dedup_threshold):
                yield qid, c
    elif tool in NATIVE_TOOLS:
        per_query: dict[str, list[Candidate]] = defaultdict(list)
        for qid, cand in NATIVE_TOOLS[tool](input_path, query_lengths):
            per_query[qid].append(cand)
        for qid, cands in per_query.items():
            for c in finalize_candidates(cands, k=top_k,
                                         dedup_threshold=dedup_threshold):
                yield qid, c
    else:
        raise ValueError(f"Unknown tool {tool!r}. Known: {sorted(ALL_TOOLS)}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _main() -> int:
    p = argparse.ArgumentParser(description="Stage 0 normalizer")
    p.add_argument("--tool", required=True, choices=sorted(ALL_TOOLS))
    p.add_argument("--input", required=True, help="Tool's native output file")
    p.add_argument("--query-lengths", required=True,
                   help="Protein FASTA or 2-col TSV with query lengths")
    p.add_argument("--top-k", type=int, default=50)
    p.add_argument("--agg-gap", type=int, default=200_000,
                   help="HSP aggregation gap (bp); ignored for native tools")
    p.add_argument("--dedup-threshold", type=float, default=0.5,
                   help="Reciprocal overlap threshold for dedup")
    p.add_argument("--hsp-combiner", choices=["sum", "max"], default="sum",
                   help="HSP-cluster score combiner (sum=BLAST default, max=B4 sensitivity)")
    p.add_argument("-o", "--out", required=True, help="Output BED path")
    args = p.parse_args()

    qlen = read_query_lengths(args.query_lengths)
    n_records = 0
    with open(args.out, "w") as fh:
        fh.write(BED_HEADER + "\n")
        for qid, cand in normalize(
                tool=args.tool,
                input_path=args.input,
                query_lengths=qlen,
                top_k=args.top_k,
                agg_gap=args.agg_gap,
                dedup_threshold=args.dedup_threshold,
                hsp_combiner=args.hsp_combiner):
            fh.write(cand.to_bed_row(name=qid) + "\n")
            n_records += 1
    sys.stderr.write(f"wrote {n_records} records to {args.out}\n")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
