"""Self-tests for Stage 0. Run with: python -m scripts.stage0.test_stage0"""

from __future__ import annotations

import os
import sys
import tempfile

from .common import (HSP, Candidate, interval_overlap, reciprocal_overlap,
                     reciprocal_overlap_frac, query_coverage, read_query_lengths)
from .aggregate import aggregate_hsps, dedup_candidates, top_k
from . import parsers
from .normalize import normalize
from .eval_locus import classify_query, load_gt_table, GTLocus


def test_overlap():
    assert interval_overlap(0, 10, 5, 15) == 5
    assert interval_overlap(0, 10, 10, 20) == 0
    # inter=50, a=b=100 → 0.5/0.5, passes at t=0.5 (inclusive), fails at t=0.51
    assert reciprocal_overlap(0, 100, 50, 150, 0.5) is True
    assert reciprocal_overlap(0, 100, 50, 150, 0.51) is False
    assert reciprocal_overlap(0, 100, 0, 1000, 0.5) is False  # 100/1000 = 0.1 < 0.5
    assert abs(reciprocal_overlap_frac(0, 100, 0, 1000) - 0.1) < 1e-9


def test_query_coverage():
    hsps = [
        HSP("c1", 0, 100, "+", 50.0, 10, 30),
        HSP("c1", 200, 300, "+", 40.0, 25, 60),  # overlapping query region
    ]
    qs, qe, frac = query_coverage(hsps, 100)
    assert qs == 10 and qe == 60
    # Union: [10,30) ∪ [25,60) = [10,60) → length 50
    assert abs(frac - 0.50) < 1e-9


def test_aggregate():
    # Same chrom, same strand, gap well under 200kb → merge
    hsps = [
        HSP("c1", 1000, 2000, "+", 80.0, 0, 30),
        HSP("c1", 5000, 6000, "+", 50.0, 30, 60),
    ]
    out = aggregate_hsps(hsps, query_len=100, gap=200_000)
    assert len(out) == 1
    c = out[0]
    assert c.start == 1000 and c.end == 6000
    assert abs(c.score - 130.0) < 1e-9
    assert abs(c.q_cov - 0.60) < 1e-9

    # Gap > threshold → split
    hsps2 = [
        HSP("c1", 1000, 2000, "+", 80.0, 0, 30),
        HSP("c1", 1_000_000, 1_001_000, "+", 50.0, 30, 60),
    ]
    out2 = aggregate_hsps(hsps2, query_len=100, gap=200_000)
    assert len(out2) == 2

    # Opposite strands never merge
    hsps3 = [
        HSP("c1", 1000, 2000, "+", 80.0, 0, 30),
        HSP("c1", 2500, 3000, "-", 50.0, 30, 60),
    ]
    out3 = aggregate_hsps(hsps3, query_len=100, gap=200_000)
    assert len(out3) == 2


def test_dedup_and_topk():
    cands = [
        Candidate("c1", 1000, 5000, "+", 100.0, 0, 50, 0.5),
        Candidate("c1", 1200, 4800, "+",  80.0, 0, 50, 0.5),  # near-dup of above
        Candidate("c2", 100, 200, "+",     40.0, 0, 10, 0.1),
    ]
    kept = dedup_candidates(cands, threshold=0.5)
    assert len(kept) == 2
    assert kept[0].score == 100.0
    top = top_k(kept, 1)
    assert len(top) == 1 and top[0].score == 100.0


def test_blast_tab_parser():
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        fh.write("q1\tchr1\t90.0\t50\t5\t0\t1\t50\t10000\t10150\t1e-20\t200.0\n")
        fh.write("q1\tchr1\t88.0\t40\t5\t0\t51\t90\t10300\t10420\t1e-15\t150.0\n")
        fh.write("q1\tchr2\t60.0\t30\t5\t0\t10\t40\t500\t590\t1e-05\t80.0\n")
        path = fh.name
    try:
        records = list(parsers.iter_hsps_blast_tab(path))
        assert len(records) == 3
        q, h = records[0]
        assert q == "q1" and h.chrom == "chr1" and h.start == 9999 and h.strand == "+"
        assert h.bitscore == 200.0 and h.q_start == 0 and h.q_end == 50
    finally:
        os.unlink(path)


def test_blat_psl_parser_skips_header():
    with tempfile.NamedTemporaryFile("w", suffix=".psl", delete=False) as fh:
        fh.write("psLayout version 3\n")
        fh.write("\n")
        fh.write("match\tmis-\trep.\tN's\tQ gap\tQ gap\tT gap\tT gap\tstrand\tQ        \tQ   \tQ    \tQ  \tT        \tT   \tT    \tT  \tblock\tblockSizes \tqStarts\t tStarts\n")
        fh.write("     \tmatch\tmatch\t   \tcount\tbases\tcount\tbases\t      \tname     \tsize\tstart\tend\tname     \tsize\tstart\tend\tcount\n")
        fh.write("-" * 100 + "\n")
        fh.write("143\t34\t0\t0\t4\t29\t6\t9864\t++\tENSP00000266743\t236\t30\t236\t10\t130530862\t88298422\t88308817\t7\t16,8,36,32,33,35,17,\t30,67,79,117,151,184,219,\t88298422,88299698,88302302,88303142,88305407,88308537,88308766,\n")
        path = fh.name
    try:
        records = list(parsers.iter_hsps_blat_psl(path))
        assert len(records) == 1, f"header rows should be skipped, got {len(records)}"
        q, h = records[0]
        assert q == "ENSP00000266743" and h.chrom == "10" and h.strand == "+"
        assert h.q_start == 30 and h.q_end == 236
        assert h.bitscore == 109.0  # 143 - 34
    finally:
        os.unlink(path)


def test_normalize_tblastn_end_to_end():
    with tempfile.TemporaryDirectory() as td:
        hits = os.path.join(td, "hits.tsv")
        prot = os.path.join(td, "p.fa")
        with open(hits, "w") as fh:
            fh.write("q1\tchr1\t90\t50\t5\t0\t1\t50\t10000\t10150\t1e-20\t200.0\n")
            fh.write("q1\tchr1\t88\t40\t5\t0\t51\t90\t10300\t10420\t1e-15\t150.0\n")
            fh.write("q1\tchr2\t60\t30\t5\t0\t10\t40\t500\t590\t1e-05\t80.0\n")
        with open(prot, "w") as fh:
            fh.write(">q1\n" + "M" * 100 + "\n")

        qlen = read_query_lengths(prot)
        out = list(normalize("tblastn", hits, qlen, top_k=10))
        # Two HSPs on chr1 within 200kb merge; chr2 stays separate.
        assert len(out) == 2
        scores = [c.score for _, c in out]
        assert scores[0] == 350.0   # 200+150 merged
        assert scores[1] == 80.0


def test_classify_query():
    cands = [
        Candidate("chr1", 1000, 5000, "+", 100.0, 0, 50, 0.5),
        Candidate("chr2", 100, 200, "+",   40.0, 0, 10, 0.1),
    ]
    gts = [GTLocus(chrom="chr1", start=2000, end=4000, strand="+", gene_id="G1")]
    v = classify_query("q1", cands, gts, threshold=0.5)
    # |cand∩GT|=2000, cand_len=4000, gt_len=2000 → 0.5/1.0 passes reciprocal
    assert v.search_ok and v.select_ok
    assert v.matched_gene_ids == ("G1",)


def main() -> int:
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"ok  {name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
