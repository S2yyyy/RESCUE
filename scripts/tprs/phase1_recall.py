"""TPRS Phase 1b — Normalize relaxed-search outputs to BED and evaluate recall.

Reads raw outputs from results/tprs/raw/, normalizes via the existing Stage 0
pipeline, then computes recall@K against paralog_broad GT.

Outputs:
  results/tprs/stage0_bed/<tool>.<species>.paralog.bed
  results/tprs/recall/phase1_recall.tsv
"""
from __future__ import annotations

import argparse
import logging
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import PROJECT_ROOT, V3_DATA  # noqa: E402

REPO = PROJECT_ROOT
RAW = REPO / "results/tprs/raw"
BED_DIR = REPO / "results/tprs/stage0_bed"
RECALL_DIR = REPO / "results/tprs/recall"
DATA = V3_DATA
GT_DIR = REPO / "results/exp2/gt"

SPECIES = ["mouse", "dog", "chicken", "zebrafish"]

# Map relaxed tool names → (normalize_tool, file_pattern)
# All paralog-only.
TOOL_SPEC = [
    ("miniprot_ultra",  "miniprot", "{raw}/miniprot_ultra.{sp}.paralog.gff"),
    ("diamond_relaxed", "diamond",  "{raw}/diamond_relaxed.{sp}.paralog.tsv"),
    ("mmseqs_relaxed",  "mmseqs",   "{raw}/mmseqs_relaxed.{sp}.paralog.tsv"),
    ("last_relaxed",    "last-swap","{raw}/last_relaxed.{sp}.paralog.tab"),
]

# Baseline tools from Exp3 for comparison
BASELINE_TOOLS = ["miniprot", "diamond", "mmseqs", "last", "genblasta",
                  "spaln3", "tblastn", "blat", "genblastg"]

RNG = np.random.default_rng(4242)


def proteins_fa(sp: str) -> Path:
    return DATA / f"{sp}_paralog" / "proteins.fa"


def normalize_one(paper_tool: str, normalize_tool: str,
                  native_input: Path, proteins: Path,
                  top_k: int, out_bed: Path) -> dict:
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    if not native_input.exists() or native_input.stat().st_size == 0:
        return {"status": "MISSING_INPUT", "n_records": 0, "wall_sec": 0,
                "input": str(native_input)}
    tmp = out_bed.with_suffix(out_bed.suffix + ".tmp")
    cmd = [
        sys.executable, "-m", "scripts.stage0.normalize",
        "--tool", normalize_tool,
        "--input", str(native_input),
        "--query-lengths", str(proteins),
        "--top-k", str(top_k),
        "-o", str(tmp),
    ]
    t = time.time()
    proc = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    wall = time.time() - t
    if proc.returncode != 0:
        err = proc.stderr.strip().splitlines()[-3:] if proc.stderr else []
        return {"status": "FAILED", "n_records": 0, "wall_sec": wall,
                "input": str(native_input), "stderr": err}
    n_records = sum(1 for _ in tmp.open()) - 1  # minus header
    tmp.replace(out_bed)
    return {"status": "DONE", "n_records": n_records, "wall_sec": wall,
            "input": str(native_input)}


def load_bed(path: Path) -> dict[str, list[tuple]]:
    """Return {qid: [(chrom, start, end, strand), ...]} preserving file order."""
    out: dict[str, list[tuple]] = defaultdict(list)
    if not path.exists():
        return out
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            out[f[3]].append((f[0], int(f[1]), int(f[2]), f[5]))
    return out


def load_gt(gt_path: Path) -> dict[str, list[tuple]]:
    out: dict[str, list[tuple]] = defaultdict(list)
    with gt_path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[idx["query_id"]]].append((
                f[idx["chrom"]],
                int(f[idx["start"]]),
                int(f[idx["end"]]),
                f[idx["strand"]],
            ))
    return out


def reciprocal_overlap_50(a, b) -> bool:
    if a[0] != b[0] or a[3] != b[3]:
        return False
    inter = max(0, min(a[2], b[2]) - max(a[1], b[1]))
    if inter == 0:
        return False
    al = a[2] - a[1]
    bl = b[2] - b[1]
    return inter / al >= 0.5 and inter / bl >= 0.5


def query_hit(cands: list[tuple], gts: list[tuple], k: int) -> bool:
    for c in cands[:k]:
        for g in gts:
            if reciprocal_overlap_50(c, g):
                return True
    return False


def bootstrap_ci(hits: np.ndarray, n_resamples: int = 10000) -> tuple[float, float]:
    if len(hits) == 0:
        return float("nan"), float("nan")
    n = len(hits)
    idx = RNG.integers(0, n, size=(n_resamples, n))
    means = hits[idx].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def union_query_hit(qid: str, all_tool_cands: dict[str, dict[str, list[tuple]]],
                    gts: dict[str, list[tuple]], k: int) -> bool:
    """Any candidate from any tool overlaps any GT locus (union of top-K per tool)."""
    union_cands: set[tuple] = set()
    for tool_cands in all_tool_cands.values():
        for c in tool_cands.get(qid, [])[:k]:
            union_cands.add(c)
    for c in union_cands:
        for g in gts.get(qid, []):
            if reciprocal_overlap_50(c, g):
                return True
    return False


def main() -> int:
    ap = argparse.ArgumentParser(description="TPRS Phase 1b: normalize + recall")
    ap.add_argument("--top-k", type=int, default=200)
    ap.add_argument("--skip-normalize", action="store_true",
                    help="Skip normalization, only compute recall")
    ap.add_argument("--tools", nargs="*", default=None)
    args = ap.parse_args()

    BED_DIR.mkdir(parents=True, exist_ok=True)
    RECALL_DIR.mkdir(parents=True, exist_ok=True)

    log_path = RECALL_DIR / "phase1.log"
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] TPRS %(message)s",
        handlers=[logging.StreamHandler(sys.stdout),
                  logging.FileHandler(log_path, mode="a")],
    )
    log = logging.getLogger("phase1")

    # Step 1: Normalize raw outputs → BED
    if not args.skip_normalize:
        for paper_tool, norm_tool, tmpl in TOOL_SPEC:
            if args.tools and paper_tool not in args.tools:
                continue
            for sp in SPECIES:
                native_input = Path(tmpl.format(raw=str(RAW), sp=sp))
                out_bed = BED_DIR / f"{paper_tool}.{sp}.paralog.bed"
                if out_bed.exists() and out_bed.stat().st_size > 0:
                    log.info("skip normalize %s %s (exists)", paper_tool, sp)
                    continue
                log.info("begin normalize %s %s", paper_tool, sp)
                result = normalize_one(paper_tool, norm_tool, native_input,
                                       proteins_fa(sp), args.top_k, out_bed)
                log.info("done normalize %s %s status=%s n=%d wall=%.1fs",
                         paper_tool, sp, result["status"],
                         result["n_records"], result["wall_sec"])
                if result["status"] == "FAILED":
                    log.error("stderr: %s", result.get("stderr", []))

    # Step 2: Compute recall
    K_VALUES = [1, 10, 50, 200]
    rows = []

    # Collect all BED files (relaxed + baseline for comparison)
    all_tool_names = [t[0] for t in TOOL_SPEC]
    if not args.tools:
        all_tool_names += BASELINE_TOOLS
    else:
        all_tool_names += [t for t in BASELINE_TOOLS if t in args.tools]

    for sp in SPECIES:
        gt = GT_DIR / f"{sp}.paralog_broad.tsv"
        gts = load_gt(gt)
        qids = sorted(set(gts))
        n_total = len(qids)
        if n_total == 0:
            continue

        # Load all tool candidates
        all_cands: dict[str, dict[str, list[tuple]]] = {}
        for tname in all_tool_names:
            bed_path = BED_DIR / f"{tname}.{sp}.paralog.bed"
            # For baseline tools, use Exp3 BED
            if tname in BASELINE_TOOLS:
                bed_path = REPO / "results/exp3/stage0_bed" / f"{tname}.{sp}.paralog.bed"
            all_cands[tname] = load_bed(bed_path)

        for tname in all_tool_names:
            cands = all_cands[tname]
            for k in K_VALUES:
                hits = np.array([
                    int(query_hit(cands.get(q, []), gts[q], k)) for q in qids
                ])
                recall = float(hits.mean())
                lo, hi = bootstrap_ci(hits)
                rows.append({
                    "tool": tname, "species": sp, "k": k,
                    "n_total": n_total, "n_hit": int(hits.sum()),
                    "recall": recall, "ci_lo": lo, "ci_hi": hi,
                })

        # Union recall: all relaxed tools combined
        relaxed_tool_names = [t[0] for t in TOOL_SPEC]
        for k in K_VALUES:
            hits = np.array([
                int(union_query_hit(q, {t: all_cands[t] for t in relaxed_tool_names
                                       if t in all_cands}, gts, k))
                for q in qids
            ])
            recall = float(hits.mean())
            lo, hi = bootstrap_ci(hits)
            rows.append({
                "tool": "union_relaxed", "species": sp, "k": k,
                "n_total": n_total, "n_hit": int(hits.sum()),
                "recall": recall, "ci_lo": lo, "ci_hi": hi,
            })
            log.info("union_relaxed %s recall@%d=%.3f [%.3f, %.3f]  hit=%d/%d",
                     sp, k, recall, lo, hi, int(hits.sum()), n_total)

    # Write results
    out_tsv = RECALL_DIR / "phase1_recall.tsv"
    cols = ["tool", "species", "k", "n_total", "n_hit", "recall", "ci_lo", "ci_hi"]
    with out_tsv.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    log.info("wrote %d rows to %s", len(rows), out_tsv)

    # Summary table: recall@50 for relaxed tools vs baselines
    print("\n=== TPRS Phase 1 Recall@50 Summary ===")
    print(f"{'Tool':20s} {'Mouse':>8s} {'Dog':>8s} {'Chicken':>8s} {'Zebrafish':>8s}")
    print("-" * 56)
    for tname in all_tool_names + ["union_relaxed"]:
        vals = []
        for sp in SPECIES:
            matching = [r for r in rows if r["tool"] == tname
                        and r["species"] == sp and r["k"] == 50]
            if matching:
                r = matching[0]
                vals.append(f"{r['recall']:.3f}")
            else:
                vals.append("   N/A")
        print(f"{tname:20s} {' '.join(f'{v:>8s}' for v in vals)}")

    # Key finding: delta of relaxed vs baseline
    print("\n=== Recall@50 Delta (relaxed - baseline) ===")
    for rt, bt in [("miniprot_ultra", "miniprot"), ("diamond_relaxed", "diamond"),
                   ("mmseqs_relaxed", "mmseqs"), ("last_relaxed", "last")]:
        deltas = []
        for sp in SPECIES:
            r_rows = [r for r in rows if r["tool"] == rt and r["species"] == sp and r["k"] == 50]
            b_rows = [r for r in rows if r["tool"] == bt and r["species"] == sp and r["k"] == 50]
            if r_rows and b_rows:
                d = r_rows[0]["recall"] - b_rows[0]["recall"]
                deltas.append(f"{d:+.3f}")
            else:
                deltas.append("   N/A")
        print(f"{rt:20s} {' '.join(f'{d:>8s}' for d in deltas)}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
