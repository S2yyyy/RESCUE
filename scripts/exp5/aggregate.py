"""Exp 5 — Block 2: aggregate per-protein TSVs into per-cell summary.

Reads results/exp5/per_protein/<tool>.<sp>.<regime>.<padding>.tsv
Writes results/exp5/cell/per_cell.tsv with columns:
  tool, species, regime, padding, n_total, n_failed, frac_failed,
  mean_exon_f1, ci_lo, ci_hi,
  mean_junction_f1, mean_cds_recovery, mean_full_length,
  mean_query_coverage,
  n_correct_locus, frac_correct_locus.

Bug-2 compliance: mean_* metrics are computed over n_total - n_failed
(skip-failures, no zero-fill). frac_failed is reported separately.
"""
from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import PROJECT_ROOT  # noqa: E402

REPO = PROJECT_ROOT
IN_DIR = REPO / "results/exp5/per_protein"
OUT = REPO / "results/exp5/cell/per_cell.tsv"


def bootstrap_ci(values, agg=np.mean, n=10000, seed=4242):
    if len(values) == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    arr = np.array(values, dtype=float)
    out = np.empty(n)
    n_v = len(arr)
    for i in range(n):
        idx = rng.integers(0, n_v, size=n_v)
        out[i] = agg(arr[idx])
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def aggregate_cell(rows: list[dict]) -> dict:
    n_total = len(rows)
    succ = [r for r in rows if r["failed"] != "True"]
    n_failed = n_total - len(succ)

    def col(name):
        out = []
        for r in succ:
            v = r.get(name)
            if v is None or v == "":
                continue
            try:
                out.append(float(v))
            except ValueError:
                continue
        return out

    f1 = col("exon_f1")
    jf1 = col("junction_f1")
    cdsr = col("cds_recovery")
    qcov = col("query_coverage")
    fl = col("full_length")  # 0/1
    cl = col("correct_locus")  # 0/1

    ci_lo, ci_hi = bootstrap_ci(f1)
    return {
        "n_total": n_total,
        "n_failed": n_failed,
        "frac_failed": n_failed / max(1, n_total),
        "mean_exon_f1": float(np.mean(f1)) if f1 else float("nan"),
        "ci_lo": ci_lo, "ci_hi": ci_hi,
        "mean_junction_f1": float(np.mean(jf1)) if jf1 else float("nan"),
        "mean_cds_recovery": float(np.mean(cdsr)) if cdsr else float("nan"),
        "mean_full_length": float(np.mean(fl)) if fl else float("nan"),
        "mean_query_coverage": float(np.mean(qcov)) if qcov else float("nan"),
        "n_correct_locus": int(sum(1 for v in cl if v >= 0.5)),
        "frac_correct_locus": (sum(1 for v in cl if v >= 0.5) / max(1, len(cl))) if cl else float("nan"),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", default=str(IN_DIR))
    ap.add_argument("--out", default=str(OUT))
    ap.add_argument("--bootstrap", type=int, default=10000)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="[%(asctime)s] B2 %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])
    log = logging.getLogger("aggregate")

    cells = []
    for tsv in sorted(Path(args.in_dir).glob("*.tsv")):
        name = tsv.stem  # tool.species.regime.padding
        parts = name.split(".")
        if len(parts) < 4:
            continue
        # tool can have underscore; let me parse from end:
        padding = int(parts[-1])
        regime = parts[-2]
        species = parts[-3]
        tool = ".".join(parts[:-3])
        rows = list(csv.DictReader(open(tsv), delimiter="\t"))
        if not rows:
            continue
        agg = aggregate_cell(rows)
        cells.append({
            "tool": tool, "species": species, "regime": regime,
            "padding": padding, **agg,
        })
        log.info("%-12s %-10s %s pad=%d  n=%d  failed=%d  exon_F1=%.4f  cl=%.4f",
                 tool, species, regime, padding,
                 agg["n_total"], agg["n_failed"],
                 agg["mean_exon_f1"], agg["frac_correct_locus"])

    cols = ["tool", "species", "regime", "padding",
            "n_total", "n_failed", "frac_failed",
            "mean_exon_f1", "ci_lo", "ci_hi",
            "mean_junction_f1", "mean_cds_recovery",
            "mean_full_length", "mean_query_coverage",
            "n_correct_locus", "frac_correct_locus"]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for c in cells:
            fh.write("\t".join(
                f"{c[k]:.4f}" if isinstance(c.get(k), float) else str(c.get(k, ""))
                for k in cols) + "\n")
    log.info("wrote %d cells to %s", len(cells), args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
