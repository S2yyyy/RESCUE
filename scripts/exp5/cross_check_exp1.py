"""B4 — paired exon F1 between Exp 5 forced-locus and Exp 1 natural-context.

For each (tool, species) shared between the two experiments, restrict to
proteins where Exp 1 reports correct_locus=1 (so the natural-context query
landed on the correct locus and was scored on it). Compute the paired
mean exon F1 difference: Δ = mean(exp5) - mean(exp1).

Acceptance: |Δ| ≤ 0.02 on ≥90% of cells.

Inputs:
  results/exp5/per_protein/<tool>.<sp>.one2one.50000.tsv
  results/exp1/per_protein.tsv (with `tool, species, regime, exon_f1, correct_locus`)

Outputs:
  results/exp5/cross/exp1_consistency.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import PROJECT_ROOT  # noqa: E402

REPO = PROJECT_ROOT
EXP5_DIR = REPO / "results/exp5/per_protein"
EXP1_TSV = REPO / "results/exp1/per_protein.tsv"
OUT = REPO / "results/exp5/cross/exp1_consistency.tsv"


SHARED_TOOLS = ["miniprot", "spaln3", "gth", "exonerate", "genewise",
                "scipio", "genblastg", "minisplice_miniprot"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--padding", type=int, default=50000)
    ap.add_argument("--bootstrap", type=int, default=10000)
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args()

    # Index Exp 1: {(tool, species, protein_id) -> (exon_f1, correct_locus)}
    e1: dict[tuple, tuple[float, int]] = {}
    e1_tool_alias = {"spaln": "spaln3"}
    with open(EXP1_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("regime") not in ("one2one", "one_to_one", "1to1"):
                continue
            tool = e1_tool_alias.get(r["tool"], r["tool"])
            try:
                f1 = float(r["exon_f1"])
                cl = int(float(r["correct_locus"]))
            except (ValueError, KeyError):
                continue
            e1[(tool, r["species"], r["protein_id"])] = (f1, cl)

    rng = np.random.default_rng(4242)
    rows = []
    for tool in SHARED_TOOLS:
        for sp in ("mouse", "dog", "chicken", "zebrafish"):
            tsv = EXP5_DIR / f"{tool}.{sp}.one2one.{args.padding}.tsv"
            if not tsv.exists():
                continue
            paired_e5: list[float] = []
            paired_e1: list[float] = []
            with tsv.open() as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r.get("failed") == "True":
                        continue
                    pid = r["protein_id"]
                    e1_rec = e1.get((tool, sp, pid))
                    if e1_rec is None:
                        continue
                    e1_f1, cl = e1_rec
                    if cl != 1:
                        continue
                    try:
                        paired_e5.append(float(r["exon_f1"]))
                        paired_e1.append(e1_f1)
                    except ValueError:
                        continue
            if not paired_e5:
                continue
            a = np.array(paired_e5, dtype=float)
            b = np.array(paired_e1, dtype=float)
            diff = a - b
            n = len(diff)
            deltas = np.empty(args.bootstrap)
            for i in range(args.bootstrap):
                idx = rng.integers(0, n, size=n)
                deltas[i] = diff[idx].mean()
            ci_lo = float(np.percentile(deltas, 2.5))
            ci_hi = float(np.percentile(deltas, 97.5))
            mean_d = float(diff.mean())
            verdict = "pass" if abs(mean_d) <= 0.02 else (
                "warn" if abs(mean_d) <= 0.05 else "fail")
            rows.append({
                "tool": tool, "species": sp, "n_paired": n,
                "mean_exon_f1_exp5": float(a.mean()),
                "mean_exon_f1_exp1": float(b.mean()),
                "delta": mean_d, "ci_lo": ci_lo, "ci_hi": ci_hi,
                "verdict": verdict,
            })

    cols = ["tool", "species", "n_paired",
            "mean_exon_f1_exp5", "mean_exon_f1_exp1",
            "delta", "ci_lo", "ci_hi", "verdict"]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                f"{r[c]:.4f}" if isinstance(r.get(c), float) else str(r.get(c, ""))
                for c in cols) + "\n")
    n_fail = sum(1 for r in rows if r["verdict"] == "fail")
    print(f"wrote {len(rows)} cells to {args.out}")
    print(f"verdict: pass={sum(1 for r in rows if r['verdict']=='pass')} "
          f"warn={sum(1 for r in rows if r['verdict']=='warn')} fail={n_fail}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
