"""Find the (A, C) pair that is universally near-optimal on all 4 paralog cells.

Two rankings:
  - mean   : highest mean F1 across 4 paralog cells (average winner)
  - minimax: smallest worst-case gap to per-cell best (most robust)
"""
from __future__ import annotations

import numpy as np

from pathlib import Path

REPO = Path(".").resolve()
HYBRID = REPO / "results/exp6/hybrid"

SPECIES = ["mouse", "dog", "chicken", "zebrafish"]
A_TOOLS = ["miniprot", "spaln3", "genblastg", "tblastn", "diamond",
           "mmseqs", "last", "genblasta"]
C_TOOLS = ["miniprot", "spaln3", "gth", "exonerate", "genewise", "scipio", "genblastg"]


def main() -> int:
    matrices = []
    cell_best = []
    for sp in SPECIES:
        npz = np.load(HYBRID / f"bootstrap_{sp}.paralog_broad.npz")
        H = npz["hybrid_orig"]    # (8, 7)
        matrices.append(H)
        cell_best.append(float(H.max()))
    M = np.stack(matrices)       # (4, 8, 7)
    best = np.array(cell_best)   # (4,)

    # For every (A, C): F1 vector across 4 cells + gap to best
    print("=" * 92)
    print("Per-cell F1 of each pair, plus mean and worst gap to per-cell best")
    print(f"{'pair':<26} {'mouse':>7} {'dog':>7} {'chick':>7} {'zfish':>7}  "
          f"{'mean':>6} {'worst-gap':>10}")
    print("-" * 92)

    rows = []
    for ai, atool in enumerate(A_TOOLS):
        for ci, ctool in enumerate(C_TOOLS):
            f1_vec = M[:, ai, ci]
            gap_vec = best - f1_vec        # per cell, ≥ 0
            mean_f1 = f1_vec.mean()
            worst_gap = gap_vec.max()
            rows.append((mean_f1, worst_gap, atool, ctool, f1_vec, gap_vec))

    print("\nTop-10 by MEAN F1 across 4 paralog cells:")
    for mean_f1, worst_gap, atool, ctool, f1_vec, gap_vec in sorted(rows, reverse=True)[:10]:
        cross = "*" if atool != ctool else " "
        print(f"  {cross} {atool}+{ctool:<14}  "
              f"{f1_vec[0]:>7.3f} {f1_vec[1]:>7.3f} {f1_vec[2]:>7.3f} {f1_vec[3]:>7.3f}  "
              f"{mean_f1:>6.3f} {worst_gap:>10.3f}")

    print("\nTop-10 by SMALLEST WORST-CASE gap to per-cell best (minimax robust):")
    for mean_f1, worst_gap, atool, ctool, f1_vec, gap_vec in sorted(
            rows, key=lambda r: r[1])[:10]:
        cross = "*" if atool != ctool else " "
        print(f"  {cross} {atool}+{ctool:<14}  "
              f"{f1_vec[0]:>7.3f} {f1_vec[1]:>7.3f} {f1_vec[2]:>7.3f} {f1_vec[3]:>7.3f}  "
              f"{mean_f1:>6.3f} {worst_gap:>10.3f}")

    print("\n(* = cross-tool pair; per-cell best = "
          f"mouse {best[0]:.3f}, dog {best[1]:.3f}, "
          f"chicken {best[2]:.3f}, zebrafish {best[3]:.3f})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
