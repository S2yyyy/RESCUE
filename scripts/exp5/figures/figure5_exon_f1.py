"""Figure 5 — exon F1 by tool, 4 species panels.

Inputs:  results/exp5/cell/per_cell.tsv (one2one, padding=50000)
Outputs: results/exp5/figures/figure5_exon_f1.{pdf,png}

Bug-2 compliance: bars show mean exon F1 over successful queries only.
A second hatched bar overlays frac_failed for transparency.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from rescue_paths import PROJECT_ROOT  # noqa: E402

REPO = PROJECT_ROOT
IN = REPO / "results/exp5/cell/per_cell.tsv"
OUT_DIR = REPO / "results/exp5/figures"

TOOLS = ["miniprot", "spaln3", "gth", "exonerate",
         "genewise", "scipio", "genblastg", "minisplice_miniprot"]
SPECIES = ["mouse", "dog", "chicken", "zebrafish"]


def main() -> int:
    rows = list(csv.DictReader(open(IN), delimiter="\t"))
    rows = [r for r in rows
            if r["regime"] == "one2one" and int(r["padding"]) == 50000]
    idx = {(r["tool"], r["species"]): r for r in rows}

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5), sharey=True)
    x = np.arange(len(TOOLS))
    bar_w = 0.8

    for ax, sp in zip(axes, SPECIES):
        means = []
        lo = []
        hi = []
        failed = []
        for t in TOOLS:
            r = idx.get((t, sp))
            if r is None:
                means.append(np.nan); lo.append(np.nan); hi.append(np.nan)
                failed.append(np.nan)
                continue
            means.append(float(r["mean_exon_f1"]))
            lo.append(float(r["ci_lo"]))
            hi.append(float(r["ci_hi"]))
            failed.append(float(r["frac_failed"]))
        means = np.array(means); lo = np.array(lo); hi = np.array(hi)
        failed = np.array(failed)

        yerr = np.vstack([means - lo, hi - means])
        ax.bar(x, means, bar_w, yerr=yerr, capsize=3,
               color="#2E7D32", label="mean exon F1 (succ.)")
        # Overlay failure-rate as hatched bar from top
        ax.bar(x, failed, bar_w, bottom=1.0, color="none",
               edgecolor="#C62828", hatch="//", label="frac_failed")
        ax.set_xticks(x)
        ax.set_xticklabels(TOOLS, rotation=30, ha="right", fontsize=8)
        ax.set_title(sp)
        ax.set_ylim(0, 1.15)
        ax.axhline(1.0, color="black", lw=0.5)
        ax.grid(alpha=0.3, axis="y")
    axes[0].set_ylabel("mean exon F1   |   frac_failed (above 1.0 line)")
    axes[0].legend(loc="lower left", fontsize=7)
    fig.suptitle("Exp 5 Figure 5 — alignment-kernel exon F1 (forced GT locus, ±50 kb pad)",
                 fontsize=11)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        out = OUT_DIR / f"figure5_exon_f1.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
