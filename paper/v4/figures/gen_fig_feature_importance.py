"""Generate feature importance bar chart for TPRS paper.

All 25 features, horizontal bars = 4-species mean, overlaid dots = per-species
values (Mouse / Dog / Chicken / Zebrafish), colour-coded by category.
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# ── Load data ────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[3]
json_path = ROOT / "results/tprs/phase5/feature_importance_phase5.json"
with open(json_path) as f:
    data = json.load(f)

feats = data["sorted"]           # [(feature, mean_gain, category, short_name), ...] descending
per_sp = data["per_species"]     # {mouse: {feat: gain}, dog: {...}, ...}
SP_ORDER = ["mouse", "dog", "chicken", "zebrafish"]
SP_LABEL = ["Mouse", "Dog", "Chicken", "Zebrafish"]

# ── Colours ──────────────────────────────────────────────────────────────────
CAT_COLORS = {
    "Stage 1":      "#2166AC",
    "Geometry":     "#B2182B",
    "Query":        "#4DAF4A",
    "Stage 2":      "#FF7F00",
    "Species ind.": "#984EA3",
}
SP_MARKERS = {"mouse": "o", "dog": "s", "chicken": "^", "zebrafish": "D"}
SP_COLORS  = {"mouse": "#d73027", "dog": "#fc8d59", "chicken": "#4575b4", "zebrafish": "#91bfdb"}

# ── Build data arrays ────────────────────────────────────────────────────────
N = len(feats)  # 25
# Strip "X | " prefix from names (category already shown by bar colour)
raw_names = [f[3] for f in feats]
names  = [n.split(" | ", 1)[1] if " | " in n else n for n in raw_names]
means  = np.array([f[1] for f in feats])
# Remap A-stage→Stage 1, C-stage→Stage 2
cat_map = {"A-stage": "Stage 1", "C-stage": "Stage 2"}
cats   = [cat_map.get(f[2], f[2]) for f in feats]
colors = [CAT_COLORS.get(c, "#999999") for c in cats]

# Per-species values in same order
sp_vals = {}
for sp in SP_ORDER:
    sp_vals[sp] = np.array([per_sp[sp][feats[i][0]] for i in range(N)])

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 4.8))  # wide for 25 vertical bars

x = np.arange(N)
bar_w = 0.58

# Bars (mean) — vertical
bars = ax.bar(x, means, width=bar_w, color=colors, edgecolor="white",
              linewidth=0.3, zorder=2, alpha=0.88)

# Per-species dots + polyline across features (each species its own line)
# Offset each species' x position slightly so dots don't fully overlap
dot_offsets = {"mouse": +0.18, "dog": +0.06, "chicken": -0.06, "zebrafish": -0.18}
sp_x = {sp: x + dot_offsets[sp] for sp in SP_ORDER}

# Polyline per species connecting its dots across all 25 features
for sp in SP_ORDER:
    ax.plot(sp_x[sp], sp_vals[sp], color=SP_COLORS[sp], linewidth=0.9,
            zorder=3, alpha=0.55, linestyle="-")
    ax.scatter(sp_x[sp], sp_vals[sp], s=22, marker=SP_MARKERS[sp],
               c=SP_COLORS[sp], edgecolors="white", linewidths=0.3,
               zorder=5, label=SP_LABEL[SP_ORDER.index(sp)], clip_on=False)

# Annotate mean values above bars
for i, (v, bar) in enumerate(zip(means, bars)):
    if v > 0.003:
        ax.text(bar.get_x() + bar.get_width() / 2, v + 0.004,
                f"{v:.3f}", ha="center", fontsize=4.5, color="0.25")

# Axes
ax.set_xticks(x)
ax.set_xticklabels(names, fontsize=5.2, family="monospace", rotation=90, ha="center")
ax.set_ylabel("Gain (fraction of total loss reduction)", fontsize=7.5)
# ylim based on max species value (not just mean), to avoid clipping dots
max_sp_val = max(max(sp_vals[sp]) for sp in SP_ORDER)
ax.set_ylim(0, max_sp_val * 1.15)

ax.grid(axis="y", alpha=0.22, linewidth=0.4, linestyle="--")
ax.set_axisbelow(True)
ax.tick_params(axis="both", labelsize=5.5)
ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))

# Legend — categories
legend_cat = [
    plt.Rectangle((0, 0), 1, 1, fc=CAT_COLORS["Stage 1"],      ec="white", lw=0.2),
    plt.Rectangle((0, 0), 1, 1, fc=CAT_COLORS["Geometry"],     ec="white", lw=0.2),
    plt.Rectangle((0, 0), 1, 1, fc=CAT_COLORS["Query"],        ec="white", lw=0.2),
    plt.Rectangle((0, 0), 1, 1, fc=CAT_COLORS["Stage 2"],      ec="white", lw=0.2),
    plt.Rectangle((0, 0), 1, 1, fc=CAT_COLORS["Species ind."], ec="white", lw=0.2),
]
leg1 = ax.legend(legend_cat,
                 ["Stage 1 (8)", "Geometry (6)", "Query (2)", "Stage 2 (8)", "Species (1)"],
                 loc="upper right", fontsize=5.8, framealpha=0.85,
                 title="Feature category", title_fontsize=6.0,
                 ncol=3, columnspacing=0.5, handlelength=0.9, handletextpad=0.3)

# Legend — species markers
legend_sp = [
    plt.Line2D([0], [0], marker="o", c="w", markerfacecolor=SP_COLORS["mouse"],
               markersize=5, markeredgewidth=0.3, markeredgecolor="white", label="Mouse"),
    plt.Line2D([0], [0], marker="s", c="w", markerfacecolor=SP_COLORS["dog"],
               markersize=5, markeredgewidth=0.3, markeredgecolor="white", label="Dog"),
    plt.Line2D([0], [0], marker="^", c="w", markerfacecolor=SP_COLORS["chicken"],
               markersize=5, markeredgewidth=0.3, markeredgecolor="white", label="Chicken"),
    plt.Line2D([0], [0], marker="D", c="w", markerfacecolor=SP_COLORS["zebrafish"],
               markersize=5, markeredgewidth=0.3, markeredgecolor="white", label="Zebrafish"),
]
leg2 = ax.legend(handles=legend_sp, loc="upper left", fontsize=5.8,
                 framealpha=0.85, title="Species", title_fontsize=6.0,
                 ncol=4, columnspacing=0.5, handletextpad=0.2)
ax.add_artist(leg1)

# Title
ax.set_title("TPRS Phase 5 — XGBoost Feature Importance (bars = mean, lines = per-species profile)",
             fontsize=8.0, fontweight="bold", pad=6)

for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)

plt.tight_layout(pad=0.5)

# ── Save ─────────────────────────────────────────────────────────────────────
out_dir = ROOT / "paper/v4/figures"
out_dir.mkdir(parents=True, exist_ok=True)

for fmt in ("pdf", "png"):
    out_path = out_dir / f"fig_feature_importance_25.{fmt}"
    fig.savefig(out_path, dpi=300, bbox_inches="tight", pad_inches=0.05)
    print(f"Saved: {out_path}")

plt.close(fig)
print("Done.")
