"""Fig 4: Ablation forest plot (cleaner version)

- Shared y labels on the left panel only
- Same component order in both panels
- Group labels (Stage I / Stage II) placed in outer margin
- Cleaner legend and lighter negligible band
"""

import os
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pathlib import Path

# 如果你自己的样式文件会覆盖太多字体/边距，可以先注释掉这一行
from paper_plot_style import *

ROOT = Path(__file__).resolve().parents[3]
with open(os.path.join(ROOT, 'results', 'tprs', 'phase5',
                       'ablation_summary_all_species.json')) as f:
    data = json.load(f)

SPECIES = ['mouse', 'dog', 'chicken', 'zebrafish']
SP_COLORS_LIST = [SPECIES_COLORS[s] for s in SPECIES]
SP_LABELS_LIST = [SPECIES_LABELS[s] for s in SPECIES]

# ── Same component order in both panels ─────────────────────────────────────
COMPONENTS = [
    ('LAST',               'E6_last',      'I'),
    ('MMseqs2',            'E6_mmseqs',    'I'),
    ('DIAMOND',            'E6_diamond',   'I'),
    ('miniprot',           'E6_miniprot',  'I'),
    ('Species protein',    'E2',           'I'),
    ('Orthology prior',    'E4',           'I'),
    ('Router',             'E3',           'II'),
    ('Alignment features', 'E5',           'II'),
]

N = len(COMPONENTS)
ypos = np.arange(N)[::-1]   # top to bottom

METRICS = [
    ('best_overlap', 'Panel A: best-overlap F1'),
    ('primary',      'Panel B: primary F1'),
]

# ── Build delta dictionaries ────────────────────────────────────────────────
def build_metric_dict(metric_name):
    out = {}
    for label, key, stage in COMPONENTS:
        sp_deltas = {}
        sp_ci_lo = {}
        sp_ci_hi = {}
        all_deltas = []

        for sp in SPECIES:
            e1 = data[sp][metric_name]['E1']['mean']
            ab = data[sp][metric_name][key]

            d = ab['mean'] - e1
            lo = ab['ci_lower'] - e1
            hi = ab['ci_upper'] - e1

            sp_deltas[sp] = d
            sp_ci_lo[sp] = lo
            sp_ci_hi[sp] = hi
            all_deltas.append(d)

        out[label] = {
            'deltas': sp_deltas,
            'ci_lo': sp_ci_lo,
            'ci_hi': sp_ci_hi,
            'mean': float(np.mean(all_deltas)),
            'stage': stage,
        }
    return out

comp_data_bo = build_metric_dict('best_overlap')
comp_data_pr = build_metric_dict('primary')

# ── Figure layout ───────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.2), sharey=True)
ax_bo, ax_pr = axes

# 留出左边距放组件名和 Stage 标签
fig.subplots_adjust(left=0.27, right=0.98, top=0.90, bottom=0.20, wspace=0.10)

# Stage separator: after Orthology prior
sep_y = ypos[5] - 0.5   # between row 5 and row 6

# y tick labels only on left panel
comp_labels = [x[0] for x in COMPONENTS]
ax_bo.set_yticks(ypos)
ax_bo.set_yticklabels(comp_labels, fontsize=8)
ax_pr.set_yticks(ypos)
ax_pr.tick_params(axis='y', left=False, labelleft=False)

# ── Outer-margin stage labels (cleaner than inside label column) ───────────
# use y-axis transform: x in axes coords, y in data coords
stage1_center = np.mean([ypos[i] for i in range(6)])
stage2_center = np.mean([ypos[i] for i in range(6, 8)])

ax_bo.text(-0.42, stage1_center, 'Stage I\nCandidate\nGeneration',
           transform=ax_bo.get_yaxis_transform(),
           ha='left', va='center', fontsize=8,
           color='#2166AC', fontweight='bold', linespacing=1.15)

ax_bo.text(-0.42, stage2_center, 'Stage II\nLearned\nSelection',
           transform=ax_bo.get_yaxis_transform(),
           ha='left', va='center', fontsize=8,
           color='#D94801', fontweight='bold', linespacing=1.15)

# ── Helper to draw one panel ────────────────────────────────────────────────
def draw_panel(ax, comp_data, title, xlabel, xlim):
    ax.set_ylim(-0.8, N - 0.2)

    # negligible band and zero line
    ax.axvspan(-0.05, 0.05, color='#BDBDBD', alpha=0.10, zorder=0)
    ax.axvline(0, color='#888888', linewidth=0.8, zorder=1)

    # stage separator
    ax.axhline(sep_y, color='#BDBDBD', linewidth=1.0, zorder=1)

    # per-component
    for i, (lbl, key, stage) in enumerate(COMPONENTS):
        yi = ypos[i]
        cd = comp_data[lbl]

        # species points + CI
        for j, sp in enumerate(SPECIES):
            d = cd['deltas'][sp]
            lo = cd['ci_lo'][sp]
            hi = cd['ci_hi'][sp]
            color = SP_COLORS_LIST[j]

            yoff = (j - 1.5) * 0.11
            ax.errorbar(
                d, yi + yoff,
                xerr=[[d - lo], [hi - d]],
                fmt='o',
                ms=4.8,
                color=color,
                ecolor=color,
                elinewidth=0.9,
                capsize=2.0,
                capthick=0.7,
                mec='white',
                mew=0.35,
                alpha=0.95,
                zorder=3
            )

        # mean diamond
        mean_d = cd['mean']
        ax.scatter(mean_d, yi, marker='D', s=48,
                   c='#3A3A3A', edgecolors='white', linewidths=0.6, zorder=4)
        ax.text(mean_d, yi - 0.28, f'{mean_d:+.2f}',
                ha='center', va='top', fontsize=6.5,
                color='#3A3A3A', fontweight='bold')

    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_title(title, fontsize=10, fontweight='bold', pad=6)

    # style
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', length=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# ── Draw both panels ────────────────────────────────────────────────────────
draw_panel(
    ax_bo, comp_data_bo,
    title='Panel A: best-overlap F1',
    xlabel='Δ best-overlap F1 from Full TPRS',
    xlim=(-0.52, 0.12)
)

draw_panel(
    ax_pr, comp_data_pr,
    title='Panel B: primary F1',
    xlabel='Δ primary F1 from Full TPRS',
    xlim=(-0.67, 0.12)
)

# ── Legend ──────────────────────────────────────────────────────────────────
legend_species = [
    plt.Line2D([0], [0],
               marker='o', color='w',
               markerfacecolor=SP_COLORS_LIST[j],
               markeredgecolor='white',
               markeredgewidth=0.5,
               markersize=6.5,
               label=SP_LABELS_LIST[j])
    for j in range(4)
]
legend_mean = plt.Line2D([0], [0],
                         marker='D', color='w',
                         markerfacecolor='#3A3A3A',
                         markeredgecolor='white',
                         markersize=7,
                         label='Mean Δ')

fig.legend(handles=legend_species + [legend_mean],
           frameon=False, fontsize=8,
           loc='lower center', ncol=5, bbox_to_anchor=(0.5, 0.03))

fig.text(0.5, 0.12, 'Gray band: |Δ| < 0.05 (negligible)',
         ha='center', fontsize=7.2, color='#7A7A7A', style='italic')

# ── Save ────────────────────────────────────────────────────────────────────
out_dir = os.path.dirname(__file__)
for fmt in ('pdf', 'png'):
    path = os.path.join(out_dir, f'fig4_ablation_forest.{fmt}')
    fig.savefig(path, dpi=300, bbox_inches='tight', pad_inches=0.03)
    print(f'Saved: {path}')

print('Fig 4 done.')
