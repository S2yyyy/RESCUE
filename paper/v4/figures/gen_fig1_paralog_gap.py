"""Fig 1: The Paralog Gap — primary F1 vs best-overlap F1 on paralogs."""
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
from pathlib import Path
from paper_plot_style import *

ROOT = Path(__file__).resolve().parents[3]

# Existing tools: primary and best_overlap F1 from bootstrap CI CSV
ci_csv = pd.read_csv(f'{ROOT}/paper/v2/figures/out/table_exp1_benchmark.csv')

# TPRS: primary and best_overlap from Phase 5 E1
with open(f'{ROOT}/results/tprs/phase5/ablation_summary_all_species.json') as f:
    tprs = json.load(f)

TOOLS_ORDER = ['miniprot', 'Spaln3', 'GTH', 'Scipio', 'genBlastG', 'mini+miniSplice']
SPECIES_ORDER = ['mouse', 'dog', 'chicken', 'zebrafish']

fig, axes = plt.subplots(1, 4, figsize=(6.8, 2.8), sharey=True)

for ax, species in zip(axes, SPECIES_ORDER):
    x = np.arange(len(TOOLS_ORDER) + 1)  # +1 for TPRS
    width = 0.30

    # Existing tools: primary (solid) and best_overlap (hatched) on paralogs
    primary_vals = []
    bo_vals = []
    for tool in TOOLS_ORDER:
        r = ci_csv[(ci_csv['tool'] == tool) & (ci_csv['sp_key'] == species)]
        primary_vals.append(r['para_prim_f1'].values[0] if len(r) > 0 else 0)
        bo_vals.append(r['para_bo_f1'].values[0] if len(r) > 0 else 0)

    x_tools = x[:-1]
    ax.bar(x_tools - width/2, primary_vals, width, color=COLORS['gray'],
           edgecolor='#333333', linewidth=0.4, label='Primary F1')
    ax.bar(x_tools + width/2, bo_vals, width, color=COLORS['light_gray'],
           edgecolor=COLORS['gray'], linewidth=0.4, hatch='///',
           label='Best-overlap F1')

    # TPRS: primary and best_overlap (both gold)
    tprs_p = tprs[species]['primary']['E1']['mean']
    tprs_bo = tprs[species]['best_overlap']['E1']['mean']
    ax.bar([x[-1] - width/2], [tprs_p], width, color='#B8860B',
           edgecolor='#8B6508', linewidth=0.6, label='TPRS primary')
    ax.bar([x[-1] + width/2], [tprs_bo], width, color='#DAA520',
           edgecolor='#8B6508', linewidth=0.4, alpha=0.7, hatch='...',
           label='TPRS best-overlap')

    ax.set_title(SPECIES_LABELS[species], fontsize=8, pad=4)
    ax.set_xticks(x)
    ax.set_xticklabels(TOOLS_ORDER + ['TPRS'], rotation=30, ha='right', fontsize=6.5)
    ax.set_ylim(0, 1.05)
    if species == 'mouse':
        ax.set_ylabel('Exon F1 (paralog)', fontsize=8)

# Single compact legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=COLORS['gray'], edgecolor='#333333', label='Primary F1'),
    Patch(facecolor=COLORS['light_gray'], edgecolor=COLORS['gray'], hatch='///',
          label='Best-overlap F1'),
    Patch(facecolor='#B8860B', edgecolor='#8B6508', label='TPRS primary'),
    Patch(facecolor='#DAA520', edgecolor='#8B6508', alpha=0.7, hatch='...',
          label='TPRS best-overlap'),
]
fig.legend(handles=legend_elements, frameon=False, fontsize=5.8,
           loc='upper center', ncol=4, bbox_to_anchor=(0.5, 1.02))

fig.tight_layout(pad=0.8, rect=(0, 0, 1, 0.88))
save_fig(fig, 'fig1_paralog_gap')
print('Fig 1 done.')
