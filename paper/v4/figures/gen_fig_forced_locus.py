"""Fig. Forced-Locus Diagnosis — RC1 (identity gap) + RC2 (search gap).

Panel (a): 4-bar groups per species.
  Left 2 bars = One-to-one reference (1:1 e2e / 1:1 forced)
  Right 2 bars = Paralog diagnosis (Paralog forced / Paralog e2e)
  1:1 e2e is background reference only; does NOT participate in gap calculation.
  RC1 = 1:1 forced − Paralog forced
  RC2 = Paralog forced − Paralog e2e
Panel (b): Correct-locus rate collapse, 1:1 vs paralog per tool per species.
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]

# ── Load data ────────────────────────────────────────────────────────────────
df_exp1 = pd.read_csv(os.path.join(ROOT, 'results', 'exp1', 'summary.tsv'), sep='\t')
df_exp5 = pd.read_csv(os.path.join(ROOT, 'results', 'exp5', 'cell', 'per_cell.tsv'), sep='\t')
df5 = df_exp5[(df_exp5['padding'] == 50000) & (df_exp5['tool'] == 'miniprot')]

SPECIES = ['mouse', 'dog', 'chicken', 'zebrafish']
SP_LABELS = {'mouse': 'Mouse', 'dog': 'Dog', 'chicken': 'Chicken', 'zebrafish': 'Zebrafish'}

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(6.8, 3.1))

x = np.arange(len(SPECIES))
bar_w = 0.18
gap_between_groups = 0.12   # extra space between 1:1 and Paralog groups

# Offset positions for 4 bars — two groups with gap
# Group 1 (One-to-one):   −gap/2 − 0.5·w   and   −gap/2 + 0.5·w
# Group 2 (Paralog):      +gap/2 − 0.5·w   and   +gap/2 + 0.5·w
g = gap_between_groups
offsets = [
    -g/2 - 0.5*bar_w,   # 1:1 e2e        (leftmost, reference)
    -g/2 + 0.5*bar_w,   # 1:1 forced
    +g/2 - 0.5*bar_w,   # Paralog forced
    +g/2 + 0.5*bar_w,   # Paralog e2e    (rightmost, reference)
]

# 1:1 e2e = muted/hatched (background ref), 1:1 forced = deep blue
# Paralog forced = orange, Paralog e2e = red
colors    = ['#BBDEF0', '#2166AC', '#F4A582', '#B2182B']
hatches   = ['//',       None,      None,      None]       # hatch 1:1 e2e only
alphas    = [0.55,       1.0,       1.0,       0.85]       # fade 1:1 e2e
edgecolor = ['#6BAED6',  '#2166AC', '#F4A582', '#B2182B']
labels    = ['1:1 e2e\n(ref.)', '1:1 forced', 'Paralog forced', 'Paralog e2e\n(ref.)']

# Store per-species values for gap annotation
gap_data = []

for i, sp in enumerate(SPECIES):
    r1 = df_exp1[(df_exp1['species'] == sp) & (df_exp1['tool'] == 'miniprot') & (df_exp1['regime'] == 'one_to_one')]
    f1_1e2e = r1['exon_f1_mean'].values[0] if len(r1) > 0 else 0
    r2 = df5[(df5['species'] == sp) & (df5['regime'] == 'one2one')]
    f1_1fl = r2['mean_exon_f1'].values[0] if len(r2) > 0 else 0
    r3 = df5[(df5['species'] == sp) & (df5['regime'] == 'paralog_broad')]
    f1_pfl = r3['mean_exon_f1'].values[0] if len(r3) > 0 else 0
    r4 = df_exp1[(df_exp1['species'] == sp) & (df_exp1['tool'] == 'miniprot') & (df_exp1['regime'] == 'paralog')]
    f1_pe2e = r4['exon_f1_mean'].values[0] if len(r4) > 0 else 0

    vals = [f1_1e2e, f1_1fl, f1_pfl, f1_pe2e]

    for j, (val, off, col, ha, al, ec, lab) in enumerate(
        zip(vals, offsets, colors, hatches, alphas, edgecolor, labels)):
        ax_a.bar(i + off, val, bar_w, color=col, alpha=al, hatch=ha,
                 edgecolor=ec, linewidth=0.6, zorder=2,
                 label=lab if i == 0 else '')

    # ── Value labels: only on diagnostic bars (1:1 forced, paralog forced, paralog e2e) ──
    for j, (val, off, show) in enumerate(zip(vals, offsets,
        [False, True, True, True])):  # skip value label on 1:1 e2e (reference)
        if not show or val <= 0:
            continue
        ax_a.text(i + off, val + 0.014, f'{val:.2f}',
                  fontsize=4.2, ha='center', va='bottom', color='0.2')

    identity_gap = f1_1fl - f1_pfl   # RC1
    search_gap  = f1_pfl - f1_pe2e   # RC2

    # ── RC1 bracket: spans between 1:1 forced (top) and Paralog forced (bottom) ──
    x_rc1 = i + (offsets[1] + offsets[2]) / 2   # midpoint between the two bars
    y_rc1_top, y_rc1_bot = f1_1fl, f1_pfl
    y_rc1_mid = (y_rc1_top + y_rc1_bot) / 2
    ax_a.plot([x_rc1, x_rc1], [y_rc1_bot, y_rc1_top],
              color='#2166AC', linewidth=1.3, zorder=6, solid_capstyle='butt')
    # Ticks at top and bottom
    tick_w = 0.05
    ax_a.plot([x_rc1 - tick_w, x_rc1], [y_rc1_bot, y_rc1_bot],
              color='#2166AC', linewidth=0.8, zorder=6)
    ax_a.plot([x_rc1 - tick_w, x_rc1], [y_rc1_top, y_rc1_top],
              color='#2166AC', linewidth=0.8, zorder=6)
    ax_a.text(x_rc1 + 0.16, y_rc1_mid, f'RC1\n{identity_gap:+.2f}',
              fontsize=5.2, ha='left', va='center', color='#2166AC', fontweight='bold')

    # ── RC2 bracket: spans between Paralog forced (top) and Paralog e2e (bottom) ──
    x_rc2 = i + (offsets[2] + offsets[3]) / 2
    y_rc2_top, y_rc2_bot = f1_pfl, f1_pe2e
    y_rc2_mid = (y_rc2_top + y_rc2_bot) / 2
    ax_a.plot([x_rc2, x_rc2], [y_rc2_bot, y_rc2_top],
              color='#B2182B', linewidth=1.3, zorder=6, solid_capstyle='butt')
    tick_w2 = 0.04
    ax_a.plot([x_rc2 - tick_w2, x_rc2], [y_rc2_bot, y_rc2_bot],
              color='#B2182B', linewidth=0.8, zorder=6)
    ax_a.plot([x_rc2 - tick_w2, x_rc2], [y_rc2_top, y_rc2_top],
              color='#B2182B', linewidth=0.8, zorder=6)
    ax_a.text(x_rc2 + 0.13, y_rc2_mid, f'RC2\n{search_gap:+.2f}',
              fontsize=5.2, ha='left', va='center', color='#B2182B', fontweight='bold')

    gap_data.append({'sp': sp, 'f1_1e2e': f1_1e2e, 'f1_1fl': f1_1fl,
                     'f1_pfl': f1_pfl, 'f1_pe2e': f1_pe2e,
                     'RC1': identity_gap, 'RC2': search_gap})

# ── Group labels above bars ──────────────────────────────────────────────────
y_group = 1.12
for s_i, sp in enumerate(SPECIES):
    left_1to1  = s_i + offsets[0] - bar_w/2
    right_1to1 = s_i + offsets[1] + bar_w/2
    left_para  = s_i + offsets[2] - bar_w/2
    right_para = s_i + offsets[3] + bar_w/2

    # bracket under group labels
    for lx, rx, txt, c in [
        (left_1to1, right_1to1, '1:1', '#2166AC'),
        (left_para, right_para, 'Paralog', '#B2182B'),
    ]:
        ax_a.plot([lx, rx], [y_group, y_group], color=c, linewidth=1.2, zorder=1)
        ax_a.plot([lx, lx], [y_group - 0.03, y_group], color=c, linewidth=0.7, zorder=1)
        ax_a.plot([rx, rx], [y_group - 0.03, y_group], color=c, linewidth=0.7, zorder=1)
        ax_a.text((lx + rx) / 2, y_group + 0.025, txt,
                  fontsize=6.0, ha='center', va='bottom', color=c, fontweight='bold')

# ── Axes & legend ────────────────────────────────────────────────────────────
ax_a.set_xticks(x)
ax_a.set_xticklabels([SP_LABELS[s] for s in SPECIES], fontsize=7)
ax_a.set_ylim(0, 1.25)
ax_a.set_ylabel('Exon F1', fontsize=8)
ax_a.legend(frameon=False, fontsize=4.8, loc='upper right', ncol=2,
            columnspacing=0.4, handletextpad=0.3, handlelength=1.0)
ax_a.set_title('(a) Forced-locus diagnosis with miniprot', fontsize=8, fontweight='bold')
ax_a.grid(axis='y', alpha=0.15, linewidth=0.4)

# Mean RC annotation box
mRC1 = np.mean([g['RC1'] for g in gap_data])
mRC2 = np.mean([g['RC2'] for g in gap_data])
ax_a.text(0.98, 0.48, f'RC1 = {mRC1:+.2f}\nRC2 = {mRC2:+.2f}',
          transform=ax_a.transAxes, fontsize=5.5, ha='right', va='center',
          bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFFDE7',
                    edgecolor='#F9A825', lw=0.5))

# ── Panel (b): Correct-locus rate collapse ────────────────────────────────────
TOOLS_DISPLAY = ['miniprot', 'Spaln3', 'GTH', 'Scipio', 'genBlastG', 'miniSplice']
TOOLS_DATA = ['miniprot', 'spaln3', 'gth', 'scipio', 'genblastg', 'minisplice_miniprot']
TOOL_COLORS = {
    'miniprot': '#2166AC', 'spaln3': '#B2182B', 'gth': '#4DAF4A',
    'scipio': '#FF7F00', 'genblastg': '#984EA3', 'minisplice_miniprot': '#F781BF',
}

x2 = np.arange(len(SPECIES))
width2 = 0.11

for j, (td, tl) in enumerate(zip(TOOLS_DATA, TOOLS_DISPLAY)):
    offset = (j - len(TOOLS_DATA) / 2 + 0.5) * width2
    o2o, pl = [], []
    for sp in SPECIES:
        r1 = df_exp1[(df_exp1['species'] == sp) & (df_exp1['tool'] == td) & (df_exp1['regime'] == 'one_to_one')]
        rp = df_exp1[(df_exp1['species'] == sp) & (df_exp1['tool'] == td) & (df_exp1['regime'] == 'paralog')]
        o2o.append(r1['correct_locus_rate'].values[0] if len(r1) > 0 else 0)
        pl.append(rp['correct_locus_rate'].values[0] if len(rp) > 0 else 0)
    c = TOOL_COLORS.get(td, '#999999')
    ax_b.bar(x2 + offset, o2o, width2, color=c, alpha=0.95, label=tl, zorder=2)
    ax_b.bar(x2 + offset, pl, width2, color=c, alpha=0.25, hatch='////',
             edgecolor=c, linewidth=0.3, zorder=2)

p1 = mpatches.Patch(facecolor='0.3', alpha=0.95, label='1:1 ortholog')
p2 = mpatches.Patch(facecolor='0.3', alpha=0.25, hatch='////', edgecolor='0.3', label='Paralog')
leg1 = ax_b.legend(handles=[p1, p2], loc='upper right', frameon=False, fontsize=5.5,
                   title='Regime', title_fontsize=6)
ax_b.add_artist(leg1)
ax_b.legend(frameon=False, fontsize=5, loc='lower left', ncol=3,
            columnspacing=0.4, handletextpad=0.3)

ax_b.set_xticks(x2)
ax_b.set_xticklabels([SP_LABELS[s] for s in SPECIES], fontsize=7)
ax_b.set_ylim(0, 1.08)
ax_b.set_ylabel('Correct-locus rate', fontsize=8)
ax_b.set_title('(b) Locus confusion', fontsize=8, fontweight='bold')
ax_b.grid(axis='y', alpha=0.15, linewidth=0.4)

fig.tight_layout(pad=0.7)

# ── Save ─────────────────────────────────────────────────────────────────────
out_dir = os.path.join(os.path.dirname(__file__))
for fmt in ('pdf', 'png'):
    path = os.path.join(out_dir, f'fig3_forced_locus_diagnosis.{fmt}')
    fig.savefig(path, dpi=300, bbox_inches='tight', pad_inches=0.05)
    print(f'Saved: {path}')
print('Done.')
