"""Fig 2: TPRS Algorithm Overview — two-stage unified diagram.
Stage 1: Candidate Generation (3 dimensions)
Stage 2: Learned Joint Selection
"""
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
from paper_plot_style import *

fig, ax = plt.subplots(1, 1, figsize=(7.0, 7.2))
ax.set_xlim(0, 12)
ax.set_ylim(0, 10.5)
ax.axis('off')

# Color palette
C_STAGE1 = '#2166AC'       # blue — Stage 1
C_STAGE1_LIGHT = '#D0E4F5'
C_STAGE2 = '#D94801'       # orange — Stage 2
C_STAGE2_LIGHT = '#FDD0A2'
C_DIM1 = '#2E7D32'         # green — tool diversity
C_DIM2 = '#7B4FB8'         # purple — query identity
C_DIM3 = '#F57F17'         # amber — prior knowledge
C_POOL = '#4DB6AC'         # teal — candidate pool
C_OUTPUT = '#2E7D32'       # green — output
C_DARK = '#333333'
C_GRAY = '#666666'
C_LIGHT = '#F5F5F5'

def add_box(ax, x, y, w, h, label, subtitle=None, facecolor='white', edgecolor=C_GRAY,
            linewidth=1.0, fontsize=7, fontweight='bold', textcolor=C_DARK, zorder=2,
            subfontsize=5.5):
    box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                          facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, zorder=zorder)
    ax.add_patch(box)
    ax.text(x + w/2, y + h - 0.18, label, fontsize=fontsize, ha='center', va='top',
            fontweight=fontweight, color=textcolor, zorder=zorder+1)
    if subtitle:
        lines = subtitle.split('\n')
        for li, line in enumerate(lines):
            ax.text(x + w/2, y + h - 0.45 - li*0.22, line, fontsize=subfontsize, ha='center',
                    va='top', color=C_GRAY, zorder=zorder+1)

def arrow(ax, x1, y1, x2, y2, color=C_GRAY, lw=0.8):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw))

# ============================================================
# STAGE 1 — Candidate Generation (upper half, y=5.5 to 10.0)
# ============================================================

# Stage 1 banner
s1_banner = FancyBboxPatch((0.3, 9.5), 11.4, 0.65, boxstyle="round,pad=0.1",
                            facecolor=C_STAGE1_LIGHT, edgecolor=C_STAGE1, linewidth=1.5, zorder=2)
ax.add_patch(s1_banner)
ax.text(6.0, 9.95, 'Stage 1: Candidate Generation', fontsize=9, ha='center', va='center',
        fontweight='bold', color=C_STAGE1, zorder=3)
ax.text(6.0, 9.68, 'Expand search space along three dimensions to ensure correct (locus, query-type) pair enters the candidate pool',
        fontsize=5.8, ha='center', va='center', color=C_DARK, style='italic', zorder=3)

# Input
add_box(ax, 0.8, 8.3, 2.4, 0.9, 'Input',
        'Human query protein\n+ target genome + CDS annotations', edgecolor=C_DARK)

# ---- Dimension 1: Tool Diversity (left) ----
add_box(ax, 4.0, 8.3, 2.4, 0.9, 'Dim 1: Tool Diversity',
        'diamond / last / miniprot / mmseqs2\nMerge overlapping hits\nChrom+strand clustering, ≤50 kb', edgecolor=C_DIM1)

# ---- Dimension 2: Query Identity (middle) ----
add_box(ax, 7.0, 8.3, 2.4, 0.9, 'Dim 2: Query Identity',
        'Retrieve species-matched protein\nRun miniprot with BOTH\nhuman query AND species protein', edgecolor=C_DIM2)

# ---- Dimension 3: Prior Knowledge (right) ----
add_box(ax, 10.0, 8.3, 1.5, 0.9, 'Dim 3:\nPrior',
        'Ensembl Compara\nco-ortholog\nrecovery\n(733/735, 99.7%)', edgecolor=C_DIM3, fontsize=6, subfontsize=4.8)

# Phase labels on top of each dim
ax.text(5.2, 9.05, 'Phase 1', fontsize=5.5, ha='center', color=C_DIM1, fontweight='bold', style='italic')
ax.text(8.2, 9.05, 'Phase 3', fontsize=5.5, ha='center', color=C_DIM2, fontweight='bold', style='italic')
ax.text(10.75, 9.05, 'Phase 2', fontsize=5.5, ha='center', color=C_DIM3, fontweight='bold', style='italic')

# Arrows: Input → each dimension
arrow(ax, 3.2, 8.75, 4.0, 8.75, C_GRAY)
arrow(ax, 3.2, 8.75, 7.0, 8.75, C_GRAY)
arrow(ax, 3.2, 8.75, 10.0, 8.75, C_GRAY)

# Candidate pool (bottom of Stage 1)
pool_box = FancyBboxPatch((3.5, 6.5), 5.5, 0.9, boxstyle="round,pad=0.1",
                           facecolor='#E0F2F1', edgecolor=C_POOL, linewidth=1.5, zorder=2)
ax.add_patch(pool_box)
ax.text(6.25, 7.2, 'Candidate Pool', fontsize=8, ha='center', fontweight='bold', color='#00695C', zorder=3)
ax.text(6.25, 6.85, r'$\mathcal{C} = \{(L_1, \mathrm{human}), (L_1, \mathrm{species}), (L_2, \mathrm{human}), \dots\}$',
        fontsize=6, ha='center', color=C_DARK, zorder=3)
ax.text(6.25, 6.6, '~9–17 (locus, query-type) pairs per protein',
        fontsize=5.2, ha='center', color=C_GRAY, style='italic', zorder=3)

# Arrows: dims → pool
arrow(ax, 5.2, 8.3, 6.25, 7.4, C_DIM1)
arrow(ax, 8.2, 8.3, 6.25, 7.4, C_DIM2)
arrow(ax, 10.75, 8.3, 9.0, 7.4, C_DIM3)

# ============================================================
# Connecting arrow between stages
# ============================================================
arrow(ax, 6.25, 6.5, 6.25, 5.4, C_DARK, lw=1.5)
ax.text(6.9, 5.95, 'Each candidate paired\nwith miniprot alignment\n→ 24-dim feature vector',
        fontsize=5.2, ha='left', va='center', color=C_GRAY, style='italic')

# ============================================================
# STAGE 2 — Learned Joint Selection (lower half, y=1.0 to 5.1)
# ============================================================

# Stage 2 banner
s2_banner = FancyBboxPatch((0.3, 5.0), 11.4, 0.65, boxstyle="round,pad=0.1",
                            facecolor=C_STAGE2_LIGHT, edgecolor=C_STAGE2, linewidth=1.5, zorder=2)
ax.add_patch(s2_banner)
ax.text(6.0, 5.45, 'Stage 2: Learned Joint Selection', fontsize=9, ha='center', va='center',
        fontweight='bold', color=C_STAGE2, zorder=3)
ax.text(6.0, 5.18, 'Score every (locus, query-type) pair jointly; select the best via learned regression',
        fontsize=5.8, ha='center', va='center', color=C_DARK, style='italic', zorder=3)

# Feature extraction — 4 blocks
feat_w = 2.2
feat_h = 1.2
feat_y = 3.5
feat_colors = [C_STAGE1, C_DIM2, C_DIM1, C_POOL]

# A-stage (8)
add_box(ax, 0.8, feat_y, feat_w, feat_h, 'A-stage (8)',
        'Tool presence ×4\nBest score / qcov\nRaw hit count\nN tools hitting', edgecolor=feat_colors[0], fontsize=6.5, subfontsize=4.8)
# Geometry (6)
add_box(ax, 3.3, feat_y, feat_w, feat_h, 'Geometry (6)',
        'Region size (lin+log)\nOn primary chrom?\nStrand agreement\nCo-ortholog overlap', edgecolor=feat_colors[1], fontsize=6.5, subfontsize=4.8)
# Query (2)
add_box(ax, 5.8, feat_y, feat_w*0.8, feat_h, 'Query (2)',
        'Prot. length\n(linear + log)', edgecolor=feat_colors[2], fontsize=6.5, subfontsize=4.8)
# C-stage (8)
add_box(ax, 7.6, feat_y, feat_w, feat_h, 'C-stage (8)',
        'c_identity, coverage\nc_n_exons, CDS length\nc_composite_score\nc_has_hit, hit_score, rank', edgecolor=feat_colors[3], fontsize=6.5, subfontsize=4.8)
# Species indicator
add_box(ax, 10.1, feat_y, 1.2, feat_h, 'Ind. (1)',
        'Species\nindicator', edgecolor=C_DARK, fontsize=6, subfontsize=4.8)

# Merge arrow → 24-dim vector
ax.annotate('', xy=(6.0, 2.8), xytext=(6.0, 3.5),
            arrowprops=dict(arrowstyle='->', color=C_DARK, lw=1.2))

ax.text(6.0, 3.15, '24-dimensional feature vector', fontsize=6, ha='center', va='center',
        color=C_DARK, fontweight='bold')

# XGBoost block
xgb_box = FancyBboxPatch((4.0, 2.0), 4.0, 0.65, boxstyle="round,pad=0.08",
                          facecolor='#FFF3E0', edgecolor=C_STAGE2, linewidth=1.3, zorder=2)
ax.add_patch(xgb_box)
ax.text(6.0, 2.5, 'XGBoost Regressor', fontsize=7.5, ha='center', fontweight='bold', color=C_STAGE2, zorder=3)
ax.text(6.0, 2.22, '200 trees, max_depth=6, lr=0.05  |  Trained on 25,732 pairs, 1,901 proteins',
        fontsize=5, ha='center', color=C_GRAY, zorder=3)

arrow(ax, 6.0, 2.8, 6.0, 2.65, C_DARK, lw=1.0)

# argmax block
add_box(ax, 4.5, 1.15, 3.0, 0.65, 'Per-Protein Argmax',
        r'(L*, t*) = arg max  $\hat{y}$(L, t)  over all candidates', edgecolor=C_STAGE2, fontsize=6.5, subfontsize=5.2)

arrow(ax, 6.0, 2.0, 6.0, 1.8, C_DARK, lw=1.0)

# Output
out_box = FancyBboxPatch((3.5, 0.15), 5.0, 0.75, boxstyle="round,pad=0.1",
                          facecolor='#E8F5E9', edgecolor=C_OUTPUT, linewidth=1.8, zorder=2)
ax.add_patch(out_box)
ax.text(6.0, 0.7, 'Final Output', fontsize=8.5, ha='center', fontweight='bold', color='#1B5E20', zorder=3)
ax.text(6.0, 0.35, 'Best spliced alignment → exon F1  |  Closes 66–90% of paralog gap',
        fontsize=5.8, ha='center', color=C_DARK, zorder=3)

arrow(ax, 6.0, 1.15, 6.0, 0.9, C_OUTPUT, lw=1.3)

# ============================================================
# Left bracket: connects Stage 1 → Stage 2
# ============================================================
ax.plot([0.05, 0.05], [9.7, 5.1], color=C_STAGE1, linewidth=3, zorder=1)
ax.text(0.0, 7.4, 'Generate', fontsize=6.5, ha='center', va='center', color=C_STAGE1,
        rotation=90, fontweight='bold', zorder=3)

ax.plot([0.05, 0.05], [4.9, 1.5], color=C_STAGE2, linewidth=3, zorder=1)
ax.text(0.0, 3.2, 'Select', fontsize=6.5, ha='center', va='center', color=C_STAGE2,
        rotation=90, fontweight='bold', zorder=3)

# ============================================================
# Key stats box (top right)
# ============================================================
stats_box = FancyBboxPatch((8.5, 1.3), 3.0, 1.7, boxstyle="round,pad=0.08",
                            facecolor='#FFFDE7', edgecolor='#F9A825', linewidth=0.8, zorder=2)
ax.add_patch(stats_box)
ax.text(10.0, 2.85, 'Key Results', fontsize=6, ha='center', fontweight='bold', color='#F57F17', zorder=3)
stats_text = ('Router vs best heuristic:\n'
              '  +0.37 primary F1\n'
              'Router oracle efficiency: 98.7%\n'
              'Candidate pool: ~9–17 loci/protein\n'
              '24 features across 5 categories')
ax.text(10.0, 2.30, stats_text, fontsize=4.5, ha='center', va='center', color=C_DARK, zorder=3,
        linespacing=1.3)

# ============================================================
# Notes
# ============================================================
ax.text(0.4, 0.05, 'Phase 1–3 = Candidate Generation (Stage 1)    |    Phase 4 = Learned Joint Selection (Stage 2)',
        fontsize=5, ha='left', va='bottom', color=C_GRAY, style='italic')

fig.tight_layout(pad=0.2)
save_fig(fig, 'fig2_tprs_pipeline')
print('Fig 2 done.')
