"""Generate V4 new tables: heuristic baselines + correct-locus rate."""
import json
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SPECIES = ['mouse', 'dog', 'chicken', 'zebrafish']
SPECIES_LABELS = {'mouse': 'Mouse', 'dog': 'Dog', 'chicken': 'Chicken', 'zebrafish': 'Zebrafish'}

# Load heuristic baselines
with open(f'{ROOT}/results/tprs/phase5/heuristic_baselines.json') as f:
    hb = json.load(f)

# Load TPRS E1 (full pipeline)
with open(f'{ROOT}/results/tprs/phase5/ablation_summary_all_species.json') as f:
    abl = json.load(f)

# Load Exp 1 for correct_locus_rate
df_exp1 = pd.read_csv(f'{ROOT}/results/exp1/summary.tsv', sep='\t')

HEURISTIC_KEYS = ['always_primary', 'max_cstage', 'max_identity', 'max_tool_consensus', 'max_coverage']
HEURISTIC_LABELS = {
    'always_primary': 'Always-primary',
    'max_cstage': 'Max C-stage composite',
    'max_identity': 'Max identity',
    'max_tool_consensus': 'Max tool consensus',
    'max_coverage': 'Max coverage',
}

# ============================================================
# Table S4: Heuristic Baselines vs TPRS Router
# ============================================================
print("=" * 80)
print("TABLE S4: Heuristic Baselines vs TPRS Router")
print("=" * 80)

# Panel A: primary F1
print("\nPanel A: primary F1")
header = f"{'Heuristic':<28} " + " ".join(f"{SPECIES_LABELS[s]:>10}" for s in SPECIES) + f"  {'Mean':>8}"
print(header)
print("-" * len(header))

for hk in HEURISTIC_KEYS:
    vals = []
    ci_strs = []
    for sp in SPECIES:
        v = hb[sp][hk]['primary_f1']['mean']
        lo = hb[sp][hk]['primary_f1']['ci_lo']
        hi = hb[sp][hk]['primary_f1']['ci_hi']
        vals.append(v)
        ci_strs.append(f"{v:.3f} ({lo:.3f}, {hi:.3f})")
    mean_val = np.mean(vals)
    row = f"{HEURISTIC_LABELS[hk]:<28} " + " ".join(f"{v:>10.3f}" for v in vals) + f"  {mean_val:>8.3f}"
    print(row)

# TPRS row
tprs_vals = []
for sp in SPECIES:
    v = abl[sp]['primary']['E1']['mean']
    tprs_vals.append(v)
mean_tprs = np.mean(tprs_vals)
row = f"{'TPRS (XGBoost router)':<28} " + " ".join(f"{v:>10.3f}" for v in tprs_vals) + f"  {mean_tprs:>8.3f}"
print(row)

# Δ row vs best heuristic
print()
for sp_i, sp in enumerate(SPECIES):
    best_heur = max(hb[sp][hk]['primary_f1']['mean'] for hk in HEURISTIC_KEYS)
    delta = tprs_vals[sp_i] - best_heur
    print(f"  {SPECIES_LABELS[sp]}: TPRS {tprs_vals[sp_i]:.3f} - best heuristic {best_heur:.3f} = Δ +{delta:.3f}")

# Panel B: best_overlap F1
print("\nPanel B: best_overlap F1")
header = f"{'Heuristic':<28} " + " ".join(f"{SPECIES_LABELS[s]:>10}" for s in SPECIES) + f"  {'Mean':>8}"
print(header)
print("-" * len(header))

for hk in HEURISTIC_KEYS:
    vals = []
    for sp in SPECIES:
        v = hb[sp][hk]['best_overlap_f1']['mean']
        vals.append(v)
    mean_val = np.mean(vals)
    row = f"{HEURISTIC_LABELS[hk]:<28} " + " ".join(f"{v:>10.3f}" for v in vals) + f"  {mean_val:>8.3f}"
    print(row)

tprs_bo_vals = []
for sp in SPECIES:
    v = abl[sp]['best_overlap']['E1']['mean']
    tprs_bo_vals.append(v)
mean_tprs_bo = np.mean(tprs_bo_vals)
row = f"{'TPRS (XGBoost router)':<28} " + " ".join(f"{v:>10.3f}" for v in tprs_bo_vals) + f"  {mean_tprs_bo:>8.3f}"
print(row)

# ============================================================
# Table S6: Correct-Locus Rate Collapse
# ============================================================
print("\n" + "=" * 80)
print("TABLE S6: Correct-Locus Rate — 1:1 vs Paralog")
print("=" * 80)

TOOLS_ORDER = ['miniprot', 'Spaln3', 'GTH', 'Scipio', 'genBlastG', 'mini+miniSplice']
TOOL_KEYS = ['miniprot', 'spaln3', 'gth', 'scipio', 'genblastg', 'minisplice_miniprot']

print(f"\n{'Tool':<22}", end="")
for sp in SPECIES:
    print(f"  {'':>6} {SPECIES_LABELS[sp]:<18}", end="")
print(f"  {'Mean':>8}")
print(f"{'':22}", end="")
for sp in SPECIES:
    print(f"  {'1:1':>6} {'Paralog':>6} {'Drop':>6}", end="")
print(f"  {'Drop':>8}")
print("-" * 115)

for ti, (tool_label, tool_key) in enumerate(zip(TOOLS_ORDER, TOOL_KEYS)):
    print(f"{tool_label:<22}", end="")
    drops = []
    for sp in SPECIES:
        r1 = df_exp1[(df_exp1['species']==sp) & (df_exp1['tool']==tool_key) & (df_exp1['regime']=='one_to_one')]
        rp = df_exp1[(df_exp1['species']==sp) & (df_exp1['tool']==tool_key) & (df_exp1['regime']=='paralog')]
        clr_1 = r1['correct_locus_rate'].values[0] if len(r1)>0 else 0
        clr_p = rp['correct_locus_rate'].values[0] if len(rp)>0 else 0
        drop = clr_1 - clr_p
        drops.append(drop)
        print(f"  {clr_1:>6.3f} {clr_p:>6.3f} {drop:>6.3f}", end="")
    mean_drop = np.mean(drops)
    print(f"  {mean_drop:>8.3f}")

print("\nDone.")
