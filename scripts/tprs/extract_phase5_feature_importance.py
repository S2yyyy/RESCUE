"""Extract Phase 5 XGBoost feature importance for all 4 species."""
import sys
import numpy as np
import pandas as pd
import json
from pathlib import Path
from xgboost import XGBRegressor

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import PROJECT_ROOT  # noqa: E402

OUTDIR = PROJECT_ROOT / "results/tprs/phase5"
SEED = 4242

FEATURE_COLS = [
    "a_diamond", "a_last", "a_miniprot", "a_mmseqs",
    "a_n_tools", "a_best_score", "a_best_qcov", "a_n_raw_hits",
    "region_size", "region_size_log",
    "on_primary_chrom", "strand_match_primary",
    "coorth_overlap_bp", "coorth_same_chrom",
    "query_length_aa", "query_length_log",
    "c_has_hit", "c_identity", "c_coverage", "c_n_exons",
    "c_cds_length", "c_composite_score", "c_hit_score", "c_rank",
    "query_type_num",
]

FEAT_SHORT = {
    "a_diamond": "A | diamond present", "a_last": "A | last present",
    "a_miniprot": "A | miniprot present", "a_mmseqs": "A | mmseqs present",
    "a_n_tools": "A | num tools hitting", "a_best_score": "A | best bitscore",
    "a_best_qcov": "A | best query cov.", "a_n_raw_hits": "A | num raw hits",
    "region_size": "G | region size", "region_size_log": "G | region size (log10)",
    "on_primary_chrom": "G | on primary chrom.", "strand_match_primary": "G | strand matches orth.",
    "coorth_overlap_bp": "G | co-orth. overlap", "coorth_same_chrom": "G | co-orth. same chrom.",
    "query_length_aa": "Q | protein length", "query_length_log": "Q | protein length (log10)",
    "c_has_hit": "C | has spliced hit", "c_identity": "C | alignment identity",
    "c_coverage": "C | query coverage", "c_n_exons": "C | num exons",
    "c_cds_length": "C | CDS length", "c_composite_score": "C | composite score",
    "c_hit_score": "C | miniprot raw score", "c_rank": "C | hit rank",
    "query_type_num": "S | query type (human / sp)",
}

FEAT_CATEGORY = {}
for f in FEATURE_COLS:
    if f.startswith("a_"):
        FEAT_CATEGORY[f] = "A-stage"
    elif f.startswith("c_"):
        FEAT_CATEGORY[f] = "C-stage"
    elif f in ("region_size", "region_size_log", "on_primary_chrom",
               "strand_match_primary", "coorth_overlap_bp", "coorth_same_chrom"):
        FEAT_CATEGORY[f] = "Geometry"
    elif f in ("query_length_aa", "query_length_log"):
        FEAT_CATEGORY[f] = "Query"
    elif f == "query_type_num":
        FEAT_CATEGORY[f] = "Species ind."

data_paths = {
    "mouse": OUTDIR / "training_data_mouse.tsv",
    "dog": OUTDIR / "dog" / "training_data.tsv",
    "chicken": OUTDIR / "chicken" / "training_data.tsv",
    "zebrafish": OUTDIR / "zebrafish" / "training_data.tsv",
}

results = {}
for sp, path in data_paths.items():
    if not path.exists():
        print(f"MISSING: {path}")
        continue
    df = pd.read_csv(path, sep="\t")
    if "query_type_num" not in df.columns:
        df["query_type_num"] = (df["query_type"] != "human").astype(int)

    X = df[FEATURE_COLS].fillna(0.0)
    y = df["label_best_overlap_f1"].fillna(0.0)

    model = XGBRegressor(
        n_estimators=200, max_depth=6, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        random_state=SEED, n_jobs=-1,
    )
    model.fit(X, y)
    imp = dict(zip(FEATURE_COLS, model.feature_importances_.tolist()))
    results[sp] = imp
    print(f"{sp}: {len(df)} rows, {df['query_id'].nunique()} queries, gain_sum={sum(imp.values()):.4f}")

# Mean across species
mean_imp = {f: float(np.mean([results[sp][f] for sp in results])) for f in FEATURE_COLS}
sorted_feats = sorted(mean_imp.items(), key=lambda x: -x[1])

print(f"\n{'Rank':<5s} {'Feature':<35s} {'Category':<13s} {'Mean':>8s}  "
      + "".join(f"{s:>9s}" for s in ["Mouse", "Dog", "Chicken", "Zfish"]))
for rank, (feat, mean_g) in enumerate(sorted_feats, 1):
    short = FEAT_SHORT.get(feat, feat)
    cat = FEAT_CATEGORY.get(feat, "?")
    vals = "  ".join(f"{results[sp][feat]:>8.4f}" for sp in ["mouse","dog","chicken","zebrafish"])
    print(f"{rank:<5d} {short:<35s} {cat:<13s} {mean_g:>8.4f}  {vals}")

json_path = OUTDIR / "feature_importance_phase5.json"
with open(json_path, "w") as f:
    json.dump({
        "mean": mean_imp,
        "per_species": {sp: results[sp] for sp in results},
        "sorted": [(f, mean_imp[f], FEAT_CATEGORY.get(f,""), FEAT_SHORT.get(f,f))
                   for f, _ in sorted_feats],
    }, f, indent=2)
print(f"\nSaved: {json_path}")
