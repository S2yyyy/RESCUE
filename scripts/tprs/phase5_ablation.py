"""Phase 5 Ablation: Multi-species pipeline component ablation (E1-E6).

Evaluates both best_overlap F1 and primary F1 for each experiment.

E1: Full Phase 5 (baseline)
E2: Remove species protein (human-only query)
E3: Remove Router (heuristic always-primary)
E4: Remove Phase 3 co-ortholog recovery
E5: Remove C-stage features (A-stage only)
E6: Remove individual A-stage tools (re-merge BED files with 3 tools)

Usage:
    python scripts/tprs/phase5_ablation.py --species mouse dog chicken zebrafish
    python scripts/tprs/phase5_ablation.py --all

Output (per species):
    results/tprs/phase5/{species}/ablation_summary.json
    results/tprs/phase5/{species}/ablation_details.json
    Stdout: formatted ablation tables (per-species + combined)

Spec anchor: docs/phase5_ablation_plan.md
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit
from xgboost import XGBRegressor

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
RANDOM_SEED = 4242
N_BOOTSTRAP = 10_000

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import PROJECT_ROOT  # noqa: E402

REPO = PROJECT_ROOT
BED_DIR = REPO / "results/tprs/stage0_bed"
OUTDIR_BASE = REPO / "results/tprs/phase5"

FEATURE_COLS_ALL = [
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

FEATURE_COLS_A_ONLY = [col for col in FEATURE_COLS_ALL if not col.startswith("c_")]

TOOL_NAMES = ["diamond", "last", "miniprot", "mmseqs"]

XGB_PARAMS = dict(
    n_estimators=200, max_depth=6, learning_rate=0.05,
    subsample=0.8, colsample_bytree=0.8,
    random_state=RANDOM_SEED, n_jobs=-1,
)

LABEL_CONFIGS = [
    ("label_best_overlap_f1", "best_overlap"),
    ("label_primary_f1", "primary"),
]

# Per-species E1 verification targets (best_overlap, approximate)
VERIFY_TARGETS = {
    "mouse": 0.836,
    "dog": 0.8535,
    "chicken": 0.8114,
    "zebrafish": 0.6498,
}

# Data paths: mouse uses legacy flat path; other species use subdirectory
def get_data_path(sp: str) -> Path:
    if sp == "mouse":
        return OUTDIR_BASE / "training_data_mouse.tsv"
    return OUTDIR_BASE / sp / "training_data.tsv"

def get_outdir(sp: str) -> Path:
    return OUTDIR_BASE / sp

def get_bed_suffix(sp: str) -> str:
    return f"relaxed.{sp}.paralog.bed"


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def paired_bootstrap_ci(values: np.ndarray, n_bootstrap: int = N_BOOTSTRAP,
                        ci: float = 0.95, seed: int = RANDOM_SEED) -> tuple[float, float, float]:
    rng = np.random.RandomState(seed)
    n = len(values)
    means = np.array([rng.choice(values, size=n, replace=True).mean()
                      for _ in range(n_bootstrap)])
    alpha = (1.0 - ci) / 2.0
    lower, upper = np.percentile(means, [100 * alpha, 100 * (1.0 - alpha)])
    return float(values.mean()), float(lower), float(upper)


def train_xgboost(X_train: pd.DataFrame, y_train: pd.Series) -> XGBRegressor:
    model = XGBRegressor(**XGB_PARAMS)
    model.fit(X_train.fillna(0.0), y_train.fillna(0.0))
    return model


def evaluate_per_query(test_df: pd.DataFrame, model, feature_cols: list[str],
                       label_col: str) -> np.ndarray:
    queries = test_df["query_id"].unique()
    f1s = np.zeros(len(queries))
    for i, qid in enumerate(queries):
        qdf = test_df[test_df["query_id"] == qid]
        if qdf.empty:
            continue
        X = qdf[feature_cols].fillna(0.0)
        preds = model.predict(X)
        best_row = qdf.iloc[np.argmax(preds)]
        f1s[i] = best_row[label_col]
    return f1s


def compute_metrics(f1s: np.ndarray, label: str, metric_name: str) -> dict:
    mean_val, lo, hi = paired_bootstrap_ci(f1s)
    return {
        "experiment": label,
        "metric": metric_name,
        "n_queries": int(len(f1s)),
        "mean": mean_val,
        "ci_lower": lo,
        "ci_upper": hi,
        "zero_f1_count": int((f1s == 0).sum()),
    }


def make_per_query_detail(test_df: pd.DataFrame, model,
                          feature_cols: list[str], label_col: str,
                          selection_method: str = "argmax_pred") -> list[dict]:
    queries = test_df["query_id"].unique()
    details = []
    for qid in queries:
        qdf = test_df[test_df["query_id"] == qid]
        if qdf.empty:
            details.append({"query_id": qid, "f1": 0.0, "region_idx": -1,
                            "selection_method": selection_method})
            continue
        if selection_method == "argmax_pred" and model is not None:
            X = qdf[feature_cols].fillna(0.0)
            preds = model.predict(X)
            best_idx = int(np.argmax(preds))
        elif selection_method == "heuristic":
            best_idx = _heuristic_select(qdf)
        else:
            best_idx = int(qdf[label_col].idxmax())
        best_row = qdf.iloc[best_idx]
        details.append({
            "query_id": qid,
            "region_idx": int(best_row["region_idx"]),
            "query_type": str(best_row["query_type"]),
            "chrom": str(best_row["chrom"]),
            "on_primary_chrom": int(best_row["on_primary_chrom"]),
            "coorth_overlap_bp": int(best_row["coorth_overlap_bp"]),
            "f1": float(best_row[label_col]),
            "label_primary_f1": float(best_row["label_primary_f1"]),
            "label_best_overlap_f1": float(best_row["label_best_overlap_f1"]),
            "selection_method": selection_method,
        })
    return details


def _heuristic_select(qdf: pd.DataFrame) -> int:
    """E3 heuristic: human rows, on_primary_chrom argmax, tiebreak coorth_overlap_bp."""
    subset = qdf[qdf["query_type"] == "human"]
    if subset.empty:
        subset = qdf
    pos = int(subset["on_primary_chrom"].values.argmax())
    label = subset.index[pos]
    if (subset["on_primary_chrom"].values == 1).sum() > 1:
        primary_sub = subset[subset["on_primary_chrom"] == 1]
        pos = int(primary_sub["coorth_overlap_bp"].values.argmax())
        label = primary_sub.index[pos]
    return qdf.index.get_loc(label)


# ---------------------------------------------------------------------------
# BED merging (for E6: 3-tool candidate region reconstruction)
# ---------------------------------------------------------------------------

def load_bed(path: Path) -> dict[str, list[tuple]]:
    """Return {query_id: [(chrom, start, end, strand, score, q_cov), ...]}."""
    out: dict[str, list[tuple]] = defaultdict(list)
    if not path.exists():
        return out
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.strip().split("\t")
            q_cov = float(f[8]) if len(f) > 8 else 0.0
            out[f[3]].append((f[0], int(f[1]), int(f[2]), f[5],
                             float(f[4]) if len(f) > 4 else 0.0, q_cov))
    return out


def merge_regions(candidates: list[tuple], slop: int = 50000) -> list[dict]:
    """Merge overlapping candidate regions per chrom+strand with slop."""
    if not candidates:
        return []
    by_key: dict[tuple, list[tuple]] = defaultdict(list)
    for c in candidates:
        by_key[(c[0], c[3])].append(c)

    merged = []
    for (chrom, strand), regs in by_key.items():
        regs_sorted = sorted(regs, key=lambda x: x[1])
        group = [regs_sorted[0]]
        for r in regs_sorted[1:]:
            cur_end = max(x[2] for x in group)
            if r[1] <= cur_end + slop:
                group.append(r)
            else:
                merged.append(_collapse_group(chrom, strand, group, slop))
                group = [r]
        if group:
            merged.append(_collapse_group(chrom, strand, group, slop))
    return merged


def _collapse_group(chrom: str, strand: str, group: list[tuple],
                    slop: int) -> dict:
    start = max(0, min(x[1] for x in group) - slop)
    end = max(x[2] for x in group) + slop
    tools = set(x[6] for x in group if len(x) > 6)
    scores = [x[4] for x in group if x[4] > 0]
    qcovs = [x[5] for x in group if x[5] > 0]
    return {
        "chrom": chrom, "start": start, "end": end, "strand": strand,
        "n_tools": len(tools),
        "tools": tools,
        "best_score": max(scores) if scores else 0.0,
        "best_qcov": max(qcovs) if qcovs else 0.0,
        "score_sum": sum(scores),
    }


def _safe_get(row, col: str, default):
    """Get column value with fallback for columns that may not exist (status, wall_s)."""
    if col in row.index:
        return row[col]
    return default


def build_3tool_dataset(df: pd.DataFrame, removed_tool: str,
                        bed_suffix: str) -> pd.DataFrame:
    """Rebuild the training dataset using only 3 out of 4 A-stage tools.

    Steps:
    1. Load BED files for the 3 remaining tools
    2. Merge regions per query (50kb slop)
    3. Match each merged region to the closest row in the existing training data
       (by query_id + chrom + reciprocal overlap > 50%)
    4. Inherit C-stage features, labels, and metadata from matched rows
    5. Compute A-stage features from the 3-tool merge

    Returns a DataFrame with the same columns as the full training data.
    """
    tool_names = [t for t in TOOL_NAMES if t != removed_tool]

    all_cands: dict[str, dict[str, list[tuple]]] = {}
    for tool in tool_names:
        bed_path = BED_DIR / f"{tool}_{bed_suffix}"
        all_cands[tool] = load_bed(bed_path)

    all_qids = set()
    for tool in tool_names:
        all_qids.update(all_cands[tool].keys())

    lookup: dict[tuple, list[tuple[int, int, int]]] = defaultdict(list)
    for idx, row in df.iterrows():
        lookup[(row["query_id"], str(row["chrom"]))].append(
            (int(row["start"]), int(row["end"]), idx))

    query_primary_info = {}
    for qid in all_qids:
        qdf = df[df["query_id"] == qid]
        if not qdf.empty:
            prim_rows = qdf[qdf["on_primary_chrom"] == 1]
            if not prim_rows.empty:
                query_primary_info[qid] = {
                    "primary_chrom": str(prim_rows.iloc[0]["chrom"]),
                    "primary_strand": str(prim_rows.iloc[0]["strand"]),
                    "coorth_overlap_bp": float(prim_rows["coorth_overlap_bp"].max()),
                    "coorth_same_chrom": int(prim_rows["coorth_same_chrom"].max()),
                    "query_length_aa": float(qdf.iloc[0]["query_length_aa"]),
                    "query_length_log": float(qdf.iloc[0]["query_length_log"]),
                }
            else:
                query_primary_info[qid] = {
                    "primary_chrom": str(qdf.iloc[0]["chrom"]),
                    "primary_strand": str(qdf.iloc[0]["strand"]),
                    "coorth_overlap_bp": 0.0,
                    "coorth_same_chrom": 0,
                    "query_length_aa": float(qdf.iloc[0]["query_length_aa"]),
                    "query_length_log": float(qdf.iloc[0]["query_length_log"]),
                }

    rows = []
    matched_count = 0
    unmatched_count = 0

    for qid in sorted(all_qids):
        raw = []
        for tool in tool_names:
            for c in all_cands.get(tool, {}).get(qid, []):
                raw.append((*c[:6], tool))

        if not raw:
            continue

        merged = merge_regions(raw, slop=50000)

        pinfo = query_primary_info.get(qid, {})
        primary_chrom = pinfo.get("primary_chrom", "")
        primary_strand = pinfo.get("primary_strand", "")
        coorth_overlap = pinfo.get("coorth_overlap_bp", 0.0)
        coorth_same = pinfo.get("coorth_same_chrom", 0)
        qlen_aa = pinfo.get("query_length_aa", 0.0)
        qlen_log = pinfo.get("query_length_log", 0.0)

        for ri, region in enumerate(merged):
            rchrom = region["chrom"]
            rstart = region["start"]
            rend = region["end"]

            candidates = lookup.get((qid, rchrom), [])
            best_overlap = 0.0
            best_row = None
            for cstart, cend, cidx in candidates:
                overlap = max(0, min(rend, cend) - max(rstart, cstart))
                region_len = rend - rstart
                cand_len = cend - cstart
                if region_len > 0 and cand_len > 0:
                    overlap_frac = overlap / max(region_len, cand_len)
                else:
                    overlap_frac = 0.0
                if overlap_frac > best_overlap and overlap_frac > 0.5:
                    best_overlap = overlap_frac
                    best_row = df.iloc[cidx]

            if best_row is not None:
                matched_count += 1
                row_dict = {
                    "query_id": qid,
                    "region_idx": ri,
                    "chrom": rchrom,
                    "start": rstart,
                    "end": rend,
                    "strand": region["strand"],
                    # A-stage: only 3 tools
                    "a_diamond": int("diamond" in region["tools"]),
                    "a_last": int("last" in region["tools"]),
                    "a_miniprot": int("miniprot" in region["tools"]),
                    "a_mmseqs": int("mmseqs" in region["tools"]),
                    "a_n_tools": region["n_tools"],
                    "a_best_score": region["best_score"],
                    "a_best_qcov": region["best_qcov"],
                    "a_n_raw_hits": len([t for t in tool_names if t in region["tools"]]),
                    "region_size": rend - rstart,
                    "region_size_log": np.log(max(1, rend - rstart)),
                    "on_primary_chrom": int(rchrom == primary_chrom),
                    "strand_match_primary": int(region["strand"] == primary_strand),
                    "coorth_overlap_bp": coorth_overlap,
                    "coorth_same_chrom": coorth_same,
                    "query_length_aa": qlen_aa,
                    "query_length_log": qlen_log,
                    # C-stage: inherited from matched row
                    "c_has_hit": best_row["c_has_hit"],
                    "c_identity": best_row["c_identity"],
                    "c_coverage": best_row["c_coverage"],
                    "c_n_exons": best_row["c_n_exons"],
                    "c_cds_length": best_row["c_cds_length"],
                    "c_composite_score": best_row["c_composite_score"],
                    "c_hit_score": best_row["c_hit_score"],
                    "c_rank": best_row["c_rank"],
                    # Labels: inherited from matched row
                    "label_primary_f1": best_row["label_primary_f1"],
                    "label_best_overlap_f1": best_row["label_best_overlap_f1"],
                    "label_primary_precision": best_row["label_primary_precision"],
                    "label_primary_recall": best_row["label_primary_recall"],
                    "label_best_gene": best_row["label_best_gene"],
                    "label_n_pred_exons": best_row["label_n_pred_exons"],
                    "label_pred_cds_bp": best_row["label_pred_cds_bp"],
                    "query_type": best_row["query_type"],
                }
                # Optional columns (mouse has status/wall_s; dog/chicken/zebrafish don't)
                if "status" in best_row.index:
                    row_dict["status"] = best_row["status"]
                if "wall_s" in best_row.index:
                    row_dict["wall_s"] = best_row["wall_s"]
            else:
                unmatched_count += 1
                row_dict = {
                    "query_id": qid, "region_idx": ri,
                    "chrom": rchrom, "start": rstart, "end": rend,
                    "strand": region["strand"],
                    "a_diamond": int("diamond" in region["tools"]),
                    "a_last": int("last" in region["tools"]),
                    "a_miniprot": int("miniprot" in region["tools"]),
                    "a_mmseqs": int("mmseqs" in region["tools"]),
                    "a_n_tools": region["n_tools"],
                    "a_best_score": region["best_score"],
                    "a_best_qcov": region["best_qcov"],
                    "a_n_raw_hits": len([t for t in tool_names if t in region["tools"]]),
                    "region_size": rend - rstart,
                    "region_size_log": np.log(max(1, rend - rstart)),
                    "on_primary_chrom": int(rchrom == primary_chrom),
                    "strand_match_primary": int(region["strand"] == primary_strand),
                    "coorth_overlap_bp": coorth_overlap,
                    "coorth_same_chrom": coorth_same,
                    "query_length_aa": qlen_aa,
                    "query_length_log": qlen_log,
                    "c_has_hit": 0.0, "c_identity": 0.0, "c_coverage": 0.0,
                    "c_n_exons": 0.0, "c_cds_length": 0.0,
                    "c_composite_score": 0.0, "c_hit_score": 0.0, "c_rank": 0.0,
                    "label_primary_f1": 0.0, "label_best_overlap_f1": 0.0,
                    "label_primary_precision": 0.0, "label_primary_recall": 0.0,
                    "label_best_gene": "", "label_n_pred_exons": 0.0,
                    "label_pred_cds_bp": 0.0,
                    "query_type": "human",
                }
            rows.append(row_dict)

    new_df = pd.DataFrame(rows)
    new_df["query_type_num"] = (new_df["query_type"] == "mouse").astype(int)
    print(f"    Re-merged: {len(new_df)} rows, {new_df['query_id'].nunique()} queries "
          f"(matched={matched_count}, new={unmatched_count})")
    return new_df


# ---------------------------------------------------------------------------
# E1: Full Phase 5
# ---------------------------------------------------------------------------

def run_e1(df: pd.DataFrame, label_col: str, metric_name: str) -> tuple[dict, list[dict], tuple]:
    df = df.copy()
    df["query_type_num"] = (df["query_type"] == "mouse").astype(int)

    gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=RANDOM_SEED)
    train_idx, test_idx = next(gss.split(df, groups=df["query_id"]))
    train_df = df.iloc[train_idx].copy()
    test_df = df.iloc[test_idx].copy()

    train_queries = set(train_df["query_id"].unique())
    test_queries = set(test_df["query_id"].unique())
    assert len(train_queries & test_queries) == 0, "Query leakage in E1!"

    X_train = train_df[FEATURE_COLS_ALL].fillna(0.0)
    y_train = train_df[label_col].fillna(0.0)
    model = train_xgboost(X_train, y_train)

    f1s = evaluate_per_query(test_df, model, FEATURE_COLS_ALL, label_col)
    metrics = compute_metrics(f1s, "E1: Full Phase 5", metric_name)
    details = make_per_query_detail(test_df, model, FEATURE_COLS_ALL,
                                    label_col, selection_method="argmax_pred")

    print(f"  E1: {metric_name}={metrics['mean']:.4f} "
          f"[{metrics['ci_lower']:.4f}, {metrics['ci_upper']:.4f}]")
    return metrics, details, (train_queries, test_queries)


# ---------------------------------------------------------------------------
# E2: Remove Species Protein (human-only)
# ---------------------------------------------------------------------------

def run_e2(df: pd.DataFrame, e1_split: tuple, label_col: str,
           metric_name: str) -> tuple[dict, list[dict]]:
    train_queries, test_queries = e1_split
    df2 = df[df["query_type"] == "human"].copy()
    feature_cols = [c for c in FEATURE_COLS_ALL if c != "query_type_num"]

    train_df = df2[df2["query_id"].isin(train_queries)].copy()
    test_df = df2[df2["query_id"].isin(test_queries)].copy()

    X_train = train_df[feature_cols].fillna(0.0)
    y_train = train_df[label_col].fillna(0.0)
    model = train_xgboost(X_train, y_train)

    f1s = evaluate_per_query(test_df, model, feature_cols, label_col)
    metrics = compute_metrics(f1s, "E2: - Species Protein", metric_name)
    details = make_per_query_detail(test_df, model, feature_cols,
                                    label_col, selection_method="argmax_pred")
    print(f"  E2: {metric_name}={metrics['mean']:.4f} "
          f"[{metrics['ci_lower']:.4f}, {metrics['ci_upper']:.4f}]")
    return metrics, details


# ---------------------------------------------------------------------------
# E3: Remove Router (heuristic)
# ---------------------------------------------------------------------------

def run_e3(df: pd.DataFrame, test_queries: set, label_col: str,
           metric_name: str) -> tuple[dict, list[dict]]:
    test_df = df[df["query_id"].isin(test_queries)].copy()
    f1s = []
    details = []
    for qid in test_queries:
        qdf = test_df[test_df["query_id"] == qid]
        if qdf.empty:
            f1s.append(0.0)
            details.append({"query_id": qid, "f1": 0.0, "region_idx": -1,
                            "selection_method": "heuristic"})
            continue
        best_idx = _heuristic_select(qdf)
        best_row = qdf.iloc[best_idx]
        f1s.append(best_row[label_col])
        details.append({
            "query_id": qid,
            "region_idx": int(best_row["region_idx"]),
            "query_type": str(best_row["query_type"]),
            "chrom": str(best_row["chrom"]),
            "on_primary_chrom": int(best_row["on_primary_chrom"]),
            "coorth_overlap_bp": int(best_row["coorth_overlap_bp"]),
            "f1": float(best_row[label_col]),
            "label_primary_f1": float(best_row["label_primary_f1"]),
            "label_best_overlap_f1": float(best_row["label_best_overlap_f1"]),
            "selection_method": "heuristic",
        })
    f1s_arr = np.array(f1s)
    metrics = compute_metrics(f1s_arr, "E3: - Router (heuristic)", metric_name)
    print(f"  E3: {metric_name}={metrics['mean']:.4f} "
          f"[{metrics['ci_lower']:.4f}, {metrics['ci_upper']:.4f}]")
    return metrics, details


# ---------------------------------------------------------------------------
# E4: Remove Phase 3 (co-ortholog)
# ---------------------------------------------------------------------------

def run_e4(df: pd.DataFrame, e1_split: tuple, label_col: str,
           metric_name: str) -> tuple[dict, list[dict]]:
    train_queries, test_queries = e1_split

    df4 = df.copy()
    coorth_mask = df4["coorth_overlap_bp"] > 0
    regions_removed = int(coorth_mask.sum())
    affected_queries = set(df4[coorth_mask]["query_id"].unique())
    before_with_regions = set(df4[df4["query_id"].isin(test_queries)]["query_id"].unique())

    df4 = df4[~coorth_mask].copy()
    df4["query_type_num"] = (df4["query_type"] == "mouse").astype(int)

    after_test = set(df4[df4["query_id"].isin(test_queries)]["query_id"].unique())
    queries_without_regions = before_with_regions - after_test

    train_df = df4[df4["query_id"].isin(train_queries)].copy()
    test_df = df4[df4["query_id"].isin(test_queries)].copy()

    X_train = train_df[FEATURE_COLS_ALL].fillna(0.0)
    y_train = train_df[label_col].fillna(0.0)
    model = train_xgboost(X_train, y_train)

    f1s = evaluate_per_query(test_df, model, FEATURE_COLS_ALL, label_col)
    missing = set(test_queries) - set(test_df["query_id"].unique())
    f1s_ext = np.concatenate([f1s, np.zeros(len(missing))])

    metrics = compute_metrics(f1s_ext, "E4: - Phase 3 (co-ortholog)", metric_name)
    metrics["regions_removed"] = regions_removed
    metrics["affected_queries"] = len(affected_queries)
    metrics["queries_without_regions"] = len(queries_without_regions)

    affected_only_f1s = []
    for qid in affected_queries & test_queries:
        qdf = test_df[test_df["query_id"] == qid]
        if qdf.empty:
            affected_only_f1s.append(0.0)
        else:
            X = qdf[FEATURE_COLS_ALL].fillna(0.0)
            preds = model.predict(X)
            best_row = qdf.iloc[np.argmax(preds)]
            affected_only_f1s.append(best_row[label_col])
    if affected_only_f1s:
        affected_arr = np.array(affected_only_f1s)
        a_mean, a_lo, a_hi = paired_bootstrap_ci(affected_arr)
        metrics["affected_subset_mean"] = a_mean
        metrics["affected_subset_ci_lower"] = a_lo
        metrics["affected_subset_ci_upper"] = a_hi
        metrics["affected_subset_n"] = len(affected_arr)

    details = make_per_query_detail(test_df, model, FEATURE_COLS_ALL,
                                    label_col, selection_method="argmax_pred")
    print(f"  E4: {metric_name}={metrics['mean']:.4f} "
          f"[{metrics['ci_lower']:.4f}, {metrics['ci_upper']:.4f}] "
          f"(removed {regions_removed} rows, {len(affected_queries)} queries, "
          f"{len(queries_without_regions)} lost all regions)")
    return metrics, details


# ---------------------------------------------------------------------------
# E5: Remove C-stage Features
# ---------------------------------------------------------------------------

def run_e5(df: pd.DataFrame, e1_split: tuple, label_col: str,
           metric_name: str) -> tuple[dict, list[dict]]:
    train_queries, test_queries = e1_split
    df5 = df.copy()
    df5["query_type_num"] = (df5["query_type"] == "mouse").astype(int)

    train_df = df5[df5["query_id"].isin(train_queries)].copy()
    test_df = df5[df5["query_id"].isin(test_queries)].copy()

    feature_cols = FEATURE_COLS_A_ONLY
    X_train = train_df[feature_cols].fillna(0.0)
    y_train = train_df[label_col].fillna(0.0)
    model = train_xgboost(X_train, y_train)

    f1s = evaluate_per_query(test_df, model, feature_cols, label_col)
    metrics = compute_metrics(f1s, "E5: - C-stage Features", metric_name)
    details = make_per_query_detail(test_df, model, feature_cols,
                                    label_col, selection_method="argmax_pred")
    print(f"  E5: {metric_name}={metrics['mean']:.4f} "
          f"[{metrics['ci_lower']:.4f}, {metrics['ci_upper']:.4f}]")
    return metrics, details


# ---------------------------------------------------------------------------
# E6: Remove Individual A-stage Tools (BED re-merge)
# ---------------------------------------------------------------------------

def run_e6(df_full: pd.DataFrame, e1_split: tuple, label_col: str,
           metric_name: str, bed_suffix: str) -> tuple[list[dict], list[list[dict]]]:
    """E6: for each tool, rebuild the dataset using only the other 3 tools' BED files."""
    train_queries, test_queries = e1_split
    all_metrics = []
    all_details = []

    for tool in TOOL_NAMES:
        print(f"  E6 ({tool}): re-merging BED files with 3 tools (no {tool})...")
        new_df = build_3tool_dataset(df_full, tool, bed_suffix)

        train_df = new_df[new_df["query_id"].isin(train_queries)].copy()
        test_df = new_df[new_df["query_id"].isin(test_queries)].copy()

        X_train = train_df[FEATURE_COLS_ALL].fillna(0.0)
        y_train = train_df[label_col].fillna(0.0)
        model = train_xgboost(X_train, y_train)

        f1s = evaluate_per_query(test_df, model, FEATURE_COLS_ALL, label_col)
        missing = set(test_queries) - set(test_df["query_id"].unique())
        f1s_ext = np.concatenate([f1s, np.zeros(len(missing))])

        metrics = compute_metrics(f1s_ext, f"E6: - {tool}", metric_name)
        metrics["tool"] = tool
        metrics["queries_lost_all_regions"] = len(missing)
        all_metrics.append(metrics)

        details = make_per_query_detail(test_df, model, FEATURE_COLS_ALL,
                                        label_col, selection_method="argmax_pred")
        all_details.append(details)

        print(f"  E6 ({tool}): {metric_name}={metrics['mean']:.4f} "
              f"[{metrics['ci_lower']:.4f}, {metrics['ci_upper']:.4f}] "
              f"({len(missing)} queries lost all regions)")
    return all_metrics, all_details


# ---------------------------------------------------------------------------
# Per-species ablation runner
# ---------------------------------------------------------------------------

def run_ablation_for_species(sp: str) -> dict:
    """Run E1-E6 ablation for a single species. Returns results dict."""
    data_path = get_data_path(sp)
    outdir = get_outdir(sp)
    bed_suffix = get_bed_suffix(sp)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'=' * 70}")
    print(f"  Phase 5 Ablation: {sp.upper()}")
    print(f"  Data: {data_path}")
    print(f"  BED suffix: {bed_suffix}")
    print(f"  Output: {outdir}")
    print(f"{'=' * 70}")

    df = pd.read_csv(data_path, sep="\t")
    print(f"Data: {len(df)} rows, {df['query_id'].nunique()} queries")

    all_results = {}
    details_output = {}

    for label_col, metric_name in LABEL_CONFIGS:
        print(f"\n  --- {metric_name.upper()} ---")

        # E1
        e1_metrics, e1_details, e1_split = run_e1(df, label_col, metric_name)
        train_queries, test_queries = e1_split

        # E2
        e2_metrics, e2_details = run_e2(df, e1_split, label_col, metric_name)
        # E3
        e3_metrics, e3_details = run_e3(df, test_queries, label_col, metric_name)
        # E4
        e4_metrics, e4_details = run_e4(df, e1_split, label_col, metric_name)
        # E5
        e5_metrics, e5_details = run_e5(df, e1_split, label_col, metric_name)
        # E6
        e6_metrics_list, e6_details_list = run_e6(df, e1_split, label_col,
                                                   metric_name, bed_suffix)

        key = metric_name
        all_results[key] = {
            "E1": e1_metrics, "E2": e2_metrics, "E3": e3_metrics,
            "E4": e4_metrics, "E5": e5_metrics,
        }
        for i, m in enumerate(e6_metrics_list):
            all_results[key][f"E6_{m['tool']}"] = m

        details_output[key] = {
            "E1": e1_details, "E2": e2_details, "E3": e3_details,
            "E4": e4_details, "E5": e5_details,
        }
        for i, d in enumerate(e6_details_list):
            details_output[key][f"E6_{e6_metrics_list[i]['tool']}"] = d

    # ---- Verification ----
    expected = VERIFY_TARGETS.get(sp)
    if expected is not None:
        e1_best = all_results["best_overlap"]["E1"]["mean"]
        diff = abs(e1_best - expected)
        status = "PASS" if diff < 0.05 else "WARN"
        print(f"\n  Verification: E1 best_overlap={e1_best:.4f} "
              f"(expected≈{expected:.4f}, diff={diff:.4f}) → {status}")

    # ---- Print table ----
    for _, metric_name in LABEL_CONFIGS:
        results = all_results[metric_name]
        e1_f1 = results["E1"]["mean"]

        print(f"\n  {'─' * 70}")
        print(f"  Pipeline Component Ablation — {sp.upper()} {metric_name.upper()}")
        print(f"  {'─' * 70}")
        print(f"  {'Experiment':<40s} {'F1':>19s} {'Delta':>8s}")
        print(f"  {'─' * 60}")

        print(f"  {'E1: Full Phase 5 pipeline':<40s} "
              f"{results['E1']['mean']:>4.2f} "
              f"[{results['E1']['ci_lower']:.2f},{results['E1']['ci_upper']:.2f}]  —")

        for label, key in [
            ("E2: - Species Protein (human only)", "E2"),
            ("E3: - Router (heuristic)", "E3"),
            ("E4: - Phase 3 co-ortholog", "E4"),
            ("E5: - C-stage features", "E5"),
        ]:
            m = results[key]
            delta = m["mean"] - e1_f1
            print(f"  {label:<40s} {m['mean']:>4.2f} "
                  f"[{m['ci_lower']:.2f},{m['ci_upper']:.2f}]  {delta:>+8.4f}")

        for tool in TOOL_NAMES:
            m = results[f"E6_{tool}"]
            delta = m["mean"] - e1_f1
            print(f"  {'E6: - ' + tool + ' (A-stage)':<40s} {m['mean']:>4.2f} "
                  f"[{m['ci_lower']:.2f},{m['ci_upper']:.2f}]  {delta:>+8.4f}")

    # ---- Save ----
    summary_path = outdir / "ablation_summary.json"
    with open(summary_path, "w") as f:
        json.dump({
            "species": sp,
            "best_overlap": {k: {kk: vv for kk, vv in v.items() if kk != "metric"}
                             for k, v in all_results["best_overlap"].items()},
            "primary": {k: {kk: vv for kk, vv in v.items() if kk != "metric"}
                        for k, v in all_results["primary"].items()},
            "metadata": {
                "n_rows": len(df),
                "n_queries": df["query_id"].nunique(),
                "random_seed": RANDOM_SEED,
                "n_bootstrap": N_BOOTSTRAP,
            },
        }, f, indent=2, default=float)
    print(f"\n  Summary saved to {summary_path}")

    details_path = outdir / "ablation_details.json"
    with open(details_path, "w") as f:
        json.dump(details_output, f, indent=2, default=float)
    print(f"  Details saved to {details_path}")

    return all_results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description="Phase 5 Multi-Species Ablation")
    parser.add_argument("--species", nargs="+",
                        choices=["mouse", "dog", "chicken", "zebrafish"],
                        help="Species to run ablation on")
    parser.add_argument("--all", action="store_true",
                        help="Run on all 4 species")
    args = parser.parse_args()

    if args.all:
        species_list = ["mouse", "dog", "chicken", "zebrafish"]
    elif args.species:
        species_list = args.species
    else:
        # Default: mouse only (backward compatible)
        species_list = ["mouse"]

    print("=" * 70)
    print(f"Phase 5 Ablation: E1-E6 on {len(species_list)} species: {species_list}")
    print("=" * 70)

    all_species_results = {}
    for sp in species_list:
        results = run_ablation_for_species(sp)
        all_species_results[sp] = results

    # ---- Combined multi-species table ----
    if len(all_species_results) > 1:
        for _, metric_name in LABEL_CONFIGS:
            print(f"\n{'=' * 90}")
            print(f"  COMBINED MULTI-SPECIES ABLATION — {metric_name.upper()}")
            print(f"{'=' * 90}")

            # Header
            header = f"  {'Experiment':<38s}"
            for sp in species_list:
                header += f" {sp:>12s}"
            header += f" {'mean':>10s}"
            print(header)
            print(f"  {'─' * 86}")

            exp_keys = [
                ("E1: Full Phase 5", "E1"),
                ("E2: - Species Protein", "E2"),
                ("E3: - Router", "E3"),
                ("E4: - Phase 3 co-ortholog", "E4"),
                ("E5: - C-stage features", "E5"),
            ] + [(f"E6: - {tool}", f"E6_{tool}") for tool in TOOL_NAMES]

            for label, key in exp_keys:
                line = f"  {label:<38s}"
                means = []
                for sp in species_list:
                    m = all_species_results[sp][metric_name][key]["mean"]
                    line += f" {m:>12.4f}"
                    means.append(m)
                line += f" {np.mean(means):>10.4f}"
                print(line)

            # Deltas
            print(f"\n  {'─' * 86}")
            print(f"  DELTA vs E1")
            print(f"  {'Experiment':<38s}", end="")
            for sp in species_list:
                print(f" {sp:>12s}", end="")
            print(f" {'mean':>10s}")
            print(f"  {'─' * 86}")

            for label, key in exp_keys:
                if key == "E1":
                    continue
                line = f"  {label:<38s}"
                deltas = []
                for sp in species_list:
                    e1 = all_species_results[sp][metric_name]["E1"]["mean"]
                    m = all_species_results[sp][metric_name][key]["mean"]
                    d = m - e1
                    line += f" {d:>+12.4f}"
                    deltas.append(d)
                line += f" {np.mean(deltas):>+10.4f}"
                print(line)

    # ---- Combined JSON ----
    if len(all_species_results) > 1:
        combined = {}
        for sp in species_list:
            combined[sp] = {
                "best_overlap": {
                    k: {kk: vv for kk, vv in v.items() if kk != "metric"}
                    for k, v in all_species_results[sp]["best_overlap"].items()
                },
                "primary": {
                    k: {kk: vv for kk, vv in v.items() if kk != "metric"}
                    for k, v in all_species_results[sp]["primary"].items()
                },
            }
        combined_path = OUTDIR_BASE / "ablation_summary_all_species.json"
        with open(combined_path, "w") as f:
            json.dump(combined, f, indent=2, default=float)
        print(f"\n  Combined summary saved to {combined_path}")

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
