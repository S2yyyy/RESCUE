"""Compute multiple heuristic baselines for E3 (-Router) ablation.

Each heuristic replaces the XGBoost Router: for each query, select one
candidate region using a fixed rule instead of the learned model.

Heuristics:
  1. always-primary  — on_primary_chrom==1, tiebreak coorth_overlap_bp (current E3)
  2. max-C-stage      — max c_composite_score
  3. max-identity     — max c_identity
  4. max-tool-consensus — max a_n_tools
  5. max-C-coverage   — max c_coverage

Uses the same 80/20 GroupShuffleSplit as the main E1-E6 experiments.
Reports both primary F1 and best_overlap F1.
"""

import json, sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import PROJECT_ROOT  # noqa: E402

DATA_DIR = PROJECT_ROOT / "results/tprs/phase5"
OUTPUT = DATA_DIR / "heuristic_baselines.json"
SEED = 4242
N_BOOTSTRAP = 10_000

SPECIES = ["mouse", "dog", "chicken", "zebrafish"]


def bootstrap_mean_ci(values, n_bootstrap=N_BOOTSTRAP, seed=SEED):
    rng = np.random.default_rng(seed)
    means = []
    for _ in range(n_bootstrap):
        sample = rng.choice(values, size=len(values), replace=True)
        means.append(np.mean(sample))
    means = np.sort(means)
    lo = means[int(n_bootstrap * 0.025)]
    hi = means[int(n_bootstrap * 0.975)]
    return float(np.mean(values)), float(lo), float(hi)


def heuristic_select(qdf, rule):
    """Apply one heuristic rule to select a row index from qdf."""
    sub = qdf[qdf["query_type"] == "human"].copy()
    if len(sub) == 0:
        return 0  # fallback: first row

    if rule == "always_primary":
        prim = sub[sub["on_primary_chrom"] == 1]
        if len(prim) > 0:
            return int(prim.sort_values("coorth_overlap_bp", ascending=False).index[0])
        return int(sub.index[0])

    elif rule == "max_cstage":
        return int(sub.sort_values("c_composite_score", ascending=False).index[0])

    elif rule == "max_identity":
        return int(sub.sort_values("c_identity", ascending=False).index[0])

    elif rule == "max_tool_consensus":
        return int(sub.sort_values("a_n_tools", ascending=False).index[0])

    elif rule == "max_coverage":
        return int(sub.sort_values("c_coverage", ascending=False).index[0])

    else:
        raise ValueError(f"Unknown rule: {rule}")


def _find_training_data(sp):
    """Try multiple paths for training data."""
    candidates = [
        DATA_DIR / f"{sp}/training_data.tsv",
        DATA_DIR / f"training_data_{sp}.tsv",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(f"No training data for {sp} at {candidates}")


def evaluate_species(sp, rules):
    path = _find_training_data(sp)
    print(f"  Loading {path}")
    df = pd.read_csv(path, sep="\t")
    queries = df["query_id"].unique()

    gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=SEED)
    train_idx, test_idx = next(gss.split(queries, groups=queries))
    test_queries = set(queries[test_idx])

    results = {}
    for rule_name in rules:
        primary_f1s = []
        best_overlap_f1s = []
        for qid in test_queries:
            qdf = df[df["query_id"] == qid]
            best_idx = heuristic_select(qdf, rule_name)
            row = qdf.loc[best_idx]
            primary_f1s.append(row["label_primary_f1"])
            best_overlap_f1s.append(row["label_best_overlap_f1"])

        p_mean, p_lo, p_hi = bootstrap_mean_ci(np.array(primary_f1s))
        b_mean, b_lo, b_hi = bootstrap_mean_ci(np.array(best_overlap_f1s))
        results[rule_name] = {
            "primary_f1": {"mean": round(p_mean, 4), "ci_lo": round(p_lo, 4), "ci_hi": round(p_hi, 4)},
            "best_overlap_f1": {"mean": round(b_mean, 4), "ci_lo": round(b_lo, 4), "ci_hi": round(b_hi, 4)},
            "n_queries": len(test_queries),
        }

    return results


def main():
    rules = ["always_primary", "max_cstage", "max_identity", "max_tool_consensus", "max_coverage"]

    all_results = {}
    for sp in SPECIES:
        print(f"Evaluating {sp}...")
        all_results[sp] = evaluate_species(sp, rules)
        for rule_name in rules:
            r = all_results[sp][rule_name]
            print(f"  {rule_name:25s}  primary={r['primary_f1']['mean']:.4f}  "
                  f"best_overlap={r['best_overlap_f1']['mean']:.4f}")

    # Compute cross-species mean for each heuristic (both metrics)
    summary = {}
    for rule_name in rules:
        summary[rule_name] = {}
        for metric in ["primary_f1", "best_overlap_f1"]:
            means = [all_results[sp][rule_name][metric]["mean"] for sp in SPECIES]
            summary[rule_name][f"mean_{metric}"] = round(np.mean(means), 4)

    all_results["_summary"] = summary

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUTPUT, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved to {OUTPUT}")
    print(f"\nCross-species means:")
    for rule_name in rules:
        s = summary[rule_name]
        print(f"  {rule_name:25s}  primary={s['mean_primary_f1']:.4f}  "
              f"best_overlap={s['mean_best_overlap_f1']:.4f}")


if __name__ == "__main__":
    main()
