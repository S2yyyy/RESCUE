"""Exp 6 Block 4 — Per-tool and per-protein stage oracles.

Reads B3's bootstrap NPZ per (sp, regime), tabulates:
  - best_single_e2e_f1          = max_T mean(E2E_T)             (point + CI)
  - oracle_e2e_f1               = mean_p max_T E2E_T(p)         (per-tool oracle)
  - best_hybrid_f1, best_hybrid_pair  = argmax over 8x7 cells   (from B3 argmax)
  - per_protein_oracle_f1       = mean_p (union_search * best_C)
  - union_search_recall         = mean_p max_A A(p)
  - mean_best_C_f1              = mean_p max_C C(p)

All means come with 95% bootstrap CIs from the same paired resamples as B3.

Output:
  results/exp6/oracles/summary.tsv
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

REPO = Path(".").resolve()
HYBRID_DIR = REPO / "results/exp6/hybrid"
ARGMAX_TSV = HYBRID_DIR / "argmax_per_cell.tsv"
OUT_DIR = REPO / "results/exp6/oracles"

SPECIES = ["mouse", "dog", "chicken", "zebrafish"]
REGIMES = ["one2one", "paralog_broad"]
E2E_TOOLS = ["miniprot", "spaln3", "gth", "scipio", "genblastg", "minisplice_miniprot"]


def ci(boot: np.ndarray) -> tuple[float, float]:
    return tuple(float(x) for x in np.percentile(boot, [2.5, 97.5]))


def best_e2e_with_ci(e2e_orig: np.ndarray, e2e_boot: np.ndarray
                     ) -> tuple[str, float, float, float]:
    """Best single e2e tool: pick argmax on original mean, return its CI."""
    ti = int(np.argmax(e2e_orig))
    lo, hi = ci(e2e_boot[ti])
    return E2E_TOOLS[ti], float(e2e_orig[ti]), lo, hi


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load B3 argmax per (sp, regime)
    argmax = {}
    with open(ARGMAX_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            argmax[(r["species"], r["regime"])] = r

    rows = []
    print("=== B4 oracles ===")
    print(f"\n{'sp':<10} {'regime':<14} {'best_e2e':<22} {'oracle_e2e':<14} "
          f"{'best_hybrid':<22} {'pp_oracle':<14} {'union_R':<10} {'best_C':<8}")

    for sp in SPECIES:
        for reg in REGIMES:
            npz = np.load(HYBRID_DIR / f"bootstrap_{sp}.{reg}.npz")
            e2e_orig = npz["e2e_orig"]
            e2e_boot = npz["e2e_boot"]
            oracle_e2e = float(npz["oracle_e2e_orig"])
            oracle_e2e_lo, oracle_e2e_hi = ci(npz["oracle_e2e_boot"])
            oracle_pp = float(npz["oracle_pp_orig"])
            oracle_pp_lo, oracle_pp_hi = ci(npz["oracle_pp_boot"])
            union_r = float(npz["union_recall_orig"])
            union_r_lo, union_r_hi = ci(npz["union_recall_boot"])
            best_C = float(npz["best_C_orig"])
            best_C_lo, best_C_hi = ci(npz["best_C_boot"])

            best_e2e_tool, best_e2e_f1, best_e2e_lo, best_e2e_hi = \
                best_e2e_with_ci(e2e_orig, e2e_boot)

            am = argmax[(sp, reg)]
            best_hybrid_pair = f"{am['A_star']}+{am['C_star']}"
            best_hybrid_f1 = float(am["hybrid_f1"])
            best_hybrid_lo = float(am["ci_lo"])
            best_hybrid_hi = float(am["ci_hi"])

            # Acceptance: oracle_e2e ≥ best single e2e mean
            assert oracle_e2e >= best_e2e_f1 - 1e-6, \
                f"oracle_e2e {oracle_e2e} < best_e2e {best_e2e_f1}"
            # Acceptance: per_protein_oracle ≥ best_hybrid (per-protein routing dominates)
            assert oracle_pp >= best_hybrid_f1 - 1e-6, \
                f"oracle_pp {oracle_pp} < best_hybrid {best_hybrid_f1}"
            # Acceptance: per_protein_oracle ≤ union_recall (upper bound)
            assert oracle_pp <= union_r + 1e-6, \
                f"oracle_pp {oracle_pp} > union_R {union_r}"

            rows.append(dict(
                species=sp, regime=reg,
                best_single_e2e_tool=best_e2e_tool,
                best_single_e2e_f1=f"{best_e2e_f1:.4f}",
                best_single_e2e_ci=f"{best_e2e_lo:.4f}|{best_e2e_hi:.4f}",
                oracle_e2e_f1=f"{oracle_e2e:.4f}",
                oracle_e2e_ci=f"{oracle_e2e_lo:.4f}|{oracle_e2e_hi:.4f}",
                best_hybrid_pair=best_hybrid_pair,
                best_hybrid_f1=f"{best_hybrid_f1:.4f}",
                best_hybrid_ci=f"{best_hybrid_lo:.4f}|{best_hybrid_hi:.4f}",
                per_protein_oracle_f1=f"{oracle_pp:.4f}",
                per_protein_oracle_ci=f"{oracle_pp_lo:.4f}|{oracle_pp_hi:.4f}",
                union_search_recall=f"{union_r:.4f}",
                union_search_recall_ci=f"{union_r_lo:.4f}|{union_r_hi:.4f}",
                mean_best_C_f1=f"{best_C:.4f}",
                mean_best_C_ci=f"{best_C_lo:.4f}|{best_C_hi:.4f}",
            ))

            print(f"  {sp:<10} {reg:<14} "
                  f"{best_e2e_tool[:18]:<18}{best_e2e_f1:>5.3f}  "
                  f"{oracle_e2e:>5.3f}±        "
                  f"{best_hybrid_pair[:20]:<20}{best_hybrid_f1:>5.3f}  "
                  f"{oracle_pp:>5.3f}±        "
                  f"{union_r:>5.3f}     {best_C:>5.3f}")

    out = OUT_DIR / "summary.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"\n  wrote {len(rows)} rows → {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
