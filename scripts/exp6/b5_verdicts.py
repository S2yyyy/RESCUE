"""Exp 6 Block 5 — H1 / H2 / H3 verdicts and recommendation table.

Reads B3 bootstrap NPZ + B4 summary, derives:

  H1: best_hybrid - best_single_e2e ≥ 0.03 with CI excluding 0
      (paired bootstrap on the same 10k resamples as B3)
  H2: per-species (A*, C*) differs across regimes for ≥ 2 species
  H3: per_protein_oracle - best_hybrid ≥ 0.05 with CI excluding 0 on ≥ 1 regime

Outputs:
  results/exp6/verdicts/h1_h2_h3.json
  results/exp6/verdicts/recommendation_table.tsv
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

REPO = Path(".").resolve()
HYBRID_DIR = REPO / "results/exp6/hybrid"
ORACLES_TSV = REPO / "results/exp6/oracles/summary.tsv"
ARGMAX_TSV = HYBRID_DIR / "argmax_per_cell.tsv"
OUT_DIR = REPO / "results/exp6/verdicts"

SPECIES = ["mouse", "dog", "chicken", "zebrafish"]
REGIMES = ["one2one", "paralog_broad"]
A_TOOLS = ["miniprot", "spaln3", "genblastg", "tblastn", "diamond",
           "mmseqs", "last", "genblasta"]
C_TOOLS = ["miniprot", "spaln3", "gth", "exonerate", "genewise", "scipio", "genblastg"]
E2E_TOOLS = ["miniprot", "spaln3", "gth", "scipio", "genblastg", "minisplice_miniprot"]

H1_DELTA_THRESHOLD = 0.03
H3_DELTA_THRESHOLD = 0.05


def percentile_ci(boot: np.ndarray) -> tuple[float, float]:
    return tuple(float(x) for x in np.percentile(boot, [2.5, 97.5]))


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    argmax = {}
    with open(ARGMAX_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            argmax[(r["species"], r["regime"])] = r

    h1_per_cell = {}    # (sp, reg) → {delta, ci, pass}
    h3_per_cell = {}
    rec_rows = []

    for sp in SPECIES:
        for reg in REGIMES:
            npz = np.load(HYBRID_DIR / f"bootstrap_{sp}.{reg}.npz")
            e2e_orig = npz["e2e_orig"]
            e2e_boot = npz["e2e_boot"]
            hybrid_orig = npz["hybrid_orig"]
            hybrid_boot = npz["hybrid_boot"]
            oracle_pp_boot = npz["oracle_pp_boot"]
            oracle_pp_orig = float(npz["oracle_pp_orig"])

            am = argmax[(sp, reg)]
            a_star = A_TOOLS.index(am["A_star"])
            c_star = C_TOOLS.index(am["C_star"])
            best_hybrid_f1 = float(hybrid_orig[a_star, c_star])
            best_hybrid_boot = hybrid_boot[a_star, c_star]   # (B,)

            # Best single e2e — argmax on original, then paired CI on Δ
            t_star = int(np.argmax(e2e_orig))
            best_e2e_tool = E2E_TOOLS[t_star]
            best_e2e_f1 = float(e2e_orig[t_star])
            best_e2e_boot = e2e_boot[t_star]

            # H1
            delta_h1 = best_hybrid_f1 - best_e2e_f1
            delta_h1_boot = best_hybrid_boot - best_e2e_boot
            h1_lo, h1_hi = percentile_ci(delta_h1_boot)
            h1_pass = (delta_h1 >= H1_DELTA_THRESHOLD) and (h1_lo > 0)
            h1_per_cell[f"{sp}.{reg}"] = dict(
                delta=round(float(delta_h1), 4),
                ci_lo=round(float(h1_lo), 4),
                ci_hi=round(float(h1_hi), 4),
                best_hybrid_pair=f"{am['A_star']}+{am['C_star']}",
                best_hybrid_f1=round(best_hybrid_f1, 4),
                best_e2e_tool=best_e2e_tool,
                best_e2e_f1=round(best_e2e_f1, 4),
                threshold=H1_DELTA_THRESHOLD,
                passed=bool(h1_pass),
            )

            # H3
            delta_h3 = oracle_pp_orig - best_hybrid_f1
            delta_h3_boot = oracle_pp_boot - best_hybrid_boot
            h3_lo, h3_hi = percentile_ci(delta_h3_boot)
            h3_pass = (delta_h3 >= H3_DELTA_THRESHOLD) and (h3_lo > 0)
            h3_per_cell[f"{sp}.{reg}"] = dict(
                delta=round(float(delta_h3), 4),
                ci_lo=round(float(h3_lo), 4),
                ci_hi=round(float(h3_hi), 4),
                per_protein_oracle_f1=round(oracle_pp_orig, 4),
                best_hybrid_f1=round(best_hybrid_f1, 4),
                threshold=H3_DELTA_THRESHOLD,
                passed=bool(h3_pass),
            )

            rec_rows.append(dict(
                regime=reg, species=sp,
                best_recipe=f"{am['A_star']}+{am['C_star']}",
                best_recipe_f1=f"{best_hybrid_f1:.4f}",
                best_single_e2e_tool=best_e2e_tool,
                best_single_e2e_f1=f"{best_e2e_f1:.4f}",
                delta_h1=f"{delta_h1:+.4f}",
                delta_h1_ci=f"{h1_lo:+.4f}|{h1_hi:+.4f}",
                h1_verdict="PASS" if h1_pass else "FAIL",
                per_protein_oracle_f1=f"{oracle_pp_orig:.4f}",
                gap_to_ceiling=f"{delta_h3:+.4f}",
                gap_to_ceiling_ci=f"{h3_lo:+.4f}|{h3_hi:+.4f}",
                h3_verdict="PASS" if h3_pass else "FAIL",
            ))

    # H2: per-species pairs differ across regimes
    h2_per_species = {}
    n_diff_species = 0
    for sp in SPECIES:
        a1 = (argmax[(sp, "one2one")]["A_star"],
              argmax[(sp, "one2one")]["C_star"])
        a2 = (argmax[(sp, "paralog_broad")]["A_star"],
              argmax[(sp, "paralog_broad")]["C_star"])
        differs = (a1 != a2)
        h2_per_species[sp] = dict(
            one2one_pair=f"{a1[0]}+{a1[1]}",
            paralog_broad_pair=f"{a2[0]}+{a2[1]}",
            differs=bool(differs),
        )
        if differs:
            n_diff_species += 1
    h2_pass = (n_diff_species >= 2)

    # H1 overall: PASS on N/8 cells; meaningful PASS = ≥ 1 regime fully PASS
    n_h1_pass = sum(1 for v in h1_per_cell.values() if v["passed"])
    by_regime_pass = {
        reg: sum(1 for sp in SPECIES if h1_per_cell[f"{sp}.{reg}"]["passed"])
        for reg in REGIMES
    }
    # spec wording: "PASS on N/8 cells, ≥ 1 regime fully PASS"
    # interpret "fully PASS" as all 4 species PASS in that regime
    h1_overall_pass = any(by_regime_pass[r] == len(SPECIES) for r in REGIMES)
    n_h3_pass = sum(1 for v in h3_per_cell.values() if v["passed"])
    h3_overall_pass = (n_h3_pass >= 1)

    verdicts = dict(
        H1=dict(
            description=(
                "best stage-level hybrid (A*,C*) beats best single e2e tool "
                f"by Δ ≥ {H1_DELTA_THRESHOLD} (CI excludes 0)"),
            per_cell=h1_per_cell,
            n_pass=n_h1_pass,
            n_total=len(h1_per_cell),
            by_regime_pass=by_regime_pass,
            overall_pass=bool(h1_overall_pass),
        ),
        H2=dict(
            description="per-species (A*, C*) differs between regimes for ≥ 2 species",
            per_species=h2_per_species,
            n_diff_species=n_diff_species,
            overall_pass=bool(h2_pass),
        ),
        H3=dict(
            description=(
                "per-protein stage oracle exceeds best fixed hybrid "
                f"by Δ ≥ {H3_DELTA_THRESHOLD} on ≥ 1 (sp, regime)"),
            per_cell=h3_per_cell,
            n_pass=n_h3_pass,
            n_total=len(h3_per_cell),
            overall_pass=bool(h3_overall_pass),
        ),
    )

    out_json = OUT_DIR / "h1_h2_h3.json"
    out_json.write_text(json.dumps(verdicts, indent=2, ensure_ascii=False))

    # Recommendation table
    rec_path = OUT_DIR / "recommendation_table.tsv"
    with open(rec_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rec_rows[0].keys()),
                           delimiter="\t")
        w.writeheader()
        # group by regime for readability
        for reg in REGIMES:
            for r in rec_rows:
                if r["regime"] == reg:
                    w.writerow(r)

    # Print summary
    print("=== B5 verdicts ===\n")
    print("H1 (hybrid > best single e2e, Δ ≥ 0.03):")
    print(f"  PASS on {n_h1_pass}/{len(h1_per_cell)} cells; "
          f"by regime: {by_regime_pass}; overall PASS = {h1_overall_pass}")
    for cell, v in h1_per_cell.items():
        marker = "✓" if v["passed"] else "✗"
        print(f"    {marker} {cell:<22} Δ={v['delta']:+.4f} "
              f"CI=[{v['ci_lo']:+.4f}, {v['ci_hi']:+.4f}]  "
              f"({v['best_hybrid_pair']} vs {v['best_e2e_tool']})")

    print("\nH2 (regime-specific recipe — pairs differ across regimes):")
    print(f"  {n_diff_species}/4 species differ; overall PASS = {h2_pass}")
    for sp, v in h2_per_species.items():
        marker = "✓" if v["differs"] else "✗"
        print(f"    {marker} {sp:<10} 1:1 → {v['one2one_pair']:<25} "
              f"paralog → {v['paralog_broad_pair']}")

    print("\nH3 (per-protein oracle > best hybrid, Δ ≥ 0.05):")
    print(f"  PASS on {n_h3_pass}/{len(h3_per_cell)} cells; "
          f"overall PASS = {h3_overall_pass}")
    for cell, v in h3_per_cell.items():
        marker = "✓" if v["passed"] else "✗"
        print(f"    {marker} {cell:<22} Δ={v['delta']:+.4f} "
              f"CI=[{v['ci_lo']:+.4f}, {v['ci_hi']:+.4f}]  "
              f"(oracle {v['per_protein_oracle_f1']} vs hybrid {v['best_hybrid_f1']})")

    print(f"\n  verdicts → {out_json}")
    print(f"  recommendation table → {rec_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
