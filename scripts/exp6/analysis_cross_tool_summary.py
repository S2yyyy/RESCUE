"""Compact summary: best cross-tool hybrid vs each single E2E per cell."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

REPO = Path(".").resolve()
HYBRID = REPO / "results/exp6/hybrid"
OUT = REPO / "results/exp6/hybrid/cross_tool_winners.tsv"

SPECIES = ["mouse", "dog", "chicken", "zebrafish"]
REGIMES = ["one2one", "paralog_broad"]
A_TOOLS = ["miniprot", "spaln3", "genblastg", "tblastn", "diamond",
           "mmseqs", "last", "genblasta"]
C_TOOLS = ["miniprot", "spaln3", "gth", "exonerate", "genewise", "scipio", "genblastg"]
E2E_TOOLS = ["miniprot", "spaln3", "gth", "scipio", "genblastg", "minisplice_miniprot"]


def main() -> int:
    rows = []
    for sp in SPECIES:
        for reg in REGIMES:
            npz = np.load(HYBRID / f"bootstrap_{sp}.{reg}.npz")
            H = npz["hybrid_orig"]
            Hb = npz["hybrid_boot"]
            E = npz["e2e_orig"]
            Eb = npz["e2e_boot"]

            cross_pairs = [(H[ai, ci], ai, ci, A_TOOLS[ai], C_TOOLS[ci])
                           for ai in range(8) for ci in range(7)
                           if A_TOOLS[ai] != C_TOOLS[ci]]
            cross_pairs.sort(reverse=True)
            best_cross_f1, best_ai, best_ci, best_a, best_c = cross_pairs[0]
            best_cross_boot = Hb[best_ai, best_ci]

            # Compare best cross-pair against each E2E
            for ei, etool in enumerate(E2E_TOOLS):
                e_orig = float(E[ei])
                e_boot = Eb[ei]
                delta = float(best_cross_f1 - e_orig)
                delta_boot = best_cross_boot - e_boot
                lo = float(np.percentile(delta_boot, 2.5))
                hi = float(np.percentile(delta_boot, 97.5))
                wins = "yes" if (delta > 0 and lo > 0) else "no"
                rows.append(dict(
                    species=sp, regime=reg,
                    best_cross_pair=f"{best_a}+{best_c}",
                    best_cross_f1=f"{best_cross_f1:.4f}",
                    e2e_tool=etool, e2e_f1=f"{e_orig:.4f}",
                    delta=f"{delta:+.4f}",
                    ci_lo=f"{lo:+.4f}", ci_hi=f"{hi:+.4f}",
                    cross_beats_e2e=wins,
                ))

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # Compact print: per cell, count of E2E tools the best cross-pair beats
    print(f"{'cell':<24} {'best cross':<28} {'F1':>6}   beats E2E")
    print("-" * 80)
    for sp in SPECIES:
        for reg in REGIMES:
            cell_rows = [r for r in rows if r["species"] == sp and r["regime"] == reg]
            best_str = f"{cell_rows[0]['best_cross_pair']}"
            f1 = cell_rows[0]['best_cross_f1']
            beaten = [r['e2e_tool'] for r in cell_rows if r['cross_beats_e2e'] == 'yes']
            n_beats = len(beaten)
            beats_str = f"{n_beats}/6: {','.join(beaten)}" if n_beats else "0/6 (none)"
            print(f"{sp+'.'+reg:<24} {best_str:<28} {f1:>6}   {beats_str}")

    print(f"\n  full table → {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
