"""Cross-tool hybrid vs every single E2E — per-cell winner table.

Question: are there cross-tool (A != C) hybrid pipelines that beat
each existing single E2E tool? We rank all 56 (A, C) cells against
the 6 single-tool E2E baselines per (sp, regime).

Per cell report:
  - best non-minisplice single E2E (the realistic "stage-decomposable" baseline)
  - best self-pair hybrid (A == C, only for tools in both panels)
  - best cross-tool hybrid (A != C)
  - which cross-tool pairs strictly beat each E2E baseline by Δ >= 0.01
    with paired-bootstrap CI excluding 0
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

REPO = Path(".").resolve()
HYBRID = REPO / "results/exp6/hybrid"

SPECIES = ["mouse", "dog", "chicken", "zebrafish"]
REGIMES = ["one2one", "paralog_broad"]
A_TOOLS = ["miniprot", "spaln3", "genblastg", "tblastn", "diamond",
           "mmseqs", "last", "genblasta"]
C_TOOLS = ["miniprot", "spaln3", "gth", "exonerate", "genewise", "scipio", "genblastg"]
E2E_TOOLS = ["miniprot", "spaln3", "gth", "scipio", "genblastg", "minisplice_miniprot"]

DELTA_THRESHOLD = 0.01    # softer than H1's 0.03


def main() -> int:
    print("=" * 80)
    print("Cross-tool hybrid vs single E2E — per-cell winners")
    print("(A != C; Δ ≥ 0.01 with paired-bootstrap 95% CI excluding 0)")
    print("=" * 80)

    for sp in SPECIES:
        for reg in REGIMES:
            npz = np.load(HYBRID / f"bootstrap_{sp}.{reg}.npz")
            H = npz["hybrid_orig"]               # (8, 7)
            Hb = npz["hybrid_boot"]              # (8, 7, 10000)
            E = npz["e2e_orig"]                  # (6,)
            Eb = npz["e2e_boot"]                 # (6, 10000)

            # Best non-minisplice single E2E
            non_ms_idx = [i for i, t in enumerate(E2E_TOOLS) if t != "minisplice_miniprot"]
            best_non_ms = max(non_ms_idx, key=lambda i: E[i])
            best_ms = E.argmax()

            # Self-pair vs cross-pair
            self_pairs = []
            cross_pairs = []
            for ai, atool in enumerate(A_TOOLS):
                for ci, ctool in enumerate(C_TOOLS):
                    item = (H[ai, ci], ai, ci, atool, ctool)
                    if atool == ctool:
                        self_pairs.append(item)
                    else:
                        cross_pairs.append(item)

            self_pairs.sort(reverse=True)
            cross_pairs.sort(reverse=True)
            best_self = self_pairs[0] if self_pairs else None
            best_cross = cross_pairs[0]

            print(f"\n--- {sp} {reg} ---")
            print(f"  best single E2E:")
            print(f"    incl. minisplice  : {E2E_TOOLS[best_ms]:<22} {E[best_ms]:.4f}")
            print(f"    excl. minisplice  : {E2E_TOOLS[best_non_ms]:<22} {E[best_non_ms]:.4f}")
            if best_self:
                print(f"  best self-pair       : {best_self[3]}+{best_self[4]:<18} "
                      f"{best_self[0]:.4f}")
            print(f"  best cross-pair      : {best_cross[3]}+{best_cross[4]:<18} "
                  f"{best_cross[0]:.4f}")

            # For each E2E baseline, which cross-tool pairs strictly beat it?
            print(f"  cross-tool pairs beating each E2E (Δ ≥ {DELTA_THRESHOLD}, CI > 0):")
            for ei, etool in enumerate(E2E_TOOLS):
                e_orig = E[ei]
                e_boot = Eb[ei]
                winners = []
                for f1, ai, ci, atool, ctool in cross_pairs:
                    delta = f1 - e_orig
                    if delta < DELTA_THRESHOLD:
                        continue
                    delta_boot = Hb[ai, ci] - e_boot
                    lo = float(np.percentile(delta_boot, 2.5))
                    if lo > 0:
                        winners.append((delta, atool, ctool, lo,
                                        float(np.percentile(delta_boot, 97.5))))
                if not winners:
                    print(f"    vs {etool:<22} {e_orig:.4f}  (none)")
                else:
                    winners.sort(reverse=True)
                    head = winners[:3]
                    extra = f" (+{len(winners) - 3} more)" if len(winners) > 3 else ""
                    print(f"    vs {etool:<22} {e_orig:.4f}  → "
                          f"{len(winners)} winners{extra}; top:")
                    for d, a, c, lo, hi in head:
                        print(f"         {a}+{c:<14} Δ=+{d:.4f} CI=[{lo:+.4f}, {hi:+.4f}]")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
