"""Exp 6 Block 6 — Anti-claims AC1 (identity confound) and AC2 (easy-find).

For each paralog cell (the regime where H1 PASSed in B5), check whether
the H1 hybrid headroom persists on filtered subsets:

  AC1: restrict to identity ≥ 90.
       Pass if H1 Δ remains ≥ 0.03 with CI excluding 0.
       Fail = headroom is identity-confounded (B5's H1 inherits the
       Exp 5 C4 identity confound).

  AC2: restrict to proteins where ≥ 4 of 8 A-tools have search_ok=1.
       Pass if H1 Δ persists.
       Fail = headroom is solely from A-stage coverage; pick a better
       search tool and the recipe doesn't add value.

Re-bootstraps on each subset (seed = B3 cell seed + 50000).

Output:
  results/exp6/verdicts/anticlaims.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(".").resolve()
JOINED = REPO / "results/exp6/joined/per_protein.tsv"
HYBRID_DIR = REPO / "results/exp6/hybrid"
OUT_DIR = REPO / "results/exp6/verdicts"

SPECIES = ["mouse", "dog", "chicken", "zebrafish"]
REGIME_TARGET = "paralog_broad"   # AC verdicts apply to H1-passing regime
A_TOOLS = ["miniprot", "spaln3", "genblastg", "tblastn", "diamond",
           "mmseqs", "last", "genblasta"]
C_TOOLS = ["miniprot", "spaln3", "gth", "exonerate", "genewise", "scipio", "genblastg"]
E2E_TOOLS = ["miniprot", "spaln3", "gth", "scipio", "genblastg", "minisplice_miniprot"]

N_BOOT = 10000
BATCH = 200
SEED_BASE = 4242
SEED_OFFSET = 50000
H1_DELTA_THRESHOLD = 0.03

ID_BIN_LO = 90.0
AC2_MIN_N_A = 4


def cell_seed(sp: str, reg: str) -> int:
    return SEED_BASE + SPECIES.index(sp) * 10 + [
        "one2one", "paralog_broad"].index(reg) + SEED_OFFSET


def load_subset(df: pd.DataFrame, sp: str, reg: str, mask_fn) -> dict:
    sub = df[(df["species"] == sp) & (df["regime"] == reg)].reset_index(drop=True)
    A = sub[[f"A_{t}" for t in A_TOOLS]].to_numpy(dtype=np.float32)
    C = np.nan_to_num(sub[[f"C_{t}" for t in C_TOOLS]].to_numpy(dtype=np.float32), nan=0.0)
    E = np.nan_to_num(sub[[f"E2E_{t}" for t in E2E_TOOLS]].to_numpy(dtype=np.float32), nan=0.0)
    identity = sub["identity"].to_numpy(dtype=np.float32)
    n_a_found = A.sum(axis=1)
    keep = mask_fn(identity=identity, n_a_found=n_a_found)
    return dict(
        n_full=len(sub), n=int(keep.sum()),
        A=A[keep], C=C[keep], E=E[keep],
    )


def hybrid_h1_delta(cell: dict, seed: int) -> dict:
    """Re-bootstrap hybrid matrix on subset; return new argmax + paired Δ vs best e2e."""
    A, C, E = cell["A"], cell["C"], cell["E"]
    n = cell["n"]
    if n < 5:
        return dict(n=n, skipped=True, reason=f"n={n} too small")

    H = (A[:, :, None] * C[:, None, :]).reshape(n, -1)   # (n, 56)
    M = np.column_stack([H, E]).astype(np.float32)       # (n, 56+6)

    rng = np.random.default_rng(seed)
    boot_means = np.empty((N_BOOT, M.shape[1]), dtype=np.float32)
    for b0 in range(0, N_BOOT, BATCH):
        b1 = min(b0 + BATCH, N_BOOT)
        idx = rng.integers(0, n, size=(b1 - b0, n))
        boot_means[b0:b1] = M[idx].mean(axis=1)

    orig = M.mean(axis=0)
    hybrid_orig = orig[:56].reshape(8, 7)
    e2e_orig = orig[56:62]

    a_star, c_star = np.unravel_index(hybrid_orig.argmax(), hybrid_orig.shape)
    t_star = int(e2e_orig.argmax())

    best_hybrid = float(hybrid_orig[a_star, c_star])
    best_e2e = float(e2e_orig[t_star])
    delta = best_hybrid - best_e2e

    hybrid_boot = boot_means[:, :56].reshape(N_BOOT, 8, 7)[:, a_star, c_star]
    e2e_boot = boot_means[:, 56 + t_star]
    delta_boot = hybrid_boot - e2e_boot
    lo, hi = (float(x) for x in np.percentile(delta_boot, [2.5, 97.5]))

    passed = (delta >= H1_DELTA_THRESHOLD) and (lo > 0)
    return dict(
        n=int(n), skipped=False,
        best_hybrid_pair=f"{A_TOOLS[a_star]}+{C_TOOLS[c_star]}",
        best_hybrid_f1=round(best_hybrid, 4),
        best_e2e_tool=E2E_TOOLS[t_star],
        best_e2e_f1=round(best_e2e, 4),
        delta=round(delta, 4),
        ci_lo=round(lo, 4), ci_hi=round(hi, 4),
        threshold=H1_DELTA_THRESHOLD,
        passed=bool(passed),
    )


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(JOINED, sep="\t")

    print("=== B6 anti-claims ===\n")

    out = dict(
        AC1_identity_ge_90=dict(
            description=("Restrict paralog cells to identity ≥ 90; "
                         "check whether H1 Δ ≥ 0.03 with CI excluding 0 persists."),
            threshold_identity=ID_BIN_LO,
            per_cell={},
        ),
        AC2_n_a_ge_4=dict(
            description=("Restrict paralog cells to proteins with ≥ 4 / 8 "
                         "A-tools finding GT; check H1 Δ persistence."),
            threshold_n_a=AC2_MIN_N_A,
            per_cell={},
        ),
    )

    for sp in SPECIES:
        cell_id = f"{sp}.{REGIME_TARGET}"

        # AC1 — identity ≥ 90
        ac1 = load_subset(df, sp, REGIME_TARGET,
                          mask_fn=lambda identity, n_a_found:
                              identity >= ID_BIN_LO)
        ac1_res = hybrid_h1_delta(ac1, cell_seed(sp, REGIME_TARGET) + 1)
        ac1_res["n_full"] = ac1["n_full"]
        out["AC1_identity_ge_90"]["per_cell"][cell_id] = ac1_res

        # AC2 — at least 4 of 8 A-tools find GT
        ac2 = load_subset(df, sp, REGIME_TARGET,
                          mask_fn=lambda identity, n_a_found:
                              n_a_found >= AC2_MIN_N_A)
        ac2_res = hybrid_h1_delta(ac2, cell_seed(sp, REGIME_TARGET) + 2)
        ac2_res["n_full"] = ac2["n_full"]
        out["AC2_n_a_ge_4"]["per_cell"][cell_id] = ac2_res

    # Print summary
    for ac_key in ("AC1_identity_ge_90", "AC2_n_a_ge_4"):
        block = out[ac_key]
        n_pass = sum(1 for v in block["per_cell"].values()
                     if not v.get("skipped") and v.get("passed"))
        print(f"{ac_key}: PASS on {n_pass}/{len(block['per_cell'])} cells")
        for cell, v in block["per_cell"].items():
            if v.get("skipped"):
                print(f"    – {cell:<22} skipped: {v.get('reason')}")
                continue
            marker = "✓" if v["passed"] else "✗"
            print(f"    {marker} {cell:<22} n={v['n']:>4}/{v['n_full']:<4} "
                  f"Δ={v['delta']:+.4f} CI=[{v['ci_lo']:+.4f}, {v['ci_hi']:+.4f}]  "
                  f"({v['best_hybrid_pair']} vs {v['best_e2e_tool']})")
        block["overall_pass"] = bool(n_pass >= 1)
        print()

    out_path = OUT_DIR / "anticlaims.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(f"  → {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
