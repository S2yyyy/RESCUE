"""Exp 6 Block 3 — Hybrid F1 matrix per (species, regime).

For each (species, regime) computes the 8x7 hybrid F1 matrix
  H[a][c] = E_p [ A_<a>(p) * (C_<c>(p) or 0) ]
plus 10k paired bootstrap CIs per cell, plus argmax (A*, C*).

Per Q3.1 (default), C-tool kernel failures (NaN) are treated as F1=0.
Per Q3.2 (default), all 56 (A, C) pairs are reported.

Bootstrap means tensors for hybrid + e2e + per-tool A/C + oracles are
persisted as NPZ so B4/B5 reuse the same paired resamples.

Outputs:
  results/exp6/hybrid/matrix_<sp>.<regime>.tsv         mean F1, 8x7
  results/exp6/hybrid/matrix_ci_<sp>.<regime>.tsv      "lo|hi" strings, 8x7
  results/exp6/hybrid/bootstrap_<sp>.<regime>.npz      boot means tensor
  results/exp6/hybrid/argmax_per_cell.tsv              one row per (sp, reg)
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(".").resolve()
JOINED = REPO / "results/exp6/joined/per_protein.tsv"
OUT_DIR = REPO / "results/exp6/hybrid"

SPECIES = ["mouse", "dog", "chicken", "zebrafish"]
REGIMES = ["one2one", "paralog_broad"]
A_TOOLS = ["miniprot", "spaln3", "genblastg", "tblastn", "diamond",
           "mmseqs", "last", "genblasta"]
C_TOOLS = ["miniprot", "spaln3", "gth", "exonerate", "genewise", "scipio", "genblastg"]
E2E_TOOLS = ["miniprot", "spaln3", "gth", "scipio", "genblastg", "minisplice_miniprot"]

N_BOOT = 10000
SEED_BASE = 4242
BATCH = 200  # resamples per chunk to keep memory < 200 MB


def cell_seed(sp: str, reg: str) -> int:
    return SEED_BASE + SPECIES.index(sp) * 10 + REGIMES.index(reg)


def load_cell(df: pd.DataFrame, sp: str, reg: str) -> dict:
    sub = df[(df["species"] == sp) & (df["regime"] == reg)].reset_index(drop=True)
    n = len(sub)
    A = sub[[f"A_{t}" for t in A_TOOLS]].to_numpy(dtype=np.float32)
    C = np.nan_to_num(
        sub[[f"C_{t}" for t in C_TOOLS]].to_numpy(dtype=np.float32), nan=0.0
    )
    E = np.nan_to_num(
        sub[[f"E2E_{t}" for t in E2E_TOOLS]].to_numpy(dtype=np.float32), nan=0.0
    )
    identity = sub["identity"].to_numpy(dtype=np.float32)
    proteins = sub["protein_id"].to_numpy()
    return dict(n=n, A=A, C=C, E=E, identity=identity, proteins=proteins)


def bootstrap_cell(cell: dict, seed: int) -> dict:
    """Compute original means + 10k paired bootstrap means for all metrics."""
    A, C, E = cell["A"], cell["C"], cell["E"]
    n = cell["n"]

    # Per-protein metrics
    H = (A[:, :, None] * C[:, None, :]).reshape(n, -1)        # (n, 56)
    union_search = A.max(axis=1)                              # (n,)
    best_C = C.max(axis=1)                                    # (n,)
    oracle_e2e_p = E.max(axis=1)                              # (n,)
    oracle_pp_p = union_search * best_C                       # (n,)
    n_a_found = A.sum(axis=1)                                 # (n,) for AC2 filter

    M = np.column_stack([
        H,                              # 0:56  hybrid
        E,                              # 56:62 e2e
        A,                              # 62:70 a-stage
        C,                              # 70:77 c-stage
        oracle_e2e_p[:, None],          # 77   per-tool oracle (e2e max)
        oracle_pp_p[:, None],           # 78   per-protein stage oracle
        union_search[:, None],          # 79   union recall
        best_C[:, None],                # 80   mean best-C
    ]).astype(np.float32)
    M_DIM = M.shape[1]

    rng = np.random.default_rng(seed)
    boot_means = np.empty((N_BOOT, M_DIM), dtype=np.float32)
    for b0 in range(0, N_BOOT, BATCH):
        b1 = min(b0 + BATCH, N_BOOT)
        idx_chunk = rng.integers(0, n, size=(b1 - b0, n))
        boot_means[b0:b1] = M[idx_chunk].mean(axis=1)

    orig_means = M.mean(axis=0)

    return dict(
        n=n,
        n_a_found=n_a_found,
        identity=cell["identity"],
        proteins=cell["proteins"],
        # original (point estimate)
        hybrid_orig=orig_means[:56].reshape(8, 7),
        e2e_orig=orig_means[56:62],
        a_orig=orig_means[62:70],
        c_orig=orig_means[70:77],
        oracle_e2e_orig=float(orig_means[77]),
        oracle_pp_orig=float(orig_means[78]),
        union_recall_orig=float(orig_means[79]),
        best_C_orig=float(orig_means[80]),
        # bootstrap means (axis=-1 is resample)
        hybrid_boot=boot_means[:, :56].reshape(N_BOOT, 8, 7).transpose(1, 2, 0),
        e2e_boot=boot_means[:, 56:62].T,
        a_boot=boot_means[:, 62:70].T,
        c_boot=boot_means[:, 70:77].T,
        oracle_e2e_boot=boot_means[:, 77],
        oracle_pp_boot=boot_means[:, 78],
        union_recall_boot=boot_means[:, 79],
        best_C_boot=boot_means[:, 80],
    )


def write_matrix_tsv(path: Path, M: np.ndarray) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["A\\C"] + C_TOOLS)
        for ai, atool in enumerate(A_TOOLS):
            w.writerow([atool] + [f"{M[ai, ci]:.4f}" for ci in range(7)])


def write_ci_tsv(path: Path, lo: np.ndarray, hi: np.ndarray) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["A\\C"] + C_TOOLS)
        for ai, atool in enumerate(A_TOOLS):
            w.writerow([atool] +
                       [f"{lo[ai, ci]:.4f}|{hi[ai, ci]:.4f}" for ci in range(7)])


def acceptance_diagonal(boot: dict, sp: str, reg: str) -> list[str]:
    """Self-pair (A=C) cells should match natural E2E within 0.03."""
    msgs = []
    self_pairs = [t for t in A_TOOLS if t in C_TOOLS and t in E2E_TOOLS]
    for t in self_pairs:
        ai, ci = A_TOOLS.index(t), C_TOOLS.index(t)
        ei = E2E_TOOLS.index(t)
        diag = boot["hybrid_orig"][ai, ci]
        e2e = boot["e2e_orig"][ei]
        diff = diag - e2e
        flag = "" if abs(diff) <= 0.03 else "  WARN"
        msgs.append(f"  {sp:<10} {reg:<14} {t:<10} hybrid={diag:.4f} "
                    f"e2e={e2e:.4f} Δ={diff:+.4f}{flag}")
    return msgs


def acceptance_argmax(boot: dict, sp: str, reg: str,
                      a_star: int, c_star: int, lo: float) -> str:
    """argmax cell ≥ best single e2e tool's mean; CI lo > 0.5."""
    argmax_f1 = boot["hybrid_orig"][a_star, c_star]
    best_e2e = boot["e2e_orig"].max()
    best_e2e_tool = E2E_TOOLS[boot["e2e_orig"].argmax()]
    geq = "" if argmax_f1 >= best_e2e else "  WARN: argmax < best e2e"
    sanity = "" if lo > 0.5 else "  WARN: CI lo <= 0.5"
    return (f"  {sp:<10} {reg:<14} argmax=({A_TOOLS[a_star]},{C_TOOLS[c_star]})"
            f" F1={argmax_f1:.4f} best_e2e=({best_e2e_tool}){best_e2e:.4f}"
            f"{geq}{sanity}")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(JOINED, sep="\t")

    argmax_rows = []
    print("=== B3 hybrid matrix ===")
    print("\nDiagonal acceptance (|hybrid_self - e2e| ≤ 0.03):")

    diag_msgs: list[str] = []
    argmax_msgs: list[str] = []

    for sp in SPECIES:
        for reg in REGIMES:
            cell = load_cell(df, sp, reg)
            seed = cell_seed(sp, reg)
            boot = bootstrap_cell(cell, seed)

            mean_M = boot["hybrid_orig"]
            ci = np.percentile(boot["hybrid_boot"], [2.5, 97.5], axis=-1)
            lo_M, hi_M = ci[0], ci[1]

            write_matrix_tsv(OUT_DIR / f"matrix_{sp}.{reg}.tsv", mean_M)
            write_ci_tsv(OUT_DIR / f"matrix_ci_{sp}.{reg}.tsv", lo_M, hi_M)

            # argmax
            flat_idx = mean_M.argmax()
            a_star, c_star = np.unravel_index(flat_idx, mean_M.shape)
            argmax_f1 = float(mean_M[a_star, c_star])
            argmax_lo = float(lo_M[a_star, c_star])
            argmax_hi = float(hi_M[a_star, c_star])
            # gap to second-best
            flat_sorted = np.sort(mean_M.ravel())[::-1]
            gap = float(flat_sorted[0] - flat_sorted[1])

            argmax_rows.append(dict(
                species=sp, regime=reg,
                A_star=A_TOOLS[a_star], C_star=C_TOOLS[c_star],
                hybrid_f1=f"{argmax_f1:.4f}",
                ci_lo=f"{argmax_lo:.4f}", ci_hi=f"{argmax_hi:.4f}",
                gap_to_second_best=f"{gap:.4f}",
            ))

            # NPZ persistence for B4/B5/B6
            np.savez_compressed(
                OUT_DIR / f"bootstrap_{sp}.{reg}.npz",
                proteins=boot["proteins"],
                identity=boot["identity"],
                n_a_found=boot["n_a_found"],
                hybrid_orig=boot["hybrid_orig"],
                e2e_orig=boot["e2e_orig"],
                a_orig=boot["a_orig"],
                c_orig=boot["c_orig"],
                oracle_e2e_orig=boot["oracle_e2e_orig"],
                oracle_pp_orig=boot["oracle_pp_orig"],
                union_recall_orig=boot["union_recall_orig"],
                best_C_orig=boot["best_C_orig"],
                hybrid_boot=boot["hybrid_boot"],
                e2e_boot=boot["e2e_boot"],
                a_boot=boot["a_boot"],
                c_boot=boot["c_boot"],
                oracle_e2e_boot=boot["oracle_e2e_boot"],
                oracle_pp_boot=boot["oracle_pp_boot"],
                union_recall_boot=boot["union_recall_boot"],
                best_C_boot=boot["best_C_boot"],
            )

            diag_msgs.extend(acceptance_diagonal(boot, sp, reg))
            argmax_msgs.append(acceptance_argmax(
                boot, sp, reg, a_star, c_star, argmax_lo))

    print("\n".join(diag_msgs))
    print("\nArgmax acceptance (≥ best e2e mean; CI lo > 0.5):")
    print("\n".join(argmax_msgs))

    # argmax summary
    argmax_path = OUT_DIR / "argmax_per_cell.tsv"
    with open(argmax_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(argmax_rows[0].keys()),
                           delimiter="\t")
        w.writeheader()
        w.writerows(argmax_rows)

    print(f"\n  argmax → {argmax_path}")
    print(f"  matrices + bootstraps → {OUT_DIR}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
