"""Exp 6 Block 2 — Build unified per-protein join table.

Joins four sources on (species, regime, protein_id):
  - Exp 1 baseline natural-pipeline exon F1 (7 e2e tools)
  - Exp 5 forced-locus exon F1 at pad=50k (7 C-tools)
  - B1 A-stage search_ok per-protein (8 A-tools)
  - protein identity from ground_truth.json

Regime mapping: Exp 1 uses `one_to_one`/`paralog`, canonicalized to
`one2one`/`paralog_broad` (Q2.1 default per EXPERIMENT_PLAN_EXP6.md).

Output:
  results/exp6/joined/per_protein.tsv   ~10k rows × ~24 cols
  results/exp6/joined/schema.md         column descriptions
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

REPO = Path(".").resolve()
EXP1_TSV = REPO / "results/exp1/per_protein.tsv"
EXP5_DIR = REPO / "results/exp5/per_protein"
B1_DIR = REPO / "results/exp6/a_search"
GT_DIR = REPO / "results/exp2/gt"
DATA_V3 = REPO / "data/v3"
OUT_DIR = REPO / "results/exp6/joined"

SPECIES = ["mouse", "dog", "chicken", "zebrafish"]
REGIMES = ["one2one", "paralog_broad"]

A_TOOLS = ["miniprot", "spaln3", "genblastg", "tblastn", "diamond",
           "mmseqs", "last", "genblasta"]
C_TOOLS = ["miniprot", "spaln3", "gth", "exonerate", "genewise", "scipio", "genblastg"]
E2E_TOOLS = ["miniprot", "spaln3", "gth", "scipio", "genblastg", "minisplice_miniprot"]

EXP1_REGIME_MAP = {"one2one": "one_to_one", "paralog_broad": "paralog"}


def load_identity(sp: str, regime: str) -> dict[str, float]:
    sub = sp if regime == "one2one" else f"{sp}_paralog"
    p = DATA_V3 / sub / "ground_truth.json"
    d = json.load(open(p))
    return {pid: float(r["identity"]) for pid, r in d.items()
            if isinstance(r, dict) and "identity" in r}


def load_exp1() -> dict[tuple[str, str, str], dict[str, float]]:
    """(tool, sp, regime_canonical) → {protein_id: exon_f1}."""
    out: dict[tuple[str, str, str], dict[str, float]] = {}
    with open(EXP1_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            # Exp 1 regime label → canonical
            inv_map = {"one_to_one": "one2one", "paralog": "paralog_broad"}
            reg = inv_map.get(r["regime"])
            if reg is None:
                continue
            key = (r["tool"], r["species"], reg)
            out.setdefault(key, {})[r["protein_id"]] = float(r["exon_f1"])
    return out


def load_exp5(tool: str, sp: str, reg: str) -> dict[str, float]:
    """{protein_id: exon_f1, NaN if failed}."""
    path = EXP5_DIR / f"{tool}.{sp}.{reg}.50000.tsv"
    if not path.exists():
        return {}
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pid = r["protein_id"]
            if r["failed"] == "True":
                out[pid] = float("nan")
            else:
                try:
                    out[pid] = float(r["exon_f1"])
                except ValueError:
                    out[pid] = float("nan")
    return out


def load_b1(sp: str, reg: str) -> dict[str, dict[str, bool]]:
    """{protein_id: {tool: search_ok}} for 8 A-tools."""
    path = B1_DIR / f"per_protein_{sp}.{reg}.tsv"
    out: dict[str, dict[str, bool]] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            pid = r["protein_id"]
            out[pid] = {t: bool(int(r[f"A_{t}"])) for t in A_TOOLS}
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    exp1_idx = load_exp1()
    out_path = OUT_DIR / "per_protein.tsv"

    cols = (["protein_id", "species", "regime", "identity"]
            + [f"A_{t}" for t in A_TOOLS]
            + [f"C_{t}" for t in C_TOOLS]
            + [f"E2E_{t}" for t in E2E_TOOLS])

    stats = []  # for logging
    with open(out_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for sp in SPECIES:
            for reg in REGIMES:
                # Canonical protein set from GT (handles duplicates from paralog_broad)
                gt_path = GT_DIR / f"{sp}.{reg}.tsv"
                qids: list[str] = []
                seen: set[str] = set()
                with open(gt_path) as gfh:
                    for r in csv.DictReader(gfh, delimiter="\t"):
                        if r["query_id"] not in seen:
                            seen.add(r["query_id"])
                            qids.append(r["query_id"])

                ident = load_identity(sp, reg)
                a = load_b1(sp, reg)
                c_data = {t: load_exp5(t, sp, reg) for t in C_TOOLS}
                e2e_data = {t: exp1_idx.get((t, sp, reg), {}) for t in E2E_TOOLS}

                n_no_a = n_no_c = n_no_e2e = 0
                for pid in qids:
                    row = [pid, sp, reg]
                    # identity
                    iv = ident.get(pid)
                    row.append(f"{iv:.4f}" if iv is not None else "")
                    # A
                    if pid not in a:
                        n_no_a += 1
                        a_vals = ["0"] * len(A_TOOLS)
                    else:
                        a_vals = [str(int(a[pid][t])) for t in A_TOOLS]
                    row.extend(a_vals)
                    # C
                    c_vals = []
                    any_c_missing = False
                    for t in C_TOOLS:
                        v = c_data[t].get(pid)
                        if v is None:
                            any_c_missing = True
                            c_vals.append("")
                        elif v != v:  # NaN
                            c_vals.append("")
                        else:
                            c_vals.append(f"{v:.4f}")
                    if any_c_missing:
                        n_no_c += 1
                    row.extend(c_vals)
                    # E2E
                    e2e_vals = []
                    any_e2e_missing = False
                    for t in E2E_TOOLS:
                        v = e2e_data[t].get(pid)
                        if v is None:
                            any_e2e_missing = True
                            e2e_vals.append("")
                        else:
                            e2e_vals.append(f"{v:.4f}")
                    if any_e2e_missing:
                        n_no_e2e += 1
                    row.extend(e2e_vals)
                    fh.write("\t".join(row) + "\n")

                stats.append(dict(species=sp, regime=reg, n_rows=len(qids),
                                  n_missing_a=n_no_a, n_missing_c=n_no_c,
                                  n_missing_e2e=n_no_e2e))

    print(f"=== B2 joined table ===")
    total = sum(s["n_rows"] for s in stats)
    print(f"  wrote {total} rows × {len(cols)} cols to {out_path}")
    print(f"\n{'sp':<10} {'regime':<14} {'n':>5} {'no_A':>5} {'no_C':>5} {'no_E2E':>6}")
    for s in stats:
        print(f"  {s['species']:<10} {s['regime']:<14} {s['n_rows']:>5}"
              f" {s['n_missing_a']:>5} {s['n_missing_c']:>5} {s['n_missing_e2e']:>6}")

    # Schema doc
    schema_path = OUT_DIR / "schema.md"
    schema_path.write_text(f"""# Exp 6 joined per-protein table — column schema

Source: `results/exp6/joined/per_protein.tsv` ({total} rows × {len(cols)} cols)

## Key columns
- `protein_id`: Ensembl protein ID (query)
- `species`: mouse / dog / chicken / zebrafish
- `regime`: one2one / paralog_broad
- `identity`: percent identity to GT (from `data/v3/<sp>{{,_paralog}}/ground_truth.json`)

## A-stage search verdicts (8 tools, boolean 0/1)
`A_<tool>` for tool ∈ {{{', '.join(A_TOOLS)}}}.
Derived in B1 from `results/exp3/stage0_bed/<tool>.<sp>.<exp3_regime>.bed`
truncated to top-K=50, reciprocal overlap 0.5 vs GT.

## C-stage forced-locus exon F1 (7 tools, float ∈ [0, 1] or empty for failed/missing)
`C_<tool>` for tool ∈ {{{', '.join(C_TOOLS)}}}.
Sourced from `results/exp5/per_protein/<tool>.<sp>.<regime>.50000.tsv`,
`failed=True` rows → empty.

## E2E natural-pipeline exon F1 (6 e2e + 1 minisplice variant, float or empty)
`E2E_<tool>` for tool ∈ {{{', '.join(E2E_TOOLS)}}}.
Sourced from `results/exp1/per_protein.tsv`.

## NaN/empty semantics
- A: always 0 or 1 (Block 1 fills 0 for proteins missing from BED).
- C: empty cell ↔ kernel run failed (Exp 5 `failed=True`) or protein had no Exp 5 row.
  Downstream (B3) treats empty as F1=0 per Q3.1 default.
- E2E: empty cell ↔ no Exp 1 row for this (tool, sp, regime, protein). Treat as F1=0.
""")
    print(f"  schema → {schema_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
