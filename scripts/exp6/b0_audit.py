"""Exp 6 Block 0 — Pre-flight audit.

Three sub-audits before any oracle computation:

1. Protein-ID inventory per (species, regime) from Exp 2 GT tables.
2. Per-source coverage: how many of the canonical proteins are present in
   each source × tool × cell combination (Exp 1 / Exp 2 / Exp 3 / Exp 5)?
3. search_ok consistency for miniprot and spaln3 — re-derive search_ok from
   Exp 3 stage0_bed via the V4 §Stage 0 protocol (eval_locus.evaluate_bed)
   and compare per-protein to Exp 2's stage0_verdicts/<tool>.<sp>.<regime>.tsv.

Outputs:
  results/exp6/audit/inventory.tsv     — n_proteins per (species, regime)
  results/exp6/audit/coverage.tsv      — long format, one row per (source, tool, sp, regime)
  results/exp6/audit/consistency.tsv   — miniprot/spaln3 agreement rates

Acceptance gates (see EXPERIMENT_PLAN_EXP6.md §Block 0):
  - For exp1/exp2/exp5 (per-protein row sources): pct_coverage ≥ 99.5%.
    These sources always emit a row per canonical protein (even failures/0-hits).
  - For exp3_bed (candidate-record source): BED file exists and contains ≥1 record.
    pct_coverage = "fraction of proteins with ≥1 candidate" is reported as an
    informational measure of tool-level hit-rate, NOT a data-quality gate
    (Block 1's evaluate_bed iterates GT and emits search_ok=0 for missing proteins).
  - miniprot/spaln3 search_ok agreement ≥ 99% per cell.
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

from scripts.stage0.eval_locus import evaluate_bed, load_gt_table

REPO = Path(".").resolve()
GT_DIR = REPO / "results/exp2/gt"
EXP1_PER_PROTEIN = REPO / "results/exp1/per_protein.tsv"
EXP2_ABC = REPO / "results/exp2/abc/per_query_abc.tsv"
EXP2_VERDICTS = REPO / "results/exp2/abc/stage0_verdicts"
EXP3_BED_DIR = REPO / "results/exp3/stage0_bed"
EXP5_DIR = REPO / "results/exp5/per_protein"
OUT_DIR = REPO / "results/exp6/audit"

SPECIES = ["mouse", "dog", "chicken", "zebrafish"]
REGIMES = ["one2one", "paralog_broad"]  # Exp 6 canonical regimes

# Exp 1 uses regime label `one_to_one` / `paralog` — map to canonical.
EXP1_REGIME_MAP = {"one2one": "one_to_one", "paralog_broad": "paralog"}
# Exp 3 stage0_bed uses `one2one` / `paralog` — map paralog_broad → paralog
EXP3_REGIME_MAP = {"one2one": "one2one", "paralog_broad": "paralog"}

E2E_TOOLS = ["miniprot", "spaln3", "gth", "scipio", "genblastg", "minisplice_miniprot"]
# Tools we expect to find search_ok for in Exp 2 (subset of E2E_TOOLS)
EXP2_TOOLS = ["miniprot", "spaln3", "gth", "scipio", "genblastg", "minisplice_miniprot"]
# A-stage panel for Exp 6 (8 tools, BLAT excluded per EXCLUSION_blat.md)
A_TOOLS = ["miniprot", "spaln3", "genblastg", "tblastn", "diamond",
           "mmseqs", "last", "genblasta"]
# C-stage panel for Exp 6 (7 tools, minisplice_miniprot deferred to B8)
C_TOOLS = ["miniprot", "spaln3", "gth", "exonerate", "genewise", "scipio", "genblastg"]


def audit_inventory() -> dict[tuple[str, str], int]:
    """Task 1: canonical protein count per (sp, regime) from GT TSV."""
    inv = {}
    for sp in SPECIES:
        for reg in REGIMES:
            path = GT_DIR / f"{sp}.{reg}.tsv"
            qids = set()
            with open(path) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    qids.add(row["query_id"])
            inv[(sp, reg)] = len(qids)
    return inv


def audit_coverage(inv: dict[tuple[str, str], int]) -> list[dict]:
    """Task 2: per-source coverage check."""
    rows = []

    # Build Exp 1 index: (tool, sp, regime) → set of protein_ids present
    exp1_idx: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    with open(EXP1_PER_PROTEIN) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            exp1_idx[(r["tool"], r["species"], r["regime"])].add(r["protein_id"])

    # Build Exp 2 index from per_query_abc
    exp2_idx: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    with open(EXP2_ABC) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            exp2_idx[(r["tool"], r["species"], r["abc_regime"])].add(r["protein_id"])

    # Helper to count proteins in a stage0_bed (column 4 = protein_id)
    def count_bed_proteins(path: Path) -> set[str]:
        qids = set()
        if not path.exists():
            return qids
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) >= 4:
                    qids.add(f[3])
        return qids

    # Helper to count proteins in Exp 5 per_protein TSV
    def count_exp5(path: Path) -> set[str]:
        qids = set()
        if not path.exists():
            return qids
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                qids.add(r["protein_id"])
        return qids

    # Iterate (source, tool, sp, regime)
    for sp in SPECIES:
        for reg in REGIMES:
            n_exp = inv[(sp, reg)]

            # Exp 1 baseline: 6 e2e tools (gth included even if shared search)
            for tool in EXP2_TOOLS:
                exp1_reg = EXP1_REGIME_MAP[reg]
                found = exp1_idx.get((tool, sp, exp1_reg), set())
                rows.append(dict(source="exp1", tool=tool, species=sp, regime=reg,
                                 n_expected=n_exp, n_found=len(found),
                                 pct_coverage=len(found) / max(1, n_exp)))

            # Exp 2 search_ok
            for tool in EXP2_TOOLS:
                found = exp2_idx.get((tool, sp, reg), set())
                rows.append(dict(source="exp2", tool=tool, species=sp, regime=reg,
                                 n_expected=n_exp, n_found=len(found),
                                 pct_coverage=len(found) / max(1, n_exp)))

            # Exp 3 stage0_bed: 8 A-tools
            for tool in A_TOOLS:
                exp3_reg = EXP3_REGIME_MAP[reg]
                bed = EXP3_BED_DIR / f"{tool}.{sp}.{exp3_reg}.bed"
                found = count_bed_proteins(bed)
                rows.append(dict(source="exp3_bed", tool=tool, species=sp, regime=reg,
                                 n_expected=n_exp, n_found=len(found),
                                 pct_coverage=len(found) / max(1, n_exp)))

            # Exp 5 forced-F1: 7 C-tools at pad=50000
            for tool in C_TOOLS:
                tsv = EXP5_DIR / f"{tool}.{sp}.{reg}.50000.tsv"
                found = count_exp5(tsv)
                rows.append(dict(source="exp5", tool=tool, species=sp, regime=reg,
                                 n_expected=n_exp, n_found=len(found),
                                 pct_coverage=len(found) / max(1, n_exp)))

    return rows


def audit_consistency() -> list[dict]:
    """Task 3: miniprot/spaln3 search_ok agreement between Exp 3 stage0_bed
    (re-derived via eval_locus) and Exp 2's stage0_verdicts/."""
    rows = []
    for tool in ["miniprot", "spaln3"]:
        for sp in SPECIES:
            for reg in REGIMES:
                exp3_reg = EXP3_REGIME_MAP[reg]
                bed = EXP3_BED_DIR / f"{tool}.{sp}.{exp3_reg}.bed"
                gt = GT_DIR / f"{sp}.{reg}.tsv"
                verdict_tsv = EXP2_VERDICTS / f"{tool}.{sp}.{reg}.tsv"

                if not bed.exists() or not gt.exists() or not verdict_tsv.exists():
                    rows.append(dict(tool=tool, species=sp, regime=reg,
                                     status="missing_input",
                                     n_compared=0, n_agree=0, agreement=float("nan"),
                                     n_disagree=0, exp3_extra=0, exp2_extra=0))
                    continue

                # Re-derive search_ok from Exp 3 BED + Exp 2 GT
                stage0_verdicts = {v.query_id: bool(v.search_ok)
                                   for v in evaluate_bed(str(bed), str(gt))}

                # Load Exp 2's canonical search_ok
                exp2_search_ok = {}
                with open(verdict_tsv) as fh:
                    for r in csv.DictReader(fh, delimiter="\t"):
                        exp2_search_ok[r["query_id"]] = bool(int(r["search_ok"]))

                common = set(stage0_verdicts) & set(exp2_search_ok)
                n_agree = sum(1 for q in common
                              if stage0_verdicts[q] == exp2_search_ok[q])
                n_disagree = len(common) - n_agree
                rows.append(dict(
                    tool=tool, species=sp, regime=reg,
                    status="ok",
                    n_compared=len(common),
                    n_agree=n_agree,
                    n_disagree=n_disagree,
                    agreement=n_agree / max(1, len(common)),
                    exp3_extra=len(set(stage0_verdicts) - set(exp2_search_ok)),
                    exp2_extra=len(set(exp2_search_ok) - set(stage0_verdicts)),
                ))
    return rows


def write_tsv(rows: list[dict], path: Path, columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(
                f"{r[c]:.4f}" if isinstance(r.get(c), float) else str(r.get(c, ""))
                for c in columns) + "\n")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Task 1: protein-ID inventory ===")
    inv = audit_inventory()
    inv_rows = [dict(species=sp, regime=reg, n_proteins=n)
                for (sp, reg), n in sorted(inv.items())]
    write_tsv(inv_rows, OUT_DIR / "inventory.tsv",
              ["species", "regime", "n_proteins"])
    for r in inv_rows:
        print(f"  {r['species']:<10} {r['regime']:<14} n={r['n_proteins']}")

    print("\n=== Task 2: per-source coverage ===")
    cov = audit_coverage(inv)
    write_tsv(cov, OUT_DIR / "coverage.tsv",
              ["source", "tool", "species", "regime",
               "n_expected", "n_found", "pct_coverage"])

    # Summary table: per source, count of cells passing the appropriate gate.
    # exp1/exp2/exp5 sources MUST have ≥99.5% (rows-per-protein semantics).
    # exp3_bed source: only check ≥1 record (gate); pct_coverage is informational.
    by_source = defaultdict(list)
    for r in cov:
        by_source[r["source"]].append(r)

    for src in ["exp1", "exp2", "exp5"]:
        rs = by_source[src]
        n_cells = len(rs)
        n_pass = sum(1 for r in rs if r["pct_coverage"] >= 0.995)
        n_low = [r for r in rs if r["pct_coverage"] < 0.995]
        print(f"  {src:<9} {n_pass}/{n_cells} cells pass ≥99.5% (row-per-protein gate)")
        for r in sorted(n_low, key=lambda x: x["pct_coverage"])[:5]:
            print(f"    LOW: {r['tool']:<14} {r['species']:<10} {r['regime']:<14}"
                  f" n_found={r['n_found']}/{r['n_expected']}"
                  f" ({r['pct_coverage']:.2%})")

    # exp3_bed gate: BED file present and ≥1 candidate
    rs3 = by_source["exp3_bed"]
    n_cells3 = len(rs3)
    n_present = sum(1 for r in rs3 if r["n_found"] >= 1)
    pcts = sorted(r["pct_coverage"] for r in rs3)
    print(f"  exp3_bed  {n_present}/{n_cells3} cells have ≥1 candidate (gate)")
    print(f"            pct_with_candidates: min={pcts[0]:.1%},"
          f" median={pcts[len(pcts)//2]:.1%}, max={pcts[-1]:.1%}")
    print(f"            (informational only — Block 1's evaluate_bed"
          f" auto-fills search_ok=0 for proteins missing from BED)")
    missing3 = [r for r in rs3 if r["n_found"] == 0]
    for r in missing3:
        print(f"    EMPTY BED: {r['tool']:<14} {r['species']:<10} {r['regime']:<14}")

    print("\n=== Task 3: miniprot/spaln3 search_ok consistency (Exp 3 BED → Exp 2 verdict) ===")
    con = audit_consistency()
    write_tsv(con, OUT_DIR / "consistency.tsv",
              ["tool", "species", "regime", "status", "n_compared",
               "n_agree", "n_disagree", "agreement",
               "exp3_extra", "exp2_extra"])

    print(f"  {'tool':<10} {'sp':<10} {'regime':<14}"
          f" {'n_cmp':>6} {'n_agree':>8} {'agree':>8} {'exp3+':>6} {'exp2+':>6}")
    n_below = 0
    for r in con:
        if r["status"] != "ok":
            print(f"  {r['tool']:<10} {r['species']:<10} {r['regime']:<14}"
                  f" {r['status']}")
            continue
        flag = " ⚠️" if r["agreement"] < 0.99 else ""
        if r["agreement"] < 0.99:
            n_below += 1
        print(f"  {r['tool']:<10} {r['species']:<10} {r['regime']:<14}"
              f" {r['n_compared']:>6} {r['n_agree']:>8}"
              f" {r['agreement']:>7.3%} {r['exp3_extra']:>6} {r['exp2_extra']:>6}{flag}")

    print(f"\n=== Verdict ===")
    # Gate 1: exp1/exp2/exp5 row-per-protein cells must be ≥99.5%
    row_sources = [r for r in cov if r["source"] in ("exp1", "exp2", "exp5")]
    row_low = [r for r in row_sources if r["pct_coverage"] < 0.995]
    print(f"  Row-per-protein coverage (exp1/exp2/exp5):"
          f" {len(row_sources) - len(row_low)}/{len(row_sources)} cells ≥ 99.5%."
          f" {'PASS' if not row_low else 'FAIL'}")

    # Gate 2: exp3_bed cells must have ≥1 record
    bed_sources = [r for r in cov if r["source"] == "exp3_bed"]
    bed_empty = [r for r in bed_sources if r["n_found"] == 0]
    print(f"  Exp 3 BED presence:"
          f" {len(bed_sources) - len(bed_empty)}/{len(bed_sources)} cells have ≥1 candidate."
          f" {'PASS' if not bed_empty else 'FAIL'}")

    # Gate 3: miniprot/spaln3 consistency ≥99%
    print(f"  miniprot/spaln3 search_ok consistency:"
          f" {len(con) - n_below}/{len(con)} cells ≥ 99% agreement."
          f" {'PASS' if n_below == 0 else 'FAIL'}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
