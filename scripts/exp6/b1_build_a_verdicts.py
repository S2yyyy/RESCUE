"""Exp 6 Block 1 — Per-protein A-stage verdicts for 8 A-tools.

For each (A-tool, species, regime) cell, derive per-protein search_ok by
truncating Exp 3 stage0_bed to top-K=50 and applying the V4 §Stage 0
reciprocal overlap rule (0.5/0.5) against Exp 2 GT.

For canonical proteins missing from BED (high-A_T tools), search_ok=0
(evaluate_bed iterates the GT union and auto-fills).

Outputs:
  results/exp6/a_search/per_protein_<sp>.<regime>.tsv   wide: protein_id, A_<8 tools>
  results/exp6/a_search/per_cell.tsv                    aggregated hit-rate vs Exp 3 reference

Cross-check: per_cell hit-rate at K=50 should match results/exp3/recall/per_cell.tsv
(the K=50 row per cell) to within rounding. Disagreement > 0.005 raises a warning.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from scripts.stage0.eval_locus import (
    classify_query, load_candidates_by_query, load_gt_table,
)

REPO = Path(".").resolve()
GT_DIR = REPO / "results/exp2/gt"
BED_DIR = REPO / "results/exp3/stage0_bed"
EXP3_RECALL = REPO / "results/exp3/recall/per_cell.tsv"
OUT_DIR = REPO / "results/exp6/a_search"

SPECIES = ["mouse", "dog", "chicken", "zebrafish"]
REGIMES = ["one2one", "paralog_broad"]
A_TOOLS = ["miniprot", "spaln3", "genblastg", "tblastn", "diamond",
           "mmseqs", "last", "genblasta"]
EXP3_REGIME_MAP = {"one2one": "one2one", "paralog_broad": "paralog"}
TOP_K = 50
OVERLAP_THRESHOLD = 0.5


def derive_search_ok(bed_path: Path, gt_path: Path) -> dict[str, bool]:
    """Per-protein search_ok at K=50."""
    cands_by_query = load_candidates_by_query(str(bed_path))
    gts_by_query = load_gt_table(str(gt_path))
    out = {}
    all_q = sorted(set(cands_by_query) | set(gts_by_query))
    for qid in all_q:
        cands = cands_by_query.get(qid, [])[:TOP_K]  # top-50 truncation
        gts = gts_by_query.get(qid, [])
        v = classify_query(qid, cands, gts,
                           threshold=OVERLAP_THRESHOLD, enforce_strand=True)
        out[qid] = bool(v.search_ok)
    return out


def load_exp3_recall_k50() -> dict[tuple[str, str, str], float]:
    """K=50 recall per (tool, species, abc_regime) from results/exp3/recall/per_cell.tsv."""
    out = {}
    with open(EXP3_RECALL) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if int(r["k"]) != 50:
                continue
            out[(r["tool"], r["species"], r["abc_regime"])] = float(r["recall"])
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    exp3_recall = load_exp3_recall_k50()
    cell_rows = []

    for sp in SPECIES:
        for reg in REGIMES:
            gt_path = GT_DIR / f"{sp}.{reg}.tsv"
            # GT canonical protein set (from GT TSV, which has multiple rows per query)
            gt_qids = set()
            with open(gt_path) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    gt_qids.add(r["query_id"])
            n_canon = len(gt_qids)

            # Per-tool: derive search_ok at K=50
            tool_verdicts: dict[str, dict[str, bool]] = {}
            for tool in A_TOOLS:
                exp3_reg = EXP3_REGIME_MAP[reg]
                bed = BED_DIR / f"{tool}.{sp}.{exp3_reg}.bed"
                v = derive_search_ok(bed, gt_path)
                # Restrict to canonical proteins (drop any BED-only stragglers)
                v = {q: v.get(q, False) for q in gt_qids}
                tool_verdicts[tool] = v

                # Aggregate hit-rate
                n_ok = sum(1 for q in gt_qids if v[q])
                rate = n_ok / max(1, n_canon)
                ref = exp3_recall.get((tool, sp, reg))
                # Match Exp 3's reference: recall@50 should equal our hit-rate
                delta = (rate - ref) if ref is not None else float("nan")
                status = ("ok" if ref is None or abs(delta) < 0.005
                          else "mismatch")
                cell_rows.append(dict(
                    tool=tool, species=sp, regime=reg,
                    n_canon=n_canon, n_search_ok=n_ok,
                    hit_rate=rate, exp3_recall_k50=ref if ref is not None else float("nan"),
                    delta=delta, status=status,
                ))

            # Write wide per-protein TSV for this (sp, regime)
            out_path = OUT_DIR / f"per_protein_{sp}.{reg}.tsv"
            cols = ["protein_id"] + [f"A_{t}" for t in A_TOOLS]
            with open(out_path, "w") as fh:
                fh.write("\t".join(cols) + "\n")
                for qid in sorted(gt_qids):
                    vals = [qid] + [str(int(tool_verdicts[t][qid])) for t in A_TOOLS]
                    fh.write("\t".join(vals) + "\n")
            print(f"  wrote {out_path.name} ({n_canon} rows, 8 A-tools)")

    # Per-cell summary
    cell_path = OUT_DIR / "per_cell.tsv"
    cols = ["tool", "species", "regime", "n_canon", "n_search_ok",
            "hit_rate", "exp3_recall_k50", "delta", "status"]
    with open(cell_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in cell_rows:
            fh.write("\t".join(
                f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c])
                for c in cols) + "\n")

    # Summary print
    print(f"\n=== B1 verdict ===")
    n_mismatch = sum(1 for r in cell_rows if r["status"] == "mismatch")
    print(f"  Per-cell hit-rate vs Exp 3 recall@50: "
          f"{len(cell_rows) - n_mismatch}/{len(cell_rows)} cells agree |Δ|<0.005."
          f" {'PASS' if n_mismatch == 0 else 'FAIL'}")
    for r in cell_rows:
        if r["status"] == "mismatch":
            print(f"    MISMATCH: {r['tool']:<10} {r['species']:<10} {r['regime']:<14}"
                  f" hit_rate={r['hit_rate']:.4f} exp3={r['exp3_recall_k50']:.4f}"
                  f" Δ={r['delta']:+.4f}")

    # Sanity: A-stage union-recall (any tool finds GT)
    print(f"\n=== Union-recall (any A-tool finds GT) per (sp, regime) ===")
    for sp in SPECIES:
        for reg in REGIMES:
            path = OUT_DIR / f"per_protein_{sp}.{reg}.tsv"
            with open(path) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                n_total = 0
                n_any = 0
                for r in reader:
                    n_total += 1
                    if any(int(r[f"A_{t}"]) for t in A_TOOLS):
                        n_any += 1
            print(f"  {sp:<10} {reg:<14} union={n_any}/{n_total} = {n_any/n_total:.4f}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
