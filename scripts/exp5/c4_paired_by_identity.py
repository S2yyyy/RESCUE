"""C4 paired-by-identity analysis.

Disambiguates the raw C4 FAIL (paralog forced-locus F1 drops 0.22-0.31 vs 1:1
forced) into "identity confound" vs "real kernel-regime sensitivity".

Why this exists: paralog protein sets have systematically lower mean identity
(20-30 points lower than 1:1 set) because Compara's broad-ortholog match pulls
in distant paralogs. Identity is the strongest C-stage F1 predictor, so raw Δ
between regimes mixes "regime effect" with "identity composition effect". This
script stratifies both regimes into identity bins and re-computes Δ within
each bin.

Reads:
  results/exp5/per_protein/<tool>.<sp>.<regime>.50000.tsv  (per-protein F1)
  data/v3/<sp>/ground_truth.json                            (1:1 identity)
  data/v3/<sp>_paralog/ground_truth.json                    (paralog identity)

Writes (stdout, and TSV if --out is given):
  Per (tool, species, identity_bin): n1:1, nPar, F1(1:1), F1(par), Δ.
  Plus a verdict table per tool at the >=90 bin (cleanest paired comparison).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

SPECIES = ["mouse", "dog", "chicken", "zebrafish"]
TOOLS = ["miniprot", "spaln3", "gth", "exonerate", "genewise", "scipio", "genblastg"]
BINS = [(">=90", lambda x: x >= 90),
        ("70-90", lambda x: 70 <= x < 90),
        ("50-70", lambda x: 50 <= x < 70)]


def load_identity(repo: Path, sp: str, regime: str) -> dict[str, float]:
    sub = sp if regime == "one2one" else f"{sp}_paralog"
    p = repo / "data" / "v3" / sub / "ground_truth.json"
    d = json.load(open(p))
    return {pid: float(r["identity"]) for pid, r in d.items()
            if isinstance(r, dict) and "identity" in r}


def load_per_protein(repo: Path, tool: str, sp: str, regime: str, pad: int):
    p = repo / "results/exp5/per_protein" / f"{tool}.{sp}.{regime}.{pad}.tsv"
    if not p.exists():
        return None
    return list(csv.DictReader(open(p), delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default=".", help="Project root (default: cwd)")
    ap.add_argument("--padding", type=int, default=50000)
    ap.add_argument("--out", default=None,
                    help="Optional TSV path to write the long-format table")
    args = ap.parse_args()
    repo = Path(args.repo).resolve()

    results = []
    for tool in TOOLS:
        for sp in SPECIES:
            ident_o = load_identity(repo, sp, "one2one")
            ident_p = load_identity(repo, sp, "paralog_broad")
            ror = load_per_protein(repo, tool, sp, "one2one", args.padding)
            rpr = load_per_protein(repo, tool, sp, "paralog_broad", args.padding)
            if ror is None or rpr is None:
                continue
            for bin_name, pred in BINS:
                f1o = [float(r["exon_f1"]) for r in ror
                       if r["failed"] != "True" and r["protein_id"] in ident_o
                       and pred(ident_o[r["protein_id"]])]
                f1p = [float(r["exon_f1"]) for r in rpr
                       if r["failed"] != "True" and r["protein_id"] in ident_p
                       and pred(ident_p[r["protein_id"]])]
                mo = sum(f1o) / len(f1o) if f1o else float("nan")
                mp = sum(f1p) / len(f1p) if f1p else float("nan")
                delta = mo - mp if f1o and f1p else float("nan")
                results.append({
                    "tool": tool, "species": sp, "bin": bin_name,
                    "n_one2one": len(f1o), "n_paralog": len(f1p),
                    "mean_f1_one2one": mo, "mean_f1_paralog": mp,
                    "delta": delta,
                })

    # Stdout: long format
    print(f"{'tool':<11} {'sp':<10} {'bin':<6} {'n1:1':>5} {'nPar':>5}"
          f" {'F1(1:1)':>8} {'F1(par)':>8} {'Δ':>8}")
    print("-" * 70)
    prev = None
    for r in results:
        if r["tool"] != prev and prev is not None:
            print()
        prev = r["tool"]
        d = r["delta"]
        d_s = f"{d:+.3f}" if d == d else "  N/A"
        mo_s = f"{r['mean_f1_one2one']:.3f}" if r["mean_f1_one2one"] == r["mean_f1_one2one"] else "N/A"
        mp_s = f"{r['mean_f1_paralog']:.3f}" if r["mean_f1_paralog"] == r["mean_f1_paralog"] else "N/A"
        flag = " PASS" if (d == d and d <= 0.05) else ""
        print(f"{r['tool']:<11} {r['species']:<10} {r['bin']:<6}"
              f" {r['n_one2one']:>5} {r['n_paralog']:>5}"
              f" {mo_s:>8} {mp_s:>8} {d_s:>8}{flag}")

    # Stdout: verdict per tool at >=90 bin
    print("\n=== C4 verdict per tool (identity-controlled at >=90 bin) ===")
    for tool in TOOLS:
        ds = [r["delta"] for r in results
              if r["tool"] == tool and r["bin"] == ">=90" and r["delta"] == r["delta"]]
        if not ds:
            continue
        mean_d = sum(ds) / len(ds)
        verdict = "PASS" if mean_d <= 0.05 else "FAIL"
        print(f"  {tool:<11} mean Δ (>=90, n_sp={len(ds)}) = {mean_d:+.3f}  {verdict}")

    # Optional TSV
    if args.out:
        cols = ["tool", "species", "bin", "n_one2one", "n_paralog",
                "mean_f1_one2one", "mean_f1_paralog", "delta"]
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in results:
                fh.write("\t".join(
                    f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k])
                    for k in cols) + "\n")
        print(f"\nWrote {len(results)} rows to {args.out}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
