"""B6 — identity-stratified exon F1.

Bin queries by V3 percent identity (≥90, 70-90, <70). Per-tool exon F1
per bin. Detect rank-flips ≥2 across bins on ≥2 species.

Inputs:
  results/exp5/per_protein/<tool>.<sp>.one2one.50000.tsv
  data/v3/<sp>/ground_truth.json (identity field)

Outputs:
  results/exp5/cell/identity_stratified.tsv
  results/exp5/cell/identity_rankflips.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import PROJECT_ROOT, V3_DATA  # noqa: E402

REPO = PROJECT_ROOT
EXP5_DIR = REPO / "results/exp5/per_protein"
V3 = V3_DATA

BINS = [(">=90", lambda x: x >= 90),
        ("70-90", lambda x: 70 <= x < 90),
        ("<70", lambda x: x < 70)]
SPECIES = ["mouse", "dog", "chicken", "zebrafish"]
TOOLS = ["miniprot", "spaln3", "gth", "exonerate", "genewise",
         "scipio", "genblastg", "minisplice_miniprot"]


def load_identity(species: str) -> dict[str, float]:
    p = V3 / species / "ground_truth.json"
    if not p.exists():
        return {}
    d = json.load(open(p))
    return {pid: float(r["identity"]) for pid, r in d.items()
            if isinstance(r, dict) and "identity" in r}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--padding", type=int, default=50000)
    ap.add_argument("--out-stratified",
                    default=str(REPO / "results/exp5/cell/identity_stratified.tsv"))
    ap.add_argument("--out-rankflips",
                    default=str(REPO / "results/exp5/cell/identity_rankflips.tsv"))
    args = ap.parse_args()

    rows = []
    for tool in TOOLS:
        for sp in SPECIES:
            tsv = EXP5_DIR / f"{tool}.{sp}.one2one.{args.padding}.tsv"
            if not tsv.exists():
                continue
            ident = load_identity(sp)
            buckets = defaultdict(list)
            with tsv.open() as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    if r.get("failed") == "True":
                        continue
                    pid = r["protein_id"]
                    iv = ident.get(pid)
                    if iv is None:
                        continue
                    bin_label = next((lbl for lbl, fn in BINS if fn(iv)), None)
                    if bin_label is None:
                        continue
                    try:
                        f1 = float(r["exon_f1"])
                    except ValueError:
                        continue
                    buckets[bin_label].append(f1)
            for bin_label, _ in BINS:
                vals = buckets[bin_label]
                rows.append({
                    "tool": tool, "species": sp, "identity_bin": bin_label,
                    "n": len(vals),
                    "mean_exon_f1": (sum(vals) / len(vals)) if vals else float("nan"),
                    "low_n": str(len(vals) < 30).lower(),
                })

    cols = ["tool", "species", "identity_bin", "n", "mean_exon_f1", "low_n"]
    Path(args.out_stratified).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_stratified, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                f"{r[c]:.4f}" if isinstance(r.get(c), float) else str(r.get(c, ""))
                for c in cols) + "\n")

    # Rank-flip detection
    flip_rows = []
    n_flip = defaultdict(int)
    for sp in SPECIES:
        tool_ranks_per_bin = {}
        for bin_label, _ in BINS:
            cell_rows = [r for r in rows
                         if r["species"] == sp and r["identity_bin"] == bin_label
                         and r["mean_exon_f1"] == r["mean_exon_f1"]
                         and r.get("low_n") == "false"]
            cell_rows.sort(key=lambda r: -r["mean_exon_f1"])
            tool_ranks_per_bin[bin_label] = {
                r["tool"]: i + 1 for i, r in enumerate(cell_rows)
            }
        for tool in TOOLS:
            ranks = [tool_ranks_per_bin.get(b, {}).get(tool) for b, _ in BINS]
            present = [r for r in ranks if r is not None]
            if len(present) < 2:
                continue
            md = max(present) - min(present)
            flip_rows.append({
                "tool": tool, "species": sp,
                "rank_>=90": ranks[0] or "",
                "rank_70-90": ranks[1] or "",
                "rank_<70": ranks[2] or "",
                "max_delta_rank": md,
            })
            if md >= 2:
                n_flip[tool] += 1

    cols_f = ["tool", "species", "rank_>=90", "rank_70-90", "rank_<70",
              "max_delta_rank"]
    with open(args.out_rankflips, "w") as fh:
        fh.write("\t".join(cols_f) + "\n")
        for r in flip_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols_f) + "\n")

    print(f"wrote {len(rows)} stratified rows to {args.out_stratified}")
    print(f"wrote {len(flip_rows)} flip rows to {args.out_rankflips}")
    print()
    print("Tools with rank-flip ≥ 2 across identity bins:")
    for tool, n_sp in sorted(n_flip.items()):
        print(f"  {tool}: {n_sp} species")
    n_pass = sum(1 for v in n_flip.values() if v >= 2)
    print(f"\nB6 acceptance (≥1 tool flips on ≥2 species): "
          f"{'PASS' if n_pass >= 1 else 'FAIL'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
