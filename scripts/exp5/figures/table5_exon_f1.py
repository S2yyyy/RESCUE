"""Table 5 — per-tool exon F1 + frac_failed at 50 kb padding.

Inputs:  results/exp5/cell/per_cell.tsv
Outputs:
  results/exp5/figures/table5_exon_f1.csv
  results/exp5/figures/table5_exon_f1.txt   (paper-ready pretty print)
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from rescue_paths import PROJECT_ROOT  # noqa: E402

REPO = PROJECT_ROOT
IN = REPO / "results/exp5/cell/per_cell.tsv"
OUT_DIR = REPO / "results/exp5/figures"

TOOLS = ["miniprot", "spaln3", "gth", "exonerate",
         "genewise", "scipio", "genblastg", "minisplice_miniprot"]
SPECIES = ["mouse", "dog", "chicken", "zebrafish"]


def main() -> int:
    rows = list(csv.DictReader(open(IN), delimiter="\t"))
    rows = [r for r in rows
            if r["regime"] == "one2one" and int(r["padding"]) == 50000]
    idx = {(r["tool"], r["species"]): r for r in rows}

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_out = OUT_DIR / "table5_exon_f1.csv"
    with csv_out.open("w") as fh:
        fh.write("tool,species,n_total,n_failed,frac_failed,"
                 "mean_exon_f1,ci_lo,ci_hi,"
                 "mean_junction_f1,mean_cds_recovery,mean_full_length\n")
        for t in TOOLS:
            for sp in SPECIES:
                r = idx.get((t, sp))
                if r is None:
                    continue
                fh.write(f"{t},{sp},{r['n_total']},{r['n_failed']},"
                         f"{float(r['frac_failed']):.4f},"
                         f"{float(r['mean_exon_f1']):.4f},"
                         f"{float(r['ci_lo']):.4f},"
                         f"{float(r['ci_hi']):.4f},"
                         f"{float(r['mean_junction_f1']):.4f},"
                         f"{float(r['mean_cds_recovery']):.4f},"
                         f"{float(r['mean_full_length']):.4f}\n")

    txt_out = OUT_DIR / "table5_exon_f1.txt"
    with txt_out.open("w") as fh:
        fh.write("Table 5 — Forced-locus exon F1 @ 50 kb padding\n")
        fh.write("Format: mean_F1 [frac_failed] (n_failed/n_total)\n\n")
        fh.write(f"{'tool':<22} ")
        for sp in SPECIES:
            fh.write(f"{sp:>22}")
        fh.write("\n")
        for t in TOOLS:
            fh.write(f"{t:<22} ")
            for sp in SPECIES:
                r = idx.get((t, sp))
                if r is None:
                    fh.write(f"{'--':>22}")
                    continue
                f1 = float(r["mean_exon_f1"])
                ff = float(r["frac_failed"])
                ent = f"{f1:.3f} [{ff:.3f}] ({r['n_failed']}/{r['n_total']})"
                fh.write(f"{ent:>22}")
            fh.write("\n")
    print(open(txt_out).read())
    print(f"wrote {csv_out}")
    print(f"wrote {txt_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
