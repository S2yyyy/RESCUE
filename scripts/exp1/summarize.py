#!/usr/bin/env python3
"""Aggregate all Exp 1 per-batch summaries into one TSV + long-form per-protein TSV."""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import PROJECT_ROOT, RESULTS_ROOT  # noqa: E402

ROOT = PROJECT_ROOT
sys.path.insert(0, str(ROOT / "scripts/exp1"))

import yaml  # noqa: E402
from common import out_path  # noqa: E402
from metrics import compute, METRICS_DIR  # noqa: E402

CFG = yaml.safe_load((ROOT / "configs/exp1_tools.yaml").read_text())
OUT = RESULTS_ROOT / "exp1"


def main():
    rows = []
    long_rows = []
    for tool in CFG["tools"]:
        for species in CFG["species"]:
            for regime in CFG["regimes"]:
                suffix = CFG["regimes"][regime]
                raw_json = out_path(tool, species, suffix)
                if not raw_json.exists():
                    continue
                try:
                    s = compute(tool, species, regime)
                except Exception as e:
                    print(f"metrics failed for {tool}.{species}.{regime}: {e}",
                          file=sys.stderr)
                    continue
                rows.append(s)
                pp = json.loads((METRICS_DIR / f"{tool}.{species}.{regime}.per_protein.json").read_text())
                for pid, ev in pp.items():
                    long_rows.append({
                        "tool": tool, "species": species, "regime": regime,
                        "protein_id": pid,
                        "exon_f1": ev.get("exon_f1", 0),
                        "junction_f1": ev.get("junction_f1", 0),
                        "cds_recovery": ev.get("cds_recovery", 0),
                        "full_length": int(bool(ev.get("full_length"))),
                        "correct_locus": int(bool(ev.get("correct_locus"))),
                        "no_hit": int(bool(ev.get("no_hit"))),
                        "strict_e1": int(bool(ev.get("e1_strict_e1", False))),
                    })

    if not rows:
        print("No results found yet.", file=sys.stderr)
        return

    keys = sorted({k for r in rows for k in r.keys()})
    with (OUT / "summary.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    lkeys = list(long_rows[0].keys())
    with (OUT / "per_protein.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=lkeys, delimiter="\t")
        w.writeheader()
        for r in long_rows:
            w.writerow(r)

    print(f"Wrote {OUT / 'summary.tsv'}  ({len(rows)} rows)")
    print(f"Wrote {OUT / 'per_protein.tsv'}  ({len(long_rows)} rows)")


if __name__ == "__main__":
    main()
