"""B5 — failure-rate stratification.

Tests anti-claim "exon F1 differences are an artefact of failure rates".

Compute mean exon F1 on:
  (a) "everyone-succeeds" set = proteins where all 7 tools produced output
  (b) full union (per-tool only on its own successes)

If tool ranking is preserved (Spearman ρ ≥ 0.85), the anti-claim is
defended. Otherwise the §3.3.3 narrative needs a failure-rate caveat.

Outputs:
  results/exp5/cell/failure_stratified.tsv
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import PROJECT_ROOT  # noqa: E402

REPO = PROJECT_ROOT
EXP5_DIR = REPO / "results/exp5/per_protein"

TOOLS = ["miniprot", "spaln3", "gth", "exonerate", "genewise",
         "scipio", "genblastg", "minisplice_miniprot"]
SPECIES = ["mouse", "dog", "chicken", "zebrafish"]


def spearman_rho(a, b):
    if len(a) != len(b) or len(a) < 2:
        return float("nan")
    def rank(xs):
        s = sorted(range(len(xs)), key=lambda i: xs[i])
        ranks = [0.0] * len(xs)
        i = 0
        while i < len(s):
            j = i
            while j + 1 < len(s) and xs[s[j+1]] == xs[s[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j+1):
                ranks[s[k]] = avg
            i = j + 1
        return ranks
    ra = rank(list(a)); rb = rank(list(b))
    n = len(ra)
    ma = sum(ra) / n; mb = sum(rb) / n
    num = sum((ra[i]-ma)*(rb[i]-mb) for i in range(n))
    den = ((sum((ra[i]-ma)**2 for i in range(n))**0.5) *
           (sum((rb[i]-mb)**2 for i in range(n))**0.5))
    return num / den if den else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--padding", type=int, default=50000)
    ap.add_argument("--out", default=str(REPO / "results/exp5/cell/failure_stratified.tsv"))
    args = ap.parse_args()

    rows_out = []
    for sp in SPECIES:
        # Load all tool rows, keyed by protein_id
        per_tool: dict[str, dict[str, dict]] = defaultdict(dict)
        for tool in TOOLS:
            tsv = EXP5_DIR / f"{tool}.{sp}.one2one.{args.padding}.tsv"
            if not tsv.exists():
                continue
            with tsv.open() as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    per_tool[tool][r["protein_id"]] = r
        if not per_tool:
            continue

        # Everyone-succeeds: protein IDs where all tools have failed=False
        all_pids = set.intersection(*[set(per_tool[t]) for t in per_tool])
        succ_pids = set()
        for pid in all_pids:
            if all(per_tool[t].get(pid, {}).get("failed") == "False"
                   for t in per_tool):
                succ_pids.add(pid)

        full_means = {}
        ev_means = {}
        for tool in per_tool:
            full_vals = [float(r["exon_f1"]) for r in per_tool[tool].values()
                         if r.get("failed") == "False" and r.get("exon_f1")]
            ev_vals = [float(per_tool[tool][p]["exon_f1"])
                       for p in succ_pids
                       if p in per_tool[tool]
                       and per_tool[tool][p].get("exon_f1")]
            full_means[tool] = (sum(full_vals)/len(full_vals)) if full_vals else float("nan")
            ev_means[tool] = (sum(ev_vals)/len(ev_vals)) if ev_vals else float("nan")
            rows_out.append({
                "species": sp, "tool": tool,
                "n_full": len(full_vals), "mean_full": full_means[tool],
                "n_ev": len(ev_vals), "mean_ev": ev_means[tool],
            })

        tool_keys = sorted(full_means)
        a = [full_means[t] for t in tool_keys
             if full_means[t] == full_means[t] and ev_means[t] == ev_means[t]]
        b = [ev_means[t] for t in tool_keys
             if full_means[t] == full_means[t] and ev_means[t] == ev_means[t]]
        rho = spearman_rho(a, b)
        rows_out.append({
            "species": sp, "tool": "_RHO_full_vs_ev",
            "n_full": "", "mean_full": "",
            "n_ev": len(succ_pids), "mean_ev": rho,
        })

    cols = ["species", "tool", "n_full", "mean_full", "n_ev", "mean_ev"]
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows_out:
            fh.write("\t".join(
                f"{r[c]:.4f}" if isinstance(r.get(c), float) else str(r.get(c, ""))
                for c in cols) + "\n")
    print(f"wrote {len(rows_out)} rows to {args.out}")
    rho_rows = [r for r in rows_out if r["tool"] == "_RHO_full_vs_ev"]
    n_pass = sum(1 for r in rho_rows if r["mean_ev"] != "" and r["mean_ev"] >= 0.85)
    print(f"B5 acceptance (rho ≥ 0.85 on every species): "
          f"{n_pass}/{len(rho_rows)} pass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
