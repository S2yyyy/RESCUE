#!/usr/bin/env python3
"""Compute Exp 1 metrics for one (tool, species, regime) result.

Failures / missing predictions are scored as 0 (failure_policy: zero).
Writes:
  results/exp1/metrics/{tool}.{species}.{regime}.per_protein.json
  results/exp1/metrics/{tool}.{species}.{regime}.summary.json
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import PROJECT_ROOT, RESULTS_ROOT  # noqa: E402

ROOT = PROJECT_ROOT
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "scripts/exp1"))

import yaml  # noqa: E402
from v3_tool_common import load_ground_truth, load_proteins  # noqa: E402
from v3_evaluate import evaluate_hit  # noqa: E402
from v3_fix_paralog_scoring import (  # noqa: E402
    load_cds_bed, get_gene_locus, classify_e1,
)


from common import out_path  # noqa: E402


def get_gene_exons(gene_cds_dict, gene_id):
    """Get CDS exon coords for a gene from CDS BED dict.

    Returns list of {chrom, start, end} dicts (1-based inclusive).
    """
    cds_list = gene_cds_dict.get(gene_id, [])
    if not cds_list:
        return []
    return [{"chrom": str(c[0]), "start": c[1], "end": c[2]} for c in cds_list]

CFG = yaml.safe_load((ROOT / "configs/exp1_tools.yaml").read_text())
METRICS_DIR = RESULTS_ROOT / "exp1/metrics"


def _empty_eval(n_gt_exons: int):
    return {
        "correct_chrom": False, "correct_locus": False,
        "exon_f1": 0.0, "exon_precision": 0.0, "exon_recall": 0.0,
        "junction_f1": 0.0, "junction_precision": 0.0, "junction_recall": 0.0,
        "cds_recovery": 0.0, "query_coverage": 0.0, "full_length": False,
        "n_pred_exons": 0, "n_gt_exons": n_gt_exons, "exon_count_diff": -n_gt_exons,
        "no_hit": True,
    }


def _species_key(species, suffix):
    return f"{species}{suffix}"


def compute(tool, species, regime):
    suffix = CFG["regimes"][regime]
    sp_key = _species_key(species, suffix)

    raw = json.loads(out_path(tool, species, suffix).read_text())
    gt = load_ground_truth(sp_key)
    prots = load_proteins(sp_key)

    is_paralog = regime == "paralog"
    gene_cds = load_cds_bed(species) if is_paralog else {}
    paralog_info = {}
    if is_paralog:
        pi_path = ROOT / f"data/v3/{sp_key}/paralog_info.json"
        if pi_path.exists():
            paralog_info = json.loads(pi_path.read_text())

    per_protein = {}
    for pid in prots:
        gt_info = gt.get(pid)
        if gt_info is None:
            continue
        plen = len(prots[pid])
        n_gt = len(gt_info.get("exons", []))
        hits = raw["results"].get(pid, [])
        top1 = hits[0] if hits else None
        ev = evaluate_hit(top1, gt_info, plen) if top1 else _empty_eval(n_gt)

        if is_paralog:
            alt_loci = []
            pi = paralog_info.get(pid, {})
            for t in pi.get("all_targets", []):
                loc = get_gene_locus(gene_cds, t.get("target_gene", ""))
                if loc:
                    loc2 = dict(loc)
                    loc2["target_gene"] = t.get("target_gene")
                    loc2["identity"] = t.get("identity", 0)
                    alt_loci.append(loc2)
            if top1:
                e1 = classify_e1(top1, gt_info, alt_loci)
            else:
                e1 = {"strict_e1": True, "broad_e1": False, "original_e1": True,
                      "overlaps_best": False, "overlaps_any_ortholog": False,
                      "no_hit": True, "matched_ortholog": None, "matched_identity": None}
            ev.update({f"e1_{k}": v for k, v in e1.items()})

            # correct_locus_any: overlaps ANY Compara ortholog (strict E1 for A-stage)
            ev["correct_locus_any"] = (not e1.get("strict_e1", True)) if top1 else False

            # best_overlapping_f1: exon F1 against matched ortholog's exons
            if top1:
                matched_gene = e1.get("matched_ortholog")
                best_gene = gt_info.get("target_gene")
                if matched_gene and matched_gene != best_gene:
                    co_exons = get_gene_exons(gene_cds, matched_gene)
                    if co_exons:
                        co_gt = {"chrom": co_exons[0]["chrom"],
                                 "exons": co_exons,
                                 "strand": gt_info.get("strand", "+")}
                        co_ev = evaluate_hit(top1, co_gt, plen)
                        ev["best_overlapping_f1"] = co_ev["exon_f1"]
                    else:
                        ev["best_overlapping_f1"] = ev["exon_f1"]
                else:
                    ev["best_overlapping_f1"] = ev["exon_f1"]
            else:
                ev["best_overlapping_f1"] = 0.0
        else:
            # 1:1 queries: no co-orthologs, these fields equal originals
            ev["correct_locus_any"] = ev["correct_locus"] if top1 else False
            ev["best_overlapping_f1"] = ev["exon_f1"] if top1 else 0.0

        per_protein[pid] = ev

    n = len(per_protein)
    def _mean(key):
        if n == 0:
            return 0.0
        return round(sum(float(v.get(key, 0) or 0) for v in per_protein.values()) / n, 4)

    summary = {
        "tool": tool, "species": species, "regime": regime, "n_proteins": n,
        "hit_rate": round(sum(1 for v in per_protein.values() if not v.get("no_hit")) / max(n, 1), 4),
        "correct_locus_rate": _mean("correct_locus"),
        "correct_locus_any_rate": _mean("correct_locus_any"),
        "exon_f1_mean": _mean("exon_f1"),
        "best_overlapping_f1_mean": _mean("best_overlapping_f1"),
        "junction_f1_mean": _mean("junction_f1"),
        "cds_recovery_mean": _mean("cds_recovery"),
        "full_length_rate": _mean("full_length"),
        "exon_f1_ge_tau": round(sum(1 for v in per_protein.values()
                                    if v.get("exon_f1", 0) >= CFG["evaluation"]["exon_f1_tau"]) / max(n, 1), 4),
    }
    if is_paralog:
        summary["strict_e1_rate"] = round(
            sum(1 for v in per_protein.values() if v.get("e1_strict_e1")) / max(n, 1), 4
        )
        summary["broad_e1_rate"] = round(
            sum(1 for v in per_protein.values() if v.get("e1_broad_e1")) / max(n, 1), 4
        )

    METRICS_DIR.mkdir(parents=True, exist_ok=True)
    stem = f"{tool}.{species}.{regime}"
    (METRICS_DIR / f"{stem}.per_protein.json").write_text(json.dumps(per_protein, indent=1))
    (METRICS_DIR / f"{stem}.summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--regime", required=True, choices=["one_to_one", "paralog"])
    args = ap.parse_args()
    print(json.dumps(compute(args.tool, args.species, args.regime), indent=2))


if __name__ == "__main__":
    main()
