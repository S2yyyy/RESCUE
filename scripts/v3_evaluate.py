#!/usr/bin/env python3
"""V3 Evaluation: compute metrics against ground truth for any tool's output.

Metrics:
  - Exon F1 (base-level overlap)
  - Junction F1 (±10bp tolerance)
  - CDS recovery (fraction of GT CDS bases covered)
  - Full-length rate (≥90% query coverage)
  - Exon count accuracy

Usage:
    python scripts/v3_evaluate.py --tool miniprot --species mouse
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from v3_tool_common import V3_DATA, V3_RESULTS, load_ground_truth, load_proteins


JUNCTION_TOLERANCE = 10  # bp


def intervals_overlap(a_start, a_end, b_start, b_end):
    """Compute overlap between two intervals."""
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def base_level_metrics(pred_cds, gt_cds):
    """Compute base-level precision, recall, F1 between predicted and GT CDS.

    Both are lists of dicts with {chrom, start, end}.
    """
    if not pred_cds or not gt_cds:
        return {"precision": 0, "recall": 0, "f1": 0}

    # Convert to sets of bases (for small exon counts this is fine)
    # For efficiency, use interval overlap
    pred_bases = 0
    gt_bases = 0
    overlap_bases = 0

    for g in gt_cds:
        gt_bases += g["end"] - g["start"]

    for p in pred_cds:
        pred_bases += p["end"] - p["start"]

    # Compute overlap: for each pred exon, sum overlap with each gt exon
    for p in pred_cds:
        for g in gt_cds:
            if p["chrom"] == g["chrom"]:
                overlap_bases += intervals_overlap(p["start"], p["end"],
                                                   g["start"], g["end"])

    precision = overlap_bases / pred_bases if pred_bases > 0 else 0
    recall = overlap_bases / gt_bases if gt_bases > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    return {"precision": round(precision, 4), "recall": round(recall, 4),
            "f1": round(f1, 4)}


def junction_metrics(pred_cds, gt_cds, tolerance=JUNCTION_TOLERANCE):
    """Compute junction-level precision, recall, F1.

    A junction is defined by the end of one exon and start of the next.
    Tolerance: ±10bp.
    """
    def get_junctions(cds_list):
        """Extract internal junctions from sorted CDS list."""
        sorted_cds = sorted(cds_list, key=lambda c: c["start"])
        junctions = []
        for i in range(len(sorted_cds) - 1):
            # Junction = (donor_end, acceptor_start) on same chrom
            if sorted_cds[i]["chrom"] == sorted_cds[i+1]["chrom"]:
                junctions.append((sorted_cds[i]["chrom"],
                                  sorted_cds[i]["end"],
                                  sorted_cds[i+1]["start"]))
        return junctions

    pred_junctions = get_junctions(pred_cds) if pred_cds else []
    gt_junctions = get_junctions(gt_cds) if gt_cds else []

    if not gt_junctions and not pred_junctions:
        return {"precision": 1.0, "recall": 1.0, "f1": 1.0, "n_gt": 0, "n_pred": 0}
    if not gt_junctions:
        return {"precision": 0, "recall": 1.0, "f1": 0, "n_gt": 0, "n_pred": len(pred_junctions)}
    if not pred_junctions:
        return {"precision": 1.0, "recall": 0, "f1": 0, "n_gt": len(gt_junctions), "n_pred": 0}

    # Match pred to gt junctions within tolerance
    matched_gt = set()
    matched_pred = set()

    for pi, pj in enumerate(pred_junctions):
        for gi, gj in enumerate(gt_junctions):
            if gi in matched_gt:
                continue
            if pj[0] == gj[0]:  # Same chrom
                if (abs(pj[1] - gj[1]) <= tolerance and
                    abs(pj[2] - gj[2]) <= tolerance):
                    matched_pred.add(pi)
                    matched_gt.add(gi)
                    break

    precision = len(matched_pred) / len(pred_junctions) if pred_junctions else 0
    recall = len(matched_gt) / len(gt_junctions) if gt_junctions else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    return {"precision": round(precision, 4), "recall": round(recall, 4),
            "f1": round(f1, 4),
            "n_gt": len(gt_junctions), "n_pred": len(pred_junctions)}


def evaluate_hit(hit, gt, protein_len):
    """Evaluate a single hit against ground truth."""
    pred_cds = hit.get("cds", [])
    gt_cds = gt.get("exons", [])  # Ground truth exons

    # Check if hit is on correct locus (same chrom)
    correct_chrom = False
    if pred_cds and gt_cds:
        pred_chrom = pred_cds[0]["chrom"]
        gt_chrom = gt_cds[0]["chrom"]
        correct_chrom = (str(pred_chrom) == str(gt_chrom))

    # Check if on correct locus (overlap with GT region)
    correct_locus = False
    if correct_chrom:
        pred_start = min(c["start"] for c in pred_cds)
        pred_end = max(c["end"] for c in pred_cds)
        gt_start = min(c["start"] for c in gt_cds)
        gt_end = max(c["end"] for c in gt_cds)
        # Overlap: any overlap means same locus
        if intervals_overlap(pred_start, pred_end, gt_start, gt_end) > 0:
            correct_locus = True

    base = base_level_metrics(pred_cds, gt_cds) if correct_chrom else \
           {"precision": 0, "recall": 0, "f1": 0}
    junc = junction_metrics(pred_cds, gt_cds) if correct_chrom else \
           {"precision": 0, "recall": 0, "f1": 0, "n_gt": len(gt_cds)-1, "n_pred": 0}

    # CDS recovery
    cds_recovery = base["recall"]

    # Query coverage — use query coordinates if available, else estimate from CDS bases
    query_coverage = 0
    if protein_len > 0 and hit.get("query_end", 0) > 0:
        covered = hit["query_end"] - hit.get("query_start", 1) + 1
        query_coverage = covered / protein_len
    elif protein_len > 0 and pred_cds:
        # Estimate from predicted CDS nucleotide length / (3 * protein_len)
        pred_nt = sum(c["end"] - c["start"] for c in pred_cds)
        query_coverage = min(1.0, pred_nt / (3 * protein_len))

    # Exon count accuracy
    n_pred_exons = len(pred_cds)
    n_gt_exons = len(gt_cds)
    exon_count_diff = n_pred_exons - n_gt_exons

    return {
        "correct_chrom": correct_chrom,
        "correct_locus": correct_locus,
        "exon_f1": base["f1"],
        "exon_precision": base["precision"],
        "exon_recall": base["recall"],
        "junction_f1": junc["f1"],
        "junction_precision": junc["precision"],
        "junction_recall": junc["recall"],
        "cds_recovery": cds_recovery,
        "query_coverage": query_coverage,
        "full_length": query_coverage >= 0.9,
        "n_pred_exons": n_pred_exons,
        "n_gt_exons": n_gt_exons,
        "exon_count_diff": exon_count_diff,
    }


def evaluate_tool(tool, species):
    """Evaluate all results for a tool × species combination."""
    result_path = V3_RESULTS / "raw" / f"{tool}_{species}.json"
    if not result_path.exists():
        print(f"  No results for {tool}/{species}")
        return None

    with open(result_path) as fh:
        data = json.load(fh)

    gt = load_ground_truth(species)
    proteins = load_proteins(species)

    tool_results = data.get("results", {})
    per_protein = {}

    for pid, hits in tool_results.items():
        if pid not in gt:
            continue

        protein_len = len(proteins.get(pid, ""))
        gt_info = gt[pid]

        if not hits:
            per_protein[pid] = {
                "n_hits": 0,
                "top1": None,
                "best": None,
            }
            continue

        # Evaluate top-1
        top1_eval = evaluate_hit(hits[0], gt_info, protein_len)

        # Evaluate all hits, find best by F1
        all_evals = [evaluate_hit(h, gt_info, protein_len) for h in hits]
        best_idx = max(range(len(all_evals)), key=lambda i: all_evals[i]["exon_f1"])
        best_eval = all_evals[best_idx]

        per_protein[pid] = {
            "n_hits": len(hits),
            "top1": top1_eval,
            "best": best_eval,
            "best_rank": best_idx + 1,
            "all_evals": all_evals,
        }

    return per_protein


def summarize(per_protein, strata=None):
    """Compute summary statistics from per-protein evaluation."""
    metrics = defaultdict(list)

    for pid, result in per_protein.items():
        if result["top1"] is None:
            metrics["exon_f1"].append(0)
            metrics["junction_f1"].append(0)
            metrics["cds_recovery"].append(0)
            metrics["correct_locus"].append(0)
            metrics["full_length"].append(0)
            metrics["no_hit"].append(1)
            continue

        top1 = result["top1"]
        metrics["exon_f1"].append(top1["exon_f1"])
        metrics["junction_f1"].append(top1["junction_f1"])
        metrics["cds_recovery"].append(top1["cds_recovery"])
        metrics["correct_locus"].append(1 if top1["correct_locus"] else 0)
        metrics["full_length"].append(1 if top1["full_length"] else 0)
        metrics["no_hit"].append(0)

        # Best hit (oracle over ranks)
        best = result["best"]
        metrics["oracle_exon_f1"].append(best["exon_f1"])
        metrics["oracle_correct_locus"].append(1 if best["correct_locus"] else 0)

    summary = {}
    for key, values in metrics.items():
        arr = np.array(values)
        summary[key] = {
            "mean": round(float(arr.mean()), 4),
            "std": round(float(arr.std()), 4),
            "n": len(arr),
        }

    return summary


def main():
    parser = argparse.ArgumentParser(description="Evaluate v3 tool results")
    parser.add_argument("--tool", required=True,
                        choices=["miniprot", "spaln", "gemoma", "exonerate"])
    parser.add_argument("--species", required=True,
                        choices=["mouse", "dog", "chicken", "zebrafish"])
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()

    print(f"Evaluating {args.tool} on {args.species}")

    per_protein = evaluate_tool(args.tool, args.species)
    if per_protein is None:
        sys.exit(1)

    summary = summarize(per_protein)

    # Load strata for stratified evaluation
    strata_path = V3_DATA / args.species / "strata.json"
    strata = json.load(open(strata_path)) if strata_path.exists() else {}

    strata_summaries = {}
    if strata:
        strata_groups = defaultdict(dict)
        for pid, labels in strata.items():
            for label in labels:
                if pid in per_protein:
                    strata_groups[label][pid] = per_protein[pid]

        for label, group in strata_groups.items():
            strata_summaries[label] = summarize(group)

    result = {
        "tool": args.tool,
        "species": args.species,
        "n_evaluated": len(per_protein),
        "summary": summary,
        "strata": strata_summaries,
        "per_protein": per_protein,
    }

    print(f"\n  Summary ({args.tool} / {args.species}):")
    for key in ["exon_f1", "junction_f1", "correct_locus", "cds_recovery", "no_hit"]:
        if key in summary:
            print(f"    {key}: {summary[key]['mean']:.4f} (±{summary[key]['std']:.4f})")

    if strata_summaries:
        print(f"\n  By stratum:")
        for label, ss in sorted(strata_summaries.items()):
            if "exon_f1" in ss:
                print(f"    {label} (n={ss['exon_f1']['n']}): "
                      f"F1={ss['exon_f1']['mean']:.3f}, "
                      f"locus={ss['correct_locus']['mean']:.3f}")

    out_path = args.output or str(V3_RESULTS / f"eval_{args.tool}_{args.species}.json")
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"\n  Saved to {out_path}")


if __name__ == "__main__":
    main()
