#!/usr/bin/env python3
"""V3 Paralog E1 scoring audit: distinguish strict vs broad vs original E1.

Problem: the original E1 (wrong locus) counts a hit as E1 if it doesn't overlap
the single best Compara ortholog. But for one2many / many2many orthologs, the
protein may legitimately map to ANY of its Compara orthologs.

This script re-evaluates E1 with three definitions:
  - Original E1: hit doesn't overlap the primary GT locus (current definition)
  - Broad E1:    hit overlaps a co-ortholog but NOT the best one
  - Strict E1:   hit overlaps NONE of the valid Compara orthologs (true misloc)

Usage:
    python scripts/v3_fix_paralog_scoring.py
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from v3_tool_common import V3_DATA, V3_RESULTS
from rescue_paths import GENE_STRUCTURES

SPECIES = ["mouse", "dog", "chicken", "zebrafish"]
TOOLS = ["miniprot", "exonerate"]
CDS_BED_DIR = GENE_STRUCTURES


def intervals_overlap(a_start, a_end, b_start, b_end):
    """Compute overlap between two intervals (1-based inclusive)."""
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def load_cds_bed(species):
    """Load CDS BED file into a dict: gene_id -> list of (chrom, start, end).

    BED format: chrom start end gene_id transcript_id strand phase
    BED coords are 0-based half-open; we convert to 1-based inclusive.
    """
    bed_path = CDS_BED_DIR / f"{species}_CDS.bed"
    if not bed_path.exists():
        print(f"  WARNING: CDS BED not found: {bed_path}")
        return {}

    gene_cds = defaultdict(list)
    with open(bed_path) as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start = int(parts[1]) + 1  # BED 0-based -> 1-based
            end = int(parts[2])        # BED half-open end = inclusive end
            gene_id = parts[3]
            gene_cds[gene_id].append((chrom, start, end))

    return gene_cds


def get_gene_locus(gene_cds_dict, gene_id):
    """Get the genomic locus (chrom, min_start, max_end) for a gene from CDS dict."""
    cds_list = gene_cds_dict.get(gene_id, [])
    if not cds_list:
        return None
    chrom = cds_list[0][0]
    start = min(c[1] for c in cds_list)
    end = max(c[2] for c in cds_list)
    return {"chrom": chrom, "start": start, "end": end}


def hit_overlaps_locus(hit, locus):
    """Check if a tool hit overlaps with a given locus."""
    pred_cds = hit.get("cds", [])
    if not pred_cds or locus is None:
        return False
    pred_chrom = str(pred_cds[0]["chrom"])
    locus_chrom = str(locus["chrom"])
    if pred_chrom != locus_chrom:
        return False
    pred_start = min(c["start"] for c in pred_cds)
    pred_end = max(c["end"] for c in pred_cds)
    return intervals_overlap(pred_start, pred_end, locus["start"], locus["end"]) > 0


def classify_e1(hit, gt_info, alt_loci):
    """Classify a hit's E1 status.

    Args:
        hit: tool alignment result (top-1)
        gt_info: ground truth dict (with exons for best ortholog)
        alt_loci: list of {chrom, start, end, target_gene, identity} for ALL
                  Compara orthologs (including the best one)

    Returns:
        dict with:
          original_e1: True if hit doesn't overlap primary GT locus
          overlaps_best: True if hit overlaps the best ortholog
          overlaps_any_ortholog: True if hit overlaps any Compara ortholog
          matched_ortholog: gene_id of the matched ortholog (if any)
          matched_identity: identity of the matched ortholog
          strict_e1: True if hit overlaps NONE of the orthologs
          broad_e1: True if hit overlaps an alt ortholog but not the best
    """
    pred_cds = hit.get("cds", [])
    if not pred_cds:
        return {
            "original_e1": False, "overlaps_best": False,
            "overlaps_any_ortholog": False, "matched_ortholog": None,
            "matched_identity": None, "strict_e1": False, "broad_e1": False,
            "no_hit": True,
        }

    # Check overlap with primary GT locus
    gt_cds = gt_info.get("exons", [])
    overlaps_best = False
    if gt_cds:
        gt_locus = {
            "chrom": gt_cds[0]["chrom"],
            "start": min(c["start"] for c in gt_cds),
            "end": max(c["end"] for c in gt_cds),
        }
        overlaps_best = hit_overlaps_locus(hit, gt_locus)

    original_e1 = not overlaps_best

    # Check overlap with any ortholog locus
    overlaps_any = False
    matched_gene = None
    matched_ident = None
    for locus in alt_loci:
        if hit_overlaps_locus(hit, locus):
            overlaps_any = True
            # Track the highest-identity match
            if matched_ident is None or locus.get("identity", 0) > matched_ident:
                matched_gene = locus.get("target_gene")
                matched_ident = locus.get("identity", 0)

    strict_e1 = not overlaps_any  # doesn't overlap ANY ortholog
    broad_e1 = overlaps_any and not overlaps_best  # overlaps alt but not best

    return {
        "original_e1": original_e1,
        "overlaps_best": overlaps_best,
        "overlaps_any_ortholog": overlaps_any,
        "matched_ortholog": matched_gene,
        "matched_identity": matched_ident,
        "strict_e1": strict_e1,
        "broad_e1": broad_e1,
        "no_hit": False,
    }


def get_orthology_type(n_targets):
    """Classify as one2many vs many2many based on target count."""
    # We don't have explicit Ensembl homology type, but:
    # one2many: human gene -> multiple target genes
    # many2many: multiple human genes -> multiple target genes
    # For simplicity, we classify by n_targets
    if n_targets <= 1:
        return "one2one"
    elif n_targets <= 5:
        return "one2few"
    else:
        return "one2many_or_m2m"


def run_audit(tool, species, gene_cds_dict):
    """Run the E1 audit for one tool/species combination."""
    # Load data
    paralog_info_path = V3_DATA / f"{species}_paralog" / "paralog_info.json"
    gt_path = V3_DATA / f"{species}_paralog" / "ground_truth.json"
    raw_path = V3_RESULTS / "raw" / f"{tool}_{species}_paralog.json"

    for p in [paralog_info_path, gt_path, raw_path]:
        if not p.exists():
            print(f"  SKIP {tool}/{species}: missing {p.name}")
            return None

    with open(paralog_info_path) as f:
        paralog_info = json.load(f)
    with open(gt_path) as f:
        gt = json.load(f)
    with open(raw_path) as f:
        raw = json.load(f)

    tool_results = raw.get("results", {})

    # Build locus info for all ortholog targets per protein
    per_protein = {}
    for pid in gt:
        if pid not in tool_results:
            continue

        gt_info = gt[pid]
        pi = paralog_info.get(pid, {})
        all_targets = pi.get("all_targets", [])

        # Build loci for all orthologs from CDS BED
        alt_loci = []
        for tgt in all_targets:
            gene_id = tgt["target_gene"]
            locus = get_gene_locus(gene_cds_dict, gene_id)
            if locus is not None:
                alt_loci.append({
                    **locus,
                    "target_gene": gene_id,
                    "identity": tgt.get("identity", 0),
                    "is_best": tgt.get("is_best", False),
                })

        hits = tool_results[pid]
        if not hits:
            per_protein[pid] = {
                "no_hit": True,
                "n_orthologs": len(all_targets),
                "n_orthologs_with_cds": len(alt_loci),
            }
            continue

        top1 = hits[0]
        classification = classify_e1(top1, gt_info, alt_loci)

        per_protein[pid] = {
            **classification,
            "n_orthologs": len(all_targets),
            "n_orthologs_with_cds": len(alt_loci),
            "hit_chrom": str(top1["cds"][0]["chrom"]) if top1.get("cds") else None,
            "hit_start": min(c["start"] for c in top1["cds"]) if top1.get("cds") else None,
            "hit_end": max(c["end"] for c in top1["cds"]) if top1.get("cds") else None,
            "gt_gene": gt_info.get("target_gene"),
            "gt_chrom": gt_info.get("chrom"),
        }

    return per_protein


def summarize_results(per_protein):
    """Compute summary rates from per-protein classifications."""
    counters = defaultdict(int)
    # Also break down by orthology type
    by_type = defaultdict(lambda: defaultdict(int))

    for pid, r in per_protein.items():
        n_ort = r.get("n_orthologs", 0)
        otype = get_orthology_type(n_ort)

        if r.get("no_hit"):
            counters["no_hit"] += 1
            by_type[otype]["no_hit"] += 1
            by_type[otype]["total"] += 1
            counters["total"] += 1
            continue

        counters["total"] += 1
        counters["has_hit"] += 1
        by_type[otype]["total"] += 1
        by_type[otype]["has_hit"] += 1

        if r.get("overlaps_best"):
            counters["overlaps_best"] += 1
            by_type[otype]["overlaps_best"] += 1

        if r.get("overlaps_any_ortholog"):
            counters["overlaps_any"] += 1
            by_type[otype]["overlaps_any"] += 1

        if r.get("original_e1"):
            counters["original_e1"] += 1
            by_type[otype]["original_e1"] += 1

        if r.get("broad_e1"):
            counters["broad_e1"] += 1
            by_type[otype]["broad_e1"] += 1

        if r.get("strict_e1"):
            counters["strict_e1"] += 1
            by_type[otype]["strict_e1"] += 1

    total = counters["total"]
    has_hit = counters["has_hit"]

    summary = {
        "total": total,
        "has_hit": has_hit,
        "no_hit": counters["no_hit"],
        "overlaps_best": counters["overlaps_best"],
        "overlaps_any": counters["overlaps_any"],
        "original_e1_count": counters["original_e1"],
        "broad_e1_count": counters["broad_e1"],
        "strict_e1_count": counters["strict_e1"],
    }

    # Rates (denominator = proteins with hits)
    if has_hit > 0:
        summary["original_e1_rate"] = round(counters["original_e1"] / has_hit, 4)
        summary["broad_e1_rate"] = round(counters["broad_e1"] / has_hit, 4)
        summary["strict_e1_rate"] = round(counters["strict_e1"] / has_hit, 4)
        summary["reclassified_count"] = counters["original_e1"] - counters["strict_e1"]
        summary["reclassified_frac"] = round(
            (counters["original_e1"] - counters["strict_e1"]) / max(1, counters["original_e1"]), 4
        )
    else:
        summary["original_e1_rate"] = 0
        summary["broad_e1_rate"] = 0
        summary["strict_e1_rate"] = 0
        summary["reclassified_count"] = 0
        summary["reclassified_frac"] = 0

    summary["by_orthology_type"] = {}
    for otype, counts in sorted(by_type.items()):
        ht = counts.get("has_hit", 0)
        summary["by_orthology_type"][otype] = {
            "total": counts["total"],
            "has_hit": ht,
            "original_e1": counts.get("original_e1", 0),
            "broad_e1": counts.get("broad_e1", 0),
            "strict_e1": counts.get("strict_e1", 0),
            "original_e1_rate": round(counts.get("original_e1", 0) / max(1, ht), 4) if ht else 0,
            "strict_e1_rate": round(counts.get("strict_e1", 0) / max(1, ht), 4) if ht else 0,
        }

    return summary


def build_audit_cases(per_protein, n=20):
    """Build top-N E1 audit cases sorted by original_e1 proteins."""
    cases = []
    for pid, r in per_protein.items():
        if r.get("no_hit"):
            continue
        if not r.get("original_e1"):
            continue
        cases.append({
            "protein_id": pid,
            "hit_location": f"{r.get('hit_chrom')}:{r.get('hit_start')}-{r.get('hit_end')}",
            "gt_gene": r.get("gt_gene"),
            "gt_chrom": r.get("gt_chrom"),
            "n_orthologs": r.get("n_orthologs"),
            "n_orthologs_with_cds": r.get("n_orthologs_with_cds"),
            "overlaps_any_ortholog": r.get("overlaps_any_ortholog"),
            "matched_ortholog": r.get("matched_ortholog"),
            "matched_identity": r.get("matched_identity"),
            "strict_e1": r.get("strict_e1"),
            "broad_e1": r.get("broad_e1"),
        })

    # Sort: strict E1 first (true mislocalization), then by n_orthologs descending
    cases.sort(key=lambda x: (not x["strict_e1"], -x["n_orthologs"]))
    return cases[:n]


def main():
    print("=" * 70)
    print("V3 PARALOG E1 SCORING AUDIT")
    print("Re-evaluating E1 with strict/broad/original definitions")
    print("=" * 70)

    all_results = {}

    for species in SPECIES:
        print(f"\n--- Loading CDS BED for {species} ---")
        gene_cds = load_cds_bed(species)
        print(f"  Loaded {len(gene_cds)} genes from CDS BED")

        for tool in TOOLS:
            key = f"{tool}/{species}"
            print(f"\n  Evaluating {key}...")

            per_protein = run_audit(tool, species, gene_cds)
            if per_protein is None:
                continue

            summary = summarize_results(per_protein)
            audit_cases = build_audit_cases(per_protein, n=20)

            all_results[key] = {
                "tool": tool,
                "species": species,
                "summary": summary,
                "audit_cases": audit_cases,
                "per_protein": per_protein,
            }

            # Print summary
            s = summary
            print(f"    Total proteins: {s['total']}, with hits: {s['has_hit']}")
            print(f"    Original E1: {s['original_e1_count']}/{s['has_hit']} "
                  f"= {s['original_e1_rate']:.4f}")
            print(f"    Broad E1 (alt ortholog hit, not best): {s['broad_e1_count']}/{s['has_hit']} "
                  f"= {s['broad_e1_rate']:.4f}")
            print(f"    Strict E1 (no ortholog overlap):       {s['strict_e1_count']}/{s['has_hit']} "
                  f"= {s['strict_e1_rate']:.4f}")
            print(f"    Reclassified: {s['reclassified_count']} "
                  f"({s['reclassified_frac']:.1%} of original E1)")

            # By orthology type
            print(f"    By orthology type:")
            for otype, tc in sorted(s["by_orthology_type"].items()):
                print(f"      {otype:>20s}: n={tc['total']:>4d}  "
                      f"orig_e1={tc['original_e1_rate']:.4f}  "
                      f"strict_e1={tc['strict_e1_rate']:.4f}")

    # Print comparison table
    print("\n" + "=" * 90)
    print("COMPARISON TABLE: Original E1 vs Strict E1")
    print("=" * 90)
    header = f"  {'Tool':<12} {'Species':<12} {'N':>5} {'Hits':>5} "
    header += f"{'Orig E1':>10} {'Broad E1':>10} {'Strict E1':>10} {'Reclass%':>10}"
    print(header)
    print("  " + "-" * 85)

    for species in SPECIES:
        for tool in TOOLS:
            key = f"{tool}/{species}"
            if key not in all_results:
                continue
            s = all_results[key]["summary"]
            print(f"  {tool:<12} {species:<12} {s['total']:>5} {s['has_hit']:>5} "
                  f"{s['original_e1_rate']:>10.4f} {s['broad_e1_rate']:>10.4f} "
                  f"{s['strict_e1_rate']:>10.4f} {s['reclassified_frac']:>9.1%}")
        print()

    # Print sample audit cases
    print("\n" + "=" * 90)
    print("SAMPLE AUDIT CASES (first tool/species with E1)")
    print("=" * 90)
    for key, data in all_results.items():
        cases = data["audit_cases"]
        if not cases:
            continue
        print(f"\n  --- {key} (top {len(cases)} original-E1 cases) ---")
        for i, c in enumerate(cases[:10]):
            tag = "STRICT" if c["strict_e1"] else "BROAD"
            match_info = ""
            if c["matched_ortholog"]:
                match_info = f" -> matched {c['matched_ortholog']} (ident={c['matched_identity']:.1f}%)"
            print(f"    [{tag:>6}] {c['protein_id']}  hit={c['hit_location']}  "
                  f"gt={c['gt_gene']}({c['gt_chrom']})  "
                  f"n_ort={c['n_orthologs']}/{c['n_orthologs_with_cds']}"
                  f"{match_info}")

    # Save results (strip per_protein for file size)
    output = {
        "description": "Paralog E1 scoring audit: strict vs broad vs original E1",
        "e1_definitions": {
            "original_e1": "hit does not overlap primary GT locus (current definition)",
            "broad_e1": "hit overlaps a co-ortholog but NOT the best one",
            "strict_e1": "hit overlaps NONE of the valid Compara orthologs (true mislocalization)",
        },
        "results": {},
    }
    for key, data in all_results.items():
        output["results"][key] = {
            "tool": data["tool"],
            "species": data["species"],
            "summary": data["summary"],
            "audit_cases": data["audit_cases"],
        }

    out_path = V3_RESULTS / "paralog_e1_audit.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(output, fh, indent=2)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
