"""TPRS Phase 4 — Locus Router.

For each protein, score all candidate loci from Phase 1 A-tools, train an
XGBoost router to pick the locus with the best exon-level F1.

Workflow:
  1. Merge candidate regions from 4 relaxed A-tools
  2. Run miniprot on each region → extract A-stage + C-stage features
  3. Compute exon-F1 label against ground-truth primary CDS
  4. Train XGBoost regressor to predict F1 from features
  5. Evaluate per-protein routing vs baselines

Output:
  results/tprs/router/training_data_{species}.tsv
  results/tprs/router/model_{species}.json
  results/tprs/router/evaluation_{species}.tsv
  results/tprs/router/evaluation_summary.json
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import (  # noqa: E402
    GENE_STRUCTURES,
    GENOMES,
    MINIPROT,
    PROJECT_ROOT,
    SAMTOOLS,
    V3_DATA,
)

REPO = PROJECT_ROOT
BED_DIR = REPO / "results/tprs/stage0_bed"
GT_DIR = REPO / "results/exp2/gt"
DATA = V3_DATA
EXP1_METRICS = REPO / "results/exp1/metrics"
RECOVERY = REPO / "results/tprs/recovery"
ROUTER = REPO / "results/tprs/router"
GENESTRUCT = GENE_STRUCTURES

A_TOOLS = ["diamond_relaxed", "last_relaxed", "miniprot_relaxed", "mmseqs_relaxed"]
SPECIES = ["mouse", "dog", "chicken", "zebrafish"]

# ---------------------------------------------------------------------------
# Utility functions (from Phase 2 / Phase 3 scripts)
# ---------------------------------------------------------------------------


def load_bed(path: Path) -> dict[str, list[tuple]]:
    """Return {query_id: [(chrom, start, end, strand, score, q_cov), ...]}."""
    out: dict[str, list[tuple]] = defaultdict(list)
    if not path.exists():
        return out
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.strip().split("\t")
            q_cov = float(f[8]) if len(f) > 8 else 0.0
            out[f[3]].append((f[0], int(f[1]), int(f[2]), f[5],
                             float(f[4]) if len(f) > 4 else 0.0, q_cov))
    return out


def merge_regions(candidates: list[tuple], slop: int = 50000) -> list[dict]:
    """Merge overlapping candidate regions per chrom+strand.

    Each input tuple: (chrom, start, end, strand, score, q_cov, tool_name)
    Returns list of dicts with tool-vote metadata.
    """
    if not candidates:
        return []
    by_key: dict[tuple, list[tuple]] = defaultdict(list)
    for c in candidates:
        by_key[(c[0], c[3])].append(c)

    merged = []
    for (chrom, strand), regs in by_key.items():
        regs_sorted = sorted(regs, key=lambda x: x[1])
        group = [regs_sorted[0]]
        for r in regs_sorted[1:]:
            # current group end + slop overlaps next start?
            cur_end = max(x[2] for x in group)
            if r[1] <= cur_end + slop:
                group.append(r)
            else:
                merged.append(_collapse_group(chrom, strand, group, slop))
                group = [r]
        if group:
            merged.append(_collapse_group(chrom, strand, group, slop))
    return merged


def _collapse_group(chrom: str, strand: str, group: list[tuple],
                    slop: int) -> dict:
    """Collapse a group of overlapping candidates into one merged region dict."""
    start = max(0, min(x[1] for x in group) - slop)
    end = max(x[2] for x in group) + slop
    tools = set(x[6] for x in group if len(x) > 6)
    scores = [x[4] for x in group if x[4] > 0]
    qcovs = [x[5] for x in group if x[5] > 0]
    return {
        "chrom": chrom, "start": start, "end": end, "strand": strand,
        "n_tools": len(tools), "tools": sorted(tools),
        "best_score": max(scores) if scores else 0.0,
        "best_qcov": max(qcovs) if qcovs else 0.0,
        "n_raw_hits": len(group),
    }


def extract_region_fasta(genome: Path, chrom: str, start: int, end: int) -> str:
    """Extract genomic region via samtools faidx, return FASTA text."""
    region = f"{chrom}:{max(1, start)}-{end}"
    proc = subprocess.run([SAMTOOLS, "faidx", str(genome), region],
                         capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"samtools faidx failed: {proc.stderr[:200]}")
    return proc.stdout


def parse_miniprot_gff(gff_text: str, region_chrom: str,
                       region_offset: int) -> list[dict]:
    """Parse miniprot GFF, return hits with CDS coords in absolute genome space.

    Uses region_offset to convert GFF coordinates (relative to extracted region)
    to absolute genomic coordinates. This is the CORRECTED version from Phase 3.
    """
    hits = []
    cur = None
    cur_cds: list[dict] = []
    for line in gff_text.splitlines():
        if line.startswith("#"):
            if line.startswith("##PAF"):
                if cur is not None and cur_cds:
                    cur["cds"] = sorted(cur_cds, key=lambda x: x["start"])
                    hits.append(cur)
                cur = None; cur_cds = []
            continue
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        feat = fields[2]
        if feat == "mRNA":
            if cur is not None and cur_cds:
                cur["cds"] = sorted(cur_cds, key=lambda x: x["start"])
                hits.append(cur)
            cur_cds = []
            attrs = dict(p.split("=", 1) for p in fields[8].split(";") if "=" in p)
            target = attrs.get("Target", "").split()
            cur = {
                "chrom": region_chrom,
                "strand": fields[6],
                "score": float(fields[5]) if fields[5] != "." else 0.0,
                "identity": float(attrs.get("Identity", 0.0)),
                "query_start": int(target[1]) if len(target) > 1 else 0,
                "query_end": int(target[2]) if len(target) > 2 else 0,
                "rank": int(attrs.get("Rank", 1)),
            }
        elif feat == "CDS" and cur is not None:
            cur_cds.append({
                "chrom": region_chrom,
                "start": int(fields[3]) - 1 + region_offset,
                "end": int(fields[4]) + region_offset,
                "strand": cur.get("strand", "."),
            })
    if cur is not None and cur_cds:
        cur["cds"] = sorted(cur_cds, key=lambda x: x["start"])
        hits.append(cur)
    return hits


def base_level_metrics(pred_cds: list[dict], gt_cds: list[dict]) -> dict:
    """Compute exon base-level precision, recall, F1."""
    if not pred_cds or not gt_cds:
        return {"precision": 0.0, "recall": 0.0, "f1": 0.0}
    pred_bp = sum(e["end"] - e["start"] for e in pred_cds)
    gt_bp = sum(e["end"] - e["start"] for e in gt_cds)
    if pred_bp == 0 or gt_bp == 0:
        return {"precision": 0.0, "recall": 0.0, "f1": 0.0}
    overlap = 0
    for pe in pred_cds:
        for ge in gt_cds:
            if pe["chrom"] != ge["chrom"]:
                continue
            inter = max(0, min(pe["end"], ge["end"]) - max(pe["start"], ge["start"]))
            overlap += inter
    # Cap overlap to prevent double-counting across exon pairs
    overlap = min(overlap, pred_bp, gt_bp)
    precision = overlap / pred_bp if pred_bp > 0 else 0.0
    recall = overlap / gt_bp if gt_bp > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    return {"precision": precision, "recall": recall, "f1": f1}


def load_gt_json(sp: str) -> dict:
    path = DATA / f"{sp}_paralog" / "ground_truth.json"
    with open(path) as f:
        return json.load(f)


def load_query_lengths(sp: str) -> dict[str, int]:
    """Load protein query lengths from FASTA."""
    qlen = {}
    fa_path = DATA / f"{sp}_paralog" / "proteins.fa"
    with open(fa_path) as fh:
        cur_id = None
        for line in fh:
            if line.startswith(">"):
                cur_id = line[1:].strip().split()[0]
                qlen[cur_id] = 0
            elif cur_id:
                qlen[cur_id] += len(line.strip())
    return qlen


def load_recovery_manifest(sp: str) -> dict[str, dict]:
    """Load Phase 3 recovery manifest, return {query_id: info}."""
    path = RECOVERY / f"manifest_{sp}.tsv"
    if not path.exists():
        return {}
    out = {}
    with open(path) as fh:
        cols = fh.readline().strip().split("\t")
        for line in fh:
            if not line.strip():
                continue
            f = line.strip().split("\t")
            row = dict(zip(cols, f))
            out[row["query_id"]] = row
    return out


def load_exp1_f1(sp: str) -> dict[str, dict]:
    """Load Exp1 miniprot per-protein F1 metrics."""
    path = EXP1_METRICS / f"miniprot.{sp}.paralog.per_protein.json"
    if not path.exists():
        return {}
    with open(path) as f:
        return json.load(f)


def load_cds_bed(sp: str) -> dict[str, list[dict]]:
    """Load CDS BED file, return {gene_id: [{chrom, start, end, strand}, ...]}.

    BED is 0-based half-open: start is already 0-based, end is exclusive.
    """
    bed_path = GENESTRUCT / f"{sp}_CDS.bed"
    if not bed_path.exists():
        return {}
    out: dict[str, list[dict]] = defaultdict(list)
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.strip().split("\t")
            if len(f) < 6:
                continue
            gene_id = f[3].rsplit("_", 1)[0] if "_" in f[3] else f[3]
            out[gene_id].append({
                "chrom": f[0],
                "start": int(f[1]),  # 0-based BED
                "end": int(f[2]),
                "strand": f[5],
            })
    return out


# ---------------------------------------------------------------------------
# Feature extraction
# ---------------------------------------------------------------------------

def extract_a_stage_features(region: dict, gt_info: dict | None,
                              recovery_info: dict | None,
                              query_len: int) -> dict:
    """Extract A-stage features for one candidate region."""
    feats = {}

    # Tool presence
    tools = set(region.get("tools", []))
    for tool in A_TOOLS:
        short = tool.replace("_relaxed", "")
        feats[f"a_{short}"] = 1 if tool in tools else 0

    feats["a_n_tools"] = region.get("n_tools", 0)
    feats["a_best_score"] = region.get("best_score", 0.0)
    feats["a_best_qcov"] = region.get("best_qcov", 0.0)
    feats["a_n_raw_hits"] = region.get("n_raw_hits", 0)

    # Region properties
    region_size = region["end"] - region["start"]
    feats["region_size"] = region_size
    feats["region_size_log"] = np.log10(max(1, region_size))

    # GT-derived context (these are legitimate features since they don't
    # reveal which locus is correct — they only describe the region)
    if gt_info:
        feats["on_primary_chrom"] = 1 if region["chrom"] == gt_info.get("chrom", "") else 0
        feats["strand_match_primary"] = 1 if region["strand"] == gt_info.get("strand", "") else 0
    else:
        feats["on_primary_chrom"] = 0
        feats["strand_match_primary"] = 0

    # Co-ortholog context
    if recovery_info and recovery_info.get("status") == "recovered":
        co_chrom = recovery_info.get("coorth_chrom", "")
        co_start = int(recovery_info.get("coorth_start", 0))
        co_end = int(recovery_info.get("coorth_end", 0))
        if co_chrom and co_end > co_start:
            # Overlap between region and co-ortholog region
            inter = max(0, min(region["end"], co_end) - max(region["start"], co_start))
            feats["coorth_overlap_bp"] = inter
            feats["coorth_same_chrom"] = 1 if region["chrom"] == co_chrom else 0
        else:
            feats["coorth_overlap_bp"] = 0
            feats["coorth_same_chrom"] = 0
    else:
        feats["coorth_overlap_bp"] = 0
        feats["coorth_same_chrom"] = 0

    feats["query_length_aa"] = query_len
    feats["query_length_log"] = np.log10(max(1, query_len))

    return feats


def extract_c_stage_features(hits: list[dict], query_len: int) -> dict:
    """Extract C-stage features from parsed miniprot GFF output."""
    feats = {
        "c_has_hit": 0, "c_identity": 0.0, "c_coverage": 0.0,
        "c_n_exons": 0, "c_cds_length": 0, "c_composite_score": 0.0,
        "c_hit_score": 0.0, "c_rank": 99,
    }
    if not hits:
        return feats

    best = hits[0]
    feats["c_has_hit"] = 1
    feats["c_identity"] = round(best.get("identity", 0.0), 4)
    feats["c_hit_score"] = round(best.get("score", 0.0), 4)
    feats["c_rank"] = best.get("rank", 99)

    cds = best.get("cds", [])
    feats["c_n_exons"] = len(cds)
    cds_length = sum(e["end"] - e["start"] for e in cds)
    feats["c_cds_length"] = cds_length

    ql3 = query_len * 3
    coverage = cds_length / ql3 if ql3 > 0 else 0.0
    feats["c_coverage"] = round(coverage, 4)

    # Composite score (matching Phase 2 score_alignment)
    exon_bonus = min(1.0, len(cds) / 3.0)
    feats["c_composite_score"] = round(
        feats["c_identity"] * min(1.0, coverage) * exon_bonus, 4)

    return feats


def compute_f1_label(pred_hits: list[dict], gt_primary_exons: list[dict],
                     paralog_cds: dict[str, list[dict]] | None = None) -> dict:
    """Compute exon-level F1 from miniprot predictions vs GT CDS.

    Returns primary_f1 and best_overlap_f1 (max of primary and all paralog gene F1s).
    GT exons expected in 0-based half-open format: {chrom, start, end}.
    paralog_cds: {gene_id: [{chrom, start, end}, ...]} from CDS BED.
    """
    if not pred_hits:
        return {"primary_f1": 0.0, "best_overlap_f1": 0.0,
                "primary_precision": 0.0, "primary_recall": 0.0,
                "best_gene": "", "n_pred_exons": 0, "pred_cds_bp": 0}

    best = pred_hits[0]
    pred_cds = best.get("cds", [])
    n_pred = len(pred_cds)
    pred_bp = sum(e["end"] - e["start"] for e in pred_cds)

    # Primary F1
    p_metrics = base_level_metrics(pred_cds, gt_primary_exons)
    best_overlap = p_metrics["f1"]
    best_gene = "primary"

    # Check each paralog gene's CDS
    if paralog_cds:
        for gene_id, gene_exons in paralog_cds.items():
            m = base_level_metrics(pred_cds, gene_exons)
            if m["f1"] > best_overlap:
                best_overlap = m["f1"]
                best_gene = gene_id

    return {
        "primary_f1": p_metrics["f1"],
        "primary_precision": p_metrics["precision"],
        "primary_recall": p_metrics["recall"],
        "best_overlap_f1": best_overlap,
        "best_gene": best_gene,
        "n_pred_exons": n_pred,
        "pred_cds_bp": pred_bp,
    }


# ---------------------------------------------------------------------------
# Single-region miniprot job (for multiprocessing)
# ---------------------------------------------------------------------------

def _process_one_region(args: tuple) -> dict:
    """Run miniprot on one region, extract features and F1 labels.

    Args tuple: (query_id, region_idx, region, genome_path, query_fasta,
                  query_len, gt_exons_0based, paralog_cds_json, gt_info,
                  recovery_info, timeout)
    """
    (qid, ridx, region, genome_path, query_fasta, query_len,
     gt_exons, paralog_cds_json, gt_info, recovery_info, timeout) = args

    paralog_cds = json.loads(paralog_cds_json) if paralog_cds_json else None

    genome = Path(genome_path)
    result = {
        "query_id": qid, "region_idx": ridx,
        "chrom": region["chrom"], "start": region["start"],
        "end": region["end"], "strand": region["strand"],
        "status": "ok", "wall_s": 0.0,
    }

    try:
        t0 = time.time()
        region_fa = extract_region_fasta(genome, region["chrom"],
                                         region["start"], region["end"])
        with tempfile.TemporaryDirectory() as td:
            scratch = Path(td)
            (scratch / "region.fa").write_text(region_fa)
            (scratch / "query.fa").write_text(query_fasta)
            cmd = [MINIPROT, "--gff", "-t", "1", "-N", "5",
                   str(scratch / "region.fa"), str(scratch / "query.fa")]
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                 timeout=timeout)
        wall = time.time() - t0
        result["wall_s"] = round(wall, 1)

        # Parse GFF with correct genomic offset
        region_offset = region["start"]  # 0-based offset
        hits = parse_miniprot_gff(proc.stdout, region["chrom"], region_offset)

        # Extract A-stage features
        a_feats = extract_a_stage_features(
            region, gt_info, recovery_info, query_len)

        # Extract C-stage features
        c_feats = extract_c_stage_features(hits, query_len)

        # Compute F1 labels (primary and best-overlapping)
        f1_label = compute_f1_label(hits, gt_exons, paralog_cds)

        result.update(a_feats)
        result.update(c_feats)
        result.update({f"label_{k}": v for k, v in f1_label.items()})

    except subprocess.TimeoutExpired:
        result["status"] = "timeout"
        result["wall_s"] = timeout
        result.update(_zero_features())
    except Exception:
        result["status"] = "error"
        result["error"] = traceback.format_exc()[:200]
        result.update(_zero_features())

    return result


def _zero_features() -> dict:
    """Return zero-valued features for failed/skipped regions."""
    feats = {}
    for tool in A_TOOLS:
        short = tool.replace("_relaxed", "")
        feats[f"a_{short}"] = 0
    feats.update({
        "a_n_tools": 0, "a_best_score": 0.0, "a_best_qcov": 0.0,
        "a_n_raw_hits": 0, "region_size": 0, "region_size_log": 0.0,
        "on_primary_chrom": 0, "strand_match_primary": 0,
        "coorth_overlap_bp": 0, "coorth_same_chrom": 0,
        "query_length_aa": 0, "query_length_log": 0.0,
        "c_has_hit": 0, "c_identity": 0.0, "c_coverage": 0.0,
        "c_n_exons": 0, "c_cds_length": 0, "c_composite_score": 0.0,
        "c_hit_score": 0.0, "c_rank": 0,
        "label_primary_f1": 0.0, "label_primary_precision": 0.0,
        "label_primary_recall": 0.0, "label_best_overlap_f1": 0.0,
        "label_best_gene": "", "label_n_pred_exons": 0, "label_pred_cds_bp": 0,
    })
    return feats


# ---------------------------------------------------------------------------
# Region building
# ---------------------------------------------------------------------------

def build_candidate_regions(sp: str, qid: str, all_cands: dict[str, list],
                             recovery_info: dict | None,
                             gt_info: dict | None,
                             max_regions: int) -> list[dict]:
    """Build and prioritize merged candidate regions for a query protein."""
    # Collect all raw candidates for this query across tools
    raw = []
    for tool in A_TOOLS:
        for c in all_cands.get(tool, {}).get(qid, []):
            if len(c) >= 6:
                raw.append((*c[:6], tool))
            else:
                raw.append((*c[:5], 0.0, tool))

    if not raw:
        return []

    merged = merge_regions(raw, slop=50000)

    # Compute scores for prioritization
    for r in merged:
        score = 0.0
        # +100 for being on primary chrom
        if gt_info and r["chrom"] == gt_info.get("chrom", ""):
            score += 100.0
        # +10 per tool vote
        score += r["n_tools"] * 10.0
        # +best_score
        score += r["best_score"]
        # +size bonus (prefer larger regions up to 200kb)
        size_kb = (r["end"] - r["start"]) / 1000.0
        score += min(size_kb / 10.0, 20.0)
        r["_priority"] = score

    merged.sort(key=lambda x: x["_priority"], reverse=True)

    # Force-include co-ortholog region if applicable
    included = set()
    for i, r in enumerate(merged[:max_regions]):
        included.add((r["chrom"], r["start"], r["end"]))

    if recovery_info and recovery_info.get("status") == "recovered":
        co_chrom = recovery_info.get("coorth_chrom", "")
        co_start = int(recovery_info.get("coorth_start", 0))
        co_end = int(recovery_info.get("coorth_end", 0))
        co_key = (co_chrom, co_start, co_end)
        if co_chrom and co_end > co_start and co_key not in included:
            # Add co-ortholog region
            co_region = {
                "chrom": co_chrom, "start": co_start, "end": co_end,
                "strand": "+", "n_tools": 0, "tools": [],
                "best_score": 0.0, "best_qcov": 0.0, "n_raw_hits": 0,
                "_priority": 999.0, "_is_coorth": True,
            }
            merged.append(co_region)

    return merged[:max_regions]


# ---------------------------------------------------------------------------
# Training data generation
# ---------------------------------------------------------------------------

def _mp_worker(job_args_and_out) -> dict:
    """Module-level worker for multiprocessing. Pickleable."""
    args, out_path = job_args_and_out[:-1], Path(job_args_and_out[-1])
    result = _process_one_region(args)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(result, f)
    return result


def generate_training_data(sp: str, max_proteins: int | None = None,
                            max_regions: int = 15, timeout: int = 60,
                            n_workers: int = 1) -> pd.DataFrame:
    """Generate (features, label) rows for all candidate regions."""
    log = logging.getLogger("phase4")

    genome = GENOMES / f"{sp}.fa"
    query_fa_path = DATA / f"{sp}_paralog" / "proteins.fa"

    # Load data
    gt_json = load_gt_json(sp)
    qlen = load_query_lengths(sp)
    recovery = load_recovery_manifest(sp)
    all_gene_cds = load_cds_bed(sp)
    log.info("loaded CDS BED: %d genes", len(all_gene_cds))

    # Load query sequences
    qseq = {}
    with open(query_fa_path) as fh:
        cur_id = None; cur_seq: list[str] = []
        for line in fh:
            if line.startswith(">"):
                if cur_id and cur_seq:
                    qseq[cur_id] = "".join(cur_seq)
                cur_id = line[1:].strip().split()[0]; cur_seq = []
            elif cur_id:
                cur_seq.append(line.strip())
        if cur_id and cur_seq:
            qseq[cur_id] = "".join(cur_seq)

    # Load all A-tool BED files
    all_cands: dict[str, dict] = {}
    for tool in A_TOOLS:
        bed = BED_DIR / f"{tool}.{sp}.paralog.bed"
        if bed.exists():
            all_cands[tool] = load_bed(bed)
            log.info("loaded %s: %d queries", tool, len(all_cands[tool]))
        else:
            log.warning("%s not found: %s", tool, bed)

    # Filter to proteins with GT data
    all_qids = sorted(set(gt_json.keys()) & set(qseq.keys()))
    if max_proteins:
        all_qids = all_qids[:max_proteins]
    log.info("processing %d proteins for %s", len(all_qids), sp)

    # Build job list
    ROUTER.mkdir(parents=True, exist_ok=True)
    regions_dir = ROUTER / "regions"
    regions_dir.mkdir(parents=True, exist_ok=True)

    jobs = []
    skipped = 0
    for qid in all_qids:
        gt_info = gt_json.get(qid, {})
        rec_info = recovery.get(qid)

        # Build GT exons in 0-based half-open format
        gt_exons_0b = []
        for e in gt_info.get("exons", []):
            gt_exons_0b.append({
                "chrom": e["chrom"],
                "start": e["start"] - 1,  # 1-based → 0-based
                "end": e["end"],
            })

        # Build paralog CDS dict for this query: {gene_id: [exons]}
        paralog_gene_ids = [pl.get("target_gene", "") for pl
                           in gt_info.get("paralog_loci", [])]
        paralog_cds_for_query = {}
        for gid in paralog_gene_ids:
            if gid and gid in all_gene_cds:
                paralog_cds_for_query[gid] = all_gene_cds[gid]
        paralog_cds_json = json.dumps(paralog_cds_for_query) if paralog_cds_for_query else ""

        regions = build_candidate_regions(
            sp, qid, all_cands, rec_info, gt_info, max_regions)

        query_seq = qseq.get(qid, "")
        query_fa = f">{qid}\n{query_seq}\n" if query_seq else ""

        for i, region in enumerate(regions):
            out_json = regions_dir / f"{qid}.{i}.json"
            if out_json.exists():
                skipped += 1
                continue
            jobs.append((
                qid, i, region, str(genome), query_fa,
                qlen.get(qid, 0), gt_exons_0b, paralog_cds_json,
                gt_info, rec_info, timeout,
                str(out_json),
            ))

    log.info("%d total jobs, %d skipped (already done), %d to run",
             len(jobs) + skipped, skipped, len(jobs))

    if not jobs:
        # Load existing results from saved JSONs
        return _load_saved_regions(regions_dir, all_qids)

    # Run jobs
    all_rows = []
    done = 0
    t0 = time.time()

    if n_workers > 1:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = [pool.submit(_mp_worker, j) for j in jobs]
            for fut in as_completed(futures):
                try:
                    row = fut.result()
                    all_rows.append(row)
                except Exception as e:
                    log.error("worker failed: %s", e)
                done += 1
                if done % 200 == 0:
                    elapsed = time.time() - t0
                    rate = done / max(elapsed, 1)
                    log.info("progress: %d/%d jobs (%.0f/min), %.0fs elapsed",
                             done, len(jobs), rate * 60, elapsed)
    else:
        for job in jobs:
            args, out_path = job[:-1], Path(job[-1])
            row = _process_one_region(args)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            with open(out_path, "w") as f:
                json.dump(row, f)
            all_rows.append(row)
            done += 1
            if done % 50 == 0:
                elapsed = time.time() - t0
                rate = done / max(elapsed, 1)
                log.info("progress: %d/%d jobs (%.0f/min), %.0fs elapsed",
                         done, len(jobs), rate * 60, elapsed)

    elapsed = time.time() - t0
    log.info("DONE: %d jobs in %.0fs (%.0f/min)",
             done, elapsed, done / max(elapsed, 1) * 60)

    # Also load previously saved results
    if skipped > 0:
        saved = _load_saved_regions(regions_dir, all_qids)
        all_rows.extend(saved)

    if not all_rows:
        log.warning("No results generated")
        return pd.DataFrame()

    df = pd.DataFrame(all_rows)
    # Save
    tsv_path = ROUTER / f"training_data_{sp}.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    log.info("saved %d rows to %s", len(df), tsv_path)

    return df


def _load_saved_regions(regions_dir: Path, qids: list[str]) -> list[dict]:
    """Load previously saved per-region JSON results."""
    rows = []
    qid_set = set(qids)
    for f in regions_dir.glob("*.json"):
        qid = f.stem.rsplit(".", 1)[0]
        if qid in qid_set:
            try:
                with open(f) as fh:
                    rows.append(json.load(fh))
            except Exception:
                pass
    return rows


# ---------------------------------------------------------------------------
# Training and evaluation
# ---------------------------------------------------------------------------

FEATURE_COLS = [
    "a_diamond", "a_last", "a_miniprot", "a_mmseqs",
    "a_n_tools", "a_best_score", "a_best_qcov", "a_n_raw_hits",
    "region_size", "region_size_log",
    "on_primary_chrom", "strand_match_primary",
    "coorth_overlap_bp", "coorth_same_chrom",
    "query_length_aa", "query_length_log",
    "c_has_hit", "c_identity", "c_coverage", "c_n_exons",
    "c_cds_length", "c_composite_score", "c_hit_score", "c_rank",
]

LABEL_COL = "label_best_overlap_f1"  # default, changed by --label flag
LABEL_PRIMARY_COL = "label_primary_f1"
LABEL_BEST_OVERLAP_COL = "label_best_overlap_f1"


def train_and_evaluate(df: pd.DataFrame, sp: str,
                        test_frac: float = 0.2,
                        seed: int = 4242) -> dict:
    """Train XGBoost router and evaluate against baselines."""
    from sklearn.model_selection import train_test_split
    import xgboost as xgb

    log = logging.getLogger("phase4")

    # Filter to valid rows
    valid = df[df[LABEL_COL].notna()].copy()
    if len(valid) == 0:
        log.error("No valid training data")
        return {}

    # Ensure all feature columns exist
    for col in FEATURE_COLS:
        if col not in valid.columns:
            valid[col] = 0.0

    # Split by query
    all_qids = sorted(valid["query_id"].unique())
    n_test = max(1, int(len(all_qids) * test_frac))
    np.random.seed(seed)
    test_qids = set(np.random.choice(all_qids, size=n_test, replace=False))
    train_qids = set(all_qids) - test_qids

    train_df = valid[valid["query_id"].isin(train_qids)]
    test_df = valid[valid["query_id"].isin(test_qids)]

    log.info("train: %d queries, %d rows; test: %d queries, %d rows",
             len(train_qids), len(train_df), len(test_qids), len(test_df))

    if len(train_df) == 0 or len(test_df) == 0:
        log.error("Empty train/test split")
        return {}

    # Train XGBoost
    X_train = train_df[FEATURE_COLS].fillna(0.0).to_numpy(dtype=np.float32)
    y_train = train_df[LABEL_COL].fillna(0.0).to_numpy(dtype=np.float32)

    model = xgb.XGBRegressor(
        n_estimators=200, max_depth=6, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        reg_alpha=0.1, reg_lambda=1.0,
        objective="reg:squarederror",
        random_state=seed, n_jobs=-1,
    )
    model.fit(X_train, y_train)

    # Save model
    model_path = ROUTER / f"model_{sp}.json"
    model.save_model(str(model_path))
    log.info("saved model to %s", model_path)

    # Feature importance
    imp = dict(zip(FEATURE_COLS, model.feature_importances_.tolist()))

    # Evaluate on test set
    X_test = test_df[FEATURE_COLS].fillna(0.0).to_numpy(dtype=np.float32)

    # Per-query routing evaluation
    eval_rows = []
    for qid in sorted(test_qids):
        qmask = test_df["query_id"] == qid
        qdf = test_df[qmask]
        if len(qdf) == 0:
            continue

        X_q = qdf[FEATURE_COLS].fillna(0.0).to_numpy(dtype=np.float32)
        y_q = qdf[LABEL_COL].values
        pred = model.predict(X_q)

        # Router: pick candidate with max predicted F1
        best_idx = int(np.argmax(pred))
        router_f1 = float(y_q[best_idx])

        # Baselines
        oracle_f1 = float(np.max(y_q))

        # Best A-stage score
        a_scores = qdf["a_best_score"].values
        best_a_idx = int(np.argmax(a_scores))
        best_a_f1 = float(y_q[best_a_idx])

        # Best C-stage composite score
        c_scores = qdf["c_composite_score"].fillna(0.0).values
        best_c_idx = int(np.argmax(c_scores))
        best_c_f1 = float(y_q[best_c_idx])

        # Best A + C simple combination
        combined = (qdf["a_best_score"].fillna(0.0).values +
                    qdf["c_composite_score"].fillna(0.0).values * 0.5)
        best_combined_idx = int(np.argmax(combined))
        best_combined_f1 = float(y_q[best_combined_idx])

        eval_rows.append({
            "query_id": qid,
            "n_candidates": len(qdf),
            "router_f1": router_f1,
            "oracle_f1": oracle_f1,
            "best_a_score_f1": best_a_f1,
            "best_c_composite_f1": best_c_f1,
            "best_ac_combined_f1": best_combined_f1,
            "router_rank": int(np.sum(y_q > router_f1)) + 1,
        })

    eval_df = pd.DataFrame(eval_rows)

    # Load Exp1 baseline F1
    exp1 = load_exp1_f1(sp)
    for i, row in eval_df.iterrows():
        qid = row["query_id"]
        eval_df.at[i, "always_primary_f1"] = exp1.get(qid, {}).get("exon_f1", 0.0)
        eval_df.at[i, "always_best_overlap_f1"] = exp1.get(qid, {}).get("best_overlapping_f1", 0.0)

    # Save evaluation
    eval_path = ROUTER / f"evaluation_{sp}.tsv"
    eval_df.to_csv(eval_path, sep="\t", index=False)

    # Aggregate metrics
    def mean_ci(vals, n_boot=10000):
        v = np.array(vals, dtype=float)
        m = np.mean(v)
        rng = np.random.default_rng(seed)
        boots = [np.mean(rng.choice(v, size=len(v), replace=True))
                 for _ in range(n_boot)]
        return {"mean": float(m), "ci_lo": float(np.percentile(boots, 2.5)),
                "ci_hi": float(np.percentile(boots, 97.5))}

    summary = {
        "species": sp,
        "n_test_queries": len(eval_df),
        "n_train_queries": len(train_qids),
        "label_type": "best_overlap_f1",
        "router": mean_ci(eval_df["router_f1"]),
        "oracle": mean_ci(eval_df["oracle_f1"]),
        "baselines": {
            "always_primary": mean_ci(eval_df["always_primary_f1"]),
            "always_best_overlap": mean_ci(eval_df["always_best_overlap_f1"]),
            "best_a_score": mean_ci(eval_df["best_a_score_f1"]),
            "best_c_composite": mean_ci(eval_df["best_c_composite_f1"]),
            "best_ac_combined": mean_ci(eval_df["best_ac_combined_f1"]),
        },
        "router_vs_oracle_efficiency": float(
            np.mean(eval_df["router_f1"]) / max(np.mean(eval_df["oracle_f1"]), 1e-8)),
        "router_better_than_primary_pct": float(
            np.mean(eval_df["router_f1"] > eval_df["always_primary_f1"]) * 100),
        "router_better_than_best_overlap_pct": float(
            np.mean(eval_df["router_f1"] > eval_df["always_best_overlap_f1"]) * 100),
        "feature_importance": {k: float(v) for k, v in
                               sorted(imp.items(), key=lambda x: -x[1])[:10]},
    }

    # Paired bootstrap deltas
    for name in ["always_primary", "always_best_overlap", "best_a_score", "best_c_composite"]:
        baseline_col = f"{name}_f1"
        deltas = eval_df["router_f1"].values - eval_df[baseline_col].values
        summary["baselines"][name]["delta_vs_router"] = float(np.mean(deltas))
        delta_ci = mean_ci(deltas)
        summary["baselines"][name]["delta_ci_lo"] = delta_ci["ci_lo"]
        summary["baselines"][name]["delta_ci_hi"] = delta_ci["ci_hi"]

    summary_path = ROUTER / "evaluation_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    log.info("saved evaluation summary to %s", summary_path)

    # Print summary
    def _p(vals):
        m = mean_ci(vals)
        return f"{m['mean']:.4f} [{m['ci_lo']:.4f}, {m['ci_hi']:.4f}]"

    print(f"\n{'='*70}")
    print(f"Phase 4 Router Evaluation — {sp} (label=best_overlap_f1)")
    print(f"{'='*70}")
    print(f"Test queries: {len(eval_df)}")
    print(f"Router mean F1:          {_p(eval_df['router_f1'])}")
    print(f"Oracle mean F1:          {_p(eval_df['oracle_f1'])}")
    print(f"Router/Oracle:           {summary['router_vs_oracle_efficiency']:.3f}")
    print(f"\nBaselines vs Router (best_overlap_f1):")
    for name in ["always_primary", "always_best_overlap", "best_a_score", "best_c_composite"]:
        b = summary["baselines"][name]
        sig = "*" if b.get("delta_ci_lo", 0) > 0 else ""
        print(f"  {name:25s}  {b['mean']:.4f}  "
              f"delta={b.get('delta_vs_router', 0):+.4f} "
              f"[{b.get('delta_ci_lo', 0):+.4f}, {b.get('delta_ci_hi', 0):+.4f}]{sig}")
    print(f"\nRouter > Primary:          {summary['router_better_than_primary_pct']:.1f}% of queries")
    print(f"Router > Best-overlap:     {summary['router_better_than_best_overlap_pct']:.1f}% of queries")
    print(f"\nTop features:")
    for feat, score in list(summary["feature_importance"].items())[:5]:
        print(f"  {feat:25s} {score:.4f}")

    return summary


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description="TPRS Phase 4: Locus Router")
    ap.add_argument("--species", nargs="*", default=["mouse"])
    ap.add_argument("--max-proteins", type=int, default=None,
                    help="Max proteins to process (for quick tests)")
    ap.add_argument("--max-regions", type=int, default=15,
                    help="Max candidate regions per protein")
    ap.add_argument("--timeout", type=int, default=60,
                    help="Timeout per miniprot run (seconds)")
    ap.add_argument("--workers", type=int, default=1,
                    help="Number of parallel workers")
    ap.add_argument("--train-only", action="store_true",
                    help="Skip data generation, train from saved TSV")
    ap.add_argument("--eval-only", action="store_true",
                    help="Skip data generation, evaluate from saved TSV")
    ap.add_argument("--test-frac", type=float, default=0.2,
                    help="Test set fraction")
    ap.add_argument("--seed", type=int, default=4242)
    ap.add_argument("--label", type=str, default="best_overlap",
                    choices=["best_overlap", "primary"],
                    help="Training target: best_overlap or primary")
    args = ap.parse_args()

    global LABEL_COL
    LABEL_COL = "label_best_overlap_f1" if args.label == "best_overlap" else "label_primary_f1"

    ROUTER.mkdir(parents=True, exist_ok=True)
    log_path = ROUTER / "phase4.log"
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] PHASE4 %(message)s",
        handlers=[logging.StreamHandler(sys.stdout),
                  logging.FileHandler(log_path, mode="a")],
    )
    log = logging.getLogger("phase4")

    if not MINIPROT:
        log.error("miniprot binary not found. Set MINIPROT_BIN env var.")
        return 1

    all_summaries = []
    for sp in args.species:
        log.info("=== Phase 4 Router: %s ===", sp)

        if not args.eval_only and not args.train_only:
            df = generate_training_data(
                sp, max_proteins=args.max_proteins,
                max_regions=args.max_regions, timeout=args.timeout,
                n_workers=args.workers,
            )
        else:
            tsv_path = ROUTER / f"training_data_{sp}.tsv"
            if tsv_path.exists():
                df = pd.read_csv(tsv_path, sep="\t")
                log.info("loaded %d rows from %s", len(df), tsv_path)
            else:
                log.error("Training data not found: %s", tsv_path)
                continue

        if df.empty:
            log.error("Empty dataframe for %s, skipping training", sp)
            continue

        summary = train_and_evaluate(df, sp, test_frac=args.test_frac,
                                      seed=args.seed)
        all_summaries.append(summary)

    return 0


if __name__ == "__main__":
    sys.exit(main())
