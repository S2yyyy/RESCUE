"""Phase 5 Step 1-3 driver: run all Phase 5 steps for a given species.

Usage: python scripts/tprs/phase5_run.py <species>
"""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupShuffleSplit
from xgboost import XGBRegressor

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import (  # noqa: E402
    GENE_STRUCTURES,
    GENOMES,
    MINIPROT,
    ORTHOLOGS,
    PROJECT_ROOT,
    PROTEOMES,
    SAMTOOLS,
    V3_DATA,
)

REPO = PROJECT_ROOT
DATA = V3_DATA
GENESTRUCT = GENE_STRUCTURES
RECOVERY = REPO / "results/tprs/recovery"
ROUTER = REPO / "results/tprs/router"
OUTDIR = REPO / "results/tprs/phase5"

RANDOM_SEED = 4242
N_BOOTSTRAP = 10_000

FEATURE_COLS = [
    "a_diamond", "a_last", "a_miniprot", "a_mmseqs",
    "a_n_tools", "a_best_score", "a_best_qcov", "a_n_raw_hits",
    "region_size", "region_size_log",
    "on_primary_chrom", "strand_match_primary",
    "coorth_overlap_bp", "coorth_same_chrom",
    "query_length_aa", "query_length_log",
    "c_has_hit", "c_identity", "c_coverage", "c_n_exons",
    "c_cds_length", "c_composite_score", "c_hit_score", "c_rank",
    "query_type_num",
]

SPECIES_ORTHOLOG = {
    "mouse": "mus_musculus",
    "dog": "canis_lupus_familiaris",
    "chicken": "gallus_gallus",
    "zebrafish": "danio_rerio",
}


# ---------------------------------------------------------------------------
# Data loading (step 1 utilities)
# ---------------------------------------------------------------------------

def load_proteome_by_gene(sp: str) -> dict[str, str]:
    path = PROTEOMES / f"{sp}_canonical.fa"
    if not path.exists():
        return {}
    gene2prot = {}
    with open(path) as fh:
        cur_gene, cur_seq = None, []
        for line in fh:
            if line.startswith(">"):
                if cur_gene and cur_seq:
                    gene2prot[cur_gene] = "".join(cur_seq)
                cur_seq = []
                header = line[1:].strip()
                gene_part = [p for p in header.split() if p.startswith("gene:")]
                cur_gene = gene_part[0].split(":")[1].split(".")[0] if gene_part else None
            elif cur_gene:
                cur_seq.append(line.strip())
        if cur_gene and cur_seq:
            gene2prot[cur_gene] = "".join(cur_seq)
    return gene2prot


def load_cds_bed(sp: str) -> dict[str, list[dict]]:
    cds = defaultdict(list)
    bed_path = GENESTRUCT / f"{sp}_CDS.bed"
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.strip().split("\t")
            gene_id = f[3].rsplit("_", 1)[0] if "_" in f[3] else f[3]
            cds[gene_id].append({
                "chrom": f[0], "start": int(f[1]), "end": int(f[2]),
                "strand": f[5],
            })
    return dict(cds)


def load_ortholog_map(sp: str) -> dict[str, str]:
    sp_name = SPECIES_ORTHOLOG.get(sp, "mus_musculus")
    path = ORTHOLOGS / f"human_vs_{sp_name}_one2one.tsv"
    if not path.exists():
        return {}
    human2gene = {}
    with open(path) as fh:
        header = next(fh).strip().split("\t")
        gene_idx = header.index("gene_stable_id")
        prot_idx = header.index("homology_protein_stable_id")
        for line in fh:
            if not line.strip():
                continue
            f = line.strip().split("\t")
            human2gene[f[prot_idx]] = f[gene_idx]
    return human2gene


def build_gene_index(cds_bed: dict) -> dict:
    index = defaultdict(list)
    for gene_id, exons in cds_bed.items():
        chrom = exons[0]["chrom"]
        gstart = min(e["start"] for e in exons)
        gend = max(e["end"] for e in exons)
        for bin_start in range(gstart // 1_000_000 * 1_000_000,
                               gend + 1_000_000, 1_000_000):
            index[(chrom, bin_start)].append((gene_id, gstart, gend))
    return index


def find_best_gene(region_chrom: str, region_start: int, region_end: int,
                   gene_index: dict, primary_gene: str = "") -> tuple[str, str]:
    best_gene, best_overlap = "", 0
    for bin_start in range(region_start // 1_000_000 * 1_000_000,
                           region_end + 1_000_000, 1_000_000):
        for gene_id, gstart, gend in gene_index.get((region_chrom, bin_start), []):
            overlap = max(0, min(region_end, gend) - max(region_start, gstart))
            if overlap > best_overlap:
                best_overlap, best_gene = overlap, gene_id
            if primary_gene and gene_id == primary_gene and overlap > 0:
                return gene_id, "primary_ortholog"
    if best_gene:
        return best_gene, "target_gene_overlap"
    return "", ""


# ---------------------------------------------------------------------------
# Miniprot execution (step 2 utilities)
# ---------------------------------------------------------------------------

def extract_region(genome: Path, chrom: str, start_0b: int, end_0b: int) -> str:
    region = f"{chrom}:{start_0b + 1}-{end_0b}"
    proc = subprocess.run(
        [SAMTOOLS, "faidx", str(genome), region],
        capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"samtools failed: {proc.stderr[:200]}")
    return proc.stdout


def run_miniprot(genome: Path, chrom: str, start_0b: int, end_0b: int,
                 query_fa: str, tmp_dir: Path,
                 extra_args: list[str] | None = None,
                 timeout: int = 120) -> tuple[str, int]:
    seq = extract_region(genome, chrom, start_0b, end_0b)
    (tmp_dir / "region.fa").write_text(seq)
    (tmp_dir / "query.fa").write_text(query_fa)
    cmd = [MINIPROT, "--gff", "-t", "1", "-N", "3"]
    if extra_args:
        cmd.extend(extra_args)
    cmd.extend([str(tmp_dir / "region.fa"), str(tmp_dir / "query.fa")])
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    return proc.stdout, start_0b


def parse_miniprot_gff(gff_text: str, region_chrom: str,
                       region_offset: int) -> list[dict]:
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


# ---------------------------------------------------------------------------
# Feature extraction & F1
# ---------------------------------------------------------------------------

def extract_c_stage_features(hits: list[dict], query_len: int) -> dict:
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
    exon_bonus = min(1.0, len(cds) / 3.0)
    feats["c_composite_score"] = round(
        feats["c_identity"] * min(1.0, coverage) * exon_bonus, 4)
    return feats


def base_level_metrics(pred_cds: list[dict], gt_cds: list[dict]) -> dict:
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
    overlap = min(overlap, pred_bp, gt_bp)
    precision = overlap / pred_bp if pred_bp > 0 else 0.0
    recall = overlap / gt_bp if gt_bp > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
    return {"precision": precision, "recall": recall, "f1": f1}


def compute_f1_label(pred_hits: list[dict], gt_primary_exons: list[dict],
                     paralog_cds: dict[str, list[dict]] | None = None) -> dict:
    if not pred_hits:
        return {"primary_f1": 0.0, "best_overlap_f1": 0.0,
                "primary_precision": 0.0, "primary_recall": 0.0,
                "best_gene": "", "n_pred_exons": 0, "pred_cds_bp": 0}
    best = pred_hits[0]
    pred_cds = best.get("cds", [])
    p_metrics = base_level_metrics(pred_cds, gt_primary_exons)
    best_overlap = p_metrics["f1"]
    best_gene = "primary"
    if paralog_cds:
        for gene_id, gene_exons in paralog_cds.items():
            m = base_level_metrics(pred_cds, gene_exons)
            if m["f1"] > best_overlap:
                best_overlap = m["f1"]
                best_gene = gene_id
    return {
        "primary_f1": p_metrics["f1"],
        "best_overlap_f1": best_overlap,
        "primary_precision": p_metrics["precision"],
        "primary_recall": p_metrics["recall"],
        "best_gene": best_gene,
        "n_pred_exons": len(pred_cds),
        "pred_cds_bp": sum(e["end"] - e["start"] for e in pred_cds),
    }


# ---------------------------------------------------------------------------
# Worker: process one region with species protein
# ---------------------------------------------------------------------------

def process_region(args: tuple) -> dict:
    (row_dict, protein_map, gene2prot, gt_data, cds_bed,
     genome_path, regions_dir, query_lens) = args

    qid = row_dict["query_id"]
    ridx = row_dict["region_idx"]
    chrom = str(row_dict["chrom"])
    rstart = int(row_dict["start"])
    rend = int(row_dict["end"])

    out_json = regions_dir / f"{qid}.{ridx}.json"
    if out_json.exists():
        with open(out_json) as f:
            cached = json.load(f)
        if cached.get("status") in ("ok", "no_protein", "no_hits", "no_hits_fallback"):
            return cached

    region_key = f"{qid}|{chrom}|{rstart}|{rend}"
    pm = protein_map.get(region_key, {})
    best_gene = pm.get("best_gene", "")
    has_protein = pm.get("has_protein", False)

    if not best_gene or not has_protein:
        result = {"query_id": qid, "region_idx": ridx, "status": "no_protein"}
        out_json.parent.mkdir(parents=True, exist_ok=True)
        with open(out_json, "w") as f:
            json.dump(result, f)
        return result

    mouse_seq = gene2prot.get(best_gene, "")
    if not mouse_seq:
        result = {"query_id": qid, "region_idx": ridx,
                  "status": "no_protein", "best_gene": best_gene}
        out_json.parent.mkdir(parents=True, exist_ok=True)
        with open(out_json, "w") as f:
            json.dump(result, f)
        return result

    gt = gt_data.get(qid, {})
    gt_primary_exons = gt.get("exons", [])
    paralog_genes = {}
    for pl in gt.get("paralog_loci", []):
        pg = pl.get("target_gene", "")
        if pg and pg in cds_bed and pg not in paralog_genes:
            paralog_genes[pg] = cds_bed[pg]

    mouse_aa_len = len(mouse_seq)
    GENE_BUFFER = 20_000
    gene_exons = cds_bed.get(best_gene, [])
    if gene_exons:
        gene_start = min(e["start"] for e in gene_exons)
        gene_end = max(e["end"] for e in gene_exons)
        extract_start = min(rstart, gene_start - GENE_BUFFER)
        extract_end = max(rend, gene_end + GENE_BUFFER)
    else:
        extract_start, extract_end = rstart, rend

    result = {
        "query_id": qid, "region_idx": ridx,
        "chrom": chrom, "start": rstart, "end": rend,
        "best_gene": best_gene,
        "mouse_protein_len_aa": mouse_aa_len,
        "extract_start": extract_start, "extract_end": extract_end,
    }

    query_fa = f">{qid}_sp@{best_gene}\n{mouse_seq}\n"
    try:
        with tempfile.TemporaryDirectory() as td:
            gff_out, offset = run_miniprot(
                genome_path, chrom, extract_start, extract_end,
                query_fa, Path(td), timeout=120)
            hits = parse_miniprot_gff(gff_out, chrom, offset)
    except subprocess.TimeoutExpired:
        result["status"] = "timeout"
        out_json.parent.mkdir(parents=True, exist_ok=True)
        with open(out_json, "w") as f:
            json.dump(result, f)
        return result
    except Exception as e:
        result["status"] = "error"
        result["error"] = str(e)[:200]
        out_json.parent.mkdir(parents=True, exist_ok=True)
        with open(out_json, "w") as f:
            json.dump(result, f)
        return result

    if not hits:
        try:
            with tempfile.TemporaryDirectory() as td:
                gff_out2, offset2 = run_miniprot(
                    genome_path, chrom, extract_start, extract_end,
                    query_fa, Path(td),
                    extra_args=["-p", "0.0"], timeout=120)
                hits2 = parse_miniprot_gff(gff_out2, chrom, offset2)
                if hits2:
                    hits = hits2
                    result["fallback_used"] = "p0.0"
        except Exception:
            pass

    c_feats = extract_c_stage_features(hits, mouse_aa_len)
    labels = compute_f1_label(hits, gt_primary_exons, paralog_genes)

    result.update({
        "status": "ok" if hits else "no_hits",
        "fallback_used": result.get("fallback_used", ""),
        **c_feats,
        **{f"label_{k}": v for k, v in labels.items()},
    })

    out_json.parent.mkdir(parents=True, exist_ok=True)
    with open(out_json, "w") as f:
        json.dump(result, f)
    return result


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def paired_bootstrap_ci(values: np.ndarray, n_bootstrap: int = N_BOOTSTRAP,
                        ci: float = 0.95, seed: int = RANDOM_SEED) -> tuple:
    rng = np.random.RandomState(seed)
    n = len(values)
    means = np.array([rng.choice(values, size=n, replace=True).mean()
                      for _ in range(n_bootstrap)])
    alpha = (1 - ci) / 2
    lower, upper = np.percentile(means, [100 * alpha, 100 * (1 - alpha)])
    return float(values.mean()), float(lower), float(upper)


def evaluate_router(model: XGBRegressor, test_df: pd.DataFrame,
                    label_col: str, eval_label_col: str) -> dict:
    queries = test_df["query_id"].unique()
    actual_f1s = []
    actual_primary_f1s = []

    for qid in queries:
        qdf = test_df[test_df["query_id"] == qid]
        if qdf.empty:
            continue
        X = qdf[FEATURE_COLS].fillna(0.0)
        preds = model.predict(X)
        best_row = qdf.iloc[np.argmax(preds)]
        actual_f1s.append(best_row[eval_label_col])
        actual_primary_f1s.append(best_row["label_primary_f1"])

    actual_f1s = np.array(actual_f1s)
    actual_primary_f1s = np.array(actual_primary_f1s)
    mean_best, lo_best, hi_best = paired_bootstrap_ci(actual_f1s)
    mean_prim, lo_prim, hi_prim = paired_bootstrap_ci(actual_primary_f1s)

    return {
        "n_queries": len(queries),
        "best_overlap_mean": mean_best,
        "best_overlap_ci_lower": lo_best,
        "best_overlap_ci_upper": hi_best,
        "primary_mean": mean_prim,
        "primary_ci_lower": lo_prim,
        "primary_ci_upper": hi_prim,
        "zero_f1_count": int((actual_f1s == 0).sum()),
    }


def evaluate_baseline(test_df: pd.DataFrame, baseline_type: str) -> dict:
    queries = test_df["query_id"].unique()
    f1s, primary_f1s = [], []

    for qid in queries:
        qdf = test_df[test_df["query_id"] == qid]
        if qdf.empty:
            continue

        if baseline_type == "always_primary":
            subset = qdf[qdf["query_type"] == "human"]
            if subset.empty:
                subset = qdf
            best_idx = subset["on_primary_chrom"].values.argmax()
            if (subset["on_primary_chrom"].values == 1).sum() > 1:
                primary_sub = subset[subset["on_primary_chrom"] == 1]
                best_idx = primary_sub["coorth_overlap_bp"].values.argmax()
                best_row = primary_sub.iloc[best_idx]
            else:
                best_row = subset.iloc[best_idx]

        elif baseline_type == "best_c_composite":
            subset = qdf[qdf["query_type"] == "human"]
            if subset.empty:
                subset = qdf
            best_row = subset.iloc[subset["c_composite_score"].values.argmax()]

        elif baseline_type == "best_a_score":
            subset = qdf[qdf["query_type"] == "human"]
            if subset.empty:
                subset = qdf
            best_row = subset.iloc[subset["a_best_score"].values.argmax()]

        elif baseline_type == "oracle":
            best_row = qdf.iloc[qdf["label_best_overlap_f1"].values.argmax()]

        else:
            continue

        f1s.append(best_row["label_best_overlap_f1"])
        primary_f1s.append(best_row["label_primary_f1"])

    f1s = np.array(f1s)
    primary_f1s = np.array(primary_f1s)
    mean_best, lo_best, hi_best = paired_bootstrap_ci(f1s)
    mean_prim, lo_prim, hi_prim = paired_bootstrap_ci(primary_f1s)

    return {
        "n_queries": len(queries),
        "best_overlap_mean": mean_best,
        "best_overlap_ci_lower": lo_best,
        "best_overlap_ci_upper": hi_best,
        "primary_mean": mean_prim,
        "primary_ci_lower": lo_prim,
        "primary_ci_upper": hi_prim,
        "zero_f1_count": int((f1s == 0).sum()),
    }


# ---------------------------------------------------------------------------
# Main per-species pipeline
# ---------------------------------------------------------------------------

def run_species(sp: str, workers: int = 8):
    sp_outdir = OUTDIR / sp
    sp_outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Phase 5: {sp}")
    print(f"{'='*60}")

    # ---- Step 1: Build candidate-to-protein mapping ----
    print("\n--- Step 1: Building candidate-to-protein mapping ---")
    cds_bed = load_cds_bed(sp)
    gene2prot = load_proteome_by_gene(sp)
    human2sp = load_ortholog_map(sp)
    phase4_path = ROUTER / f"training_data_{sp}.tsv"
    df = pd.read_csv(phase4_path, sep="\t")

    print(f"CDS genes: {len(cds_bed)}, proteome: {len(gene2prot)}, "
          f"ortholog pairs: {len(human2sp)}")
    print(f"Training data: {len(df)} rows, {df['query_id'].nunique()} queries")

    gt_path = DATA / f"{sp}_paralog" / "ground_truth.json"
    with open(gt_path) as f:
        gt_data = json.load(f)

    gene_index = build_gene_index(cds_bed)

    region_cols = ["query_id", "chrom", "start", "end"]
    regions = df[region_cols].drop_duplicates()
    print(f"Unique regions: {len(regions)}")

    mapping = {}
    stats = {"total": len(regions), "gene_found": 0, "protein_found": 0}

    for _, row in regions.iterrows():
        qid = row["query_id"]
        chrom = str(row["chrom"])
        rstart = int(row["start"])
        rend = int(row["end"])
        region_key = f"{qid}|{chrom}|{rstart}|{rend}"

        gt = gt_data.get(qid, {})
        primary_gene = gt.get("target_gene", "")

        best_gene, source = find_best_gene(
            chrom, rstart, rend, gene_index, primary_gene)

        if not best_gene:
            ortholog_gene = human2sp.get(qid, "")
            if ortholog_gene and ortholog_gene in cds_bed:
                og_exons = cds_bed[ortholog_gene]
                og_start = min(e["start"] for e in og_exons)
                og_end = max(e["end"] for e in og_exons)
                if max(0, min(rend, og_end) - max(rstart, og_start)) > 0:
                    best_gene = ortholog_gene
                    source = "primary_ortholog"

        protein_seq = ""
        if best_gene:
            stats["gene_found"] += 1
            protein_seq = gene2prot.get(best_gene, "")
            if protein_seq:
                stats["protein_found"] += 1
            else:
                ortholog_gene = human2sp.get(qid, "")
                if ortholog_gene and ortholog_gene != best_gene:
                    protein_seq = gene2prot.get(ortholog_gene, "")
                    if protein_seq:
                        best_gene = ortholog_gene
                        source = "primary_ortholog_fallback"
                        stats["protein_found"] += 1

        mapping[region_key] = {
            "query_id": qid, "chrom": chrom, "start": rstart,
            "end": rend, "best_gene": best_gene,
            "protein_source": source,
            "has_protein": len(protein_seq) > 0,
            "protein_length": len(protein_seq),
        }

    print(f"Gene overlap: {stats['gene_found']}/{stats['total']} "
          f"({100*stats['gene_found']/stats['total']:.1f}%)")
    print(f"Protein found: {stats['protein_found']}/{stats['total']} "
          f"({100*stats['protein_found']/stats['total']:.1f}%)")

    with open(sp_outdir / "candidate_protein_map.json", "w") as f:
        json.dump({"stats": stats, "mapping": mapping}, f, indent=2)

    # ---- Step 2: Run miniprot with species-matched proteins ----
    print("\n--- Step 2: Running miniprot with species-matched proteins ---")
    genome_path = GENOMES / f"{sp}.fa"
    regions_dir = sp_outdir / "regions"
    regions_dir.mkdir(parents=True, exist_ok=True)

    row_dicts = df.to_dict("records")
    tasks = [
        (row_dicts[i], mapping, gene2prot, gt_data, cds_bed,
         genome_path, regions_dir, None)
        for i in range(len(row_dicts))
    ]

    print(f"Processing {len(tasks)} regions with {workers} workers...")
    t0 = time.time()

    with Pool(workers) as pool:
        results = pool.map(process_region, tasks)

    elapsed = time.time() - t0
    ok_count = sum(1 for r in results if r.get("status") == "ok")
    print(f"Done in {elapsed:.0f}s. ok={ok_count}/{len(results)}")

    # Build combined training data
    a_cols = [c for c in df.columns if c.startswith("a_")]
    region_cols_list = ["region_size", "region_size_log", "on_primary_chrom",
                        "strand_match_primary", "coorth_overlap_bp", "coorth_same_chrom"]
    query_cols_list = ["query_length_aa", "query_length_log"]
    c_cols = [c for c in df.columns if c.startswith("c_")]
    label_cols = [c for c in df.columns if c.startswith("label_")]
    base_cols = ["query_id", "region_idx", "chrom", "start", "end", "strand"]

    human_df = df[base_cols + a_cols + region_cols_list + query_cols_list + c_cols + label_cols].copy()
    human_df["query_type"] = "human"

    result_by_key = {}
    for r in results:
        result_by_key[(r.get("query_id"), r.get("region_idx"))] = r

    mouse_rows = []
    for i, row in df.iterrows():
        qid = row["query_id"]
        ridx = row["region_idx"]
        r = result_by_key.get((qid, ridx), {})

        if r.get("status") not in ("ok", "no_hits", "no_hits_fallback"):
            mouse_rows.append({
                "query_id": qid, "region_idx": ridx,
                "chrom": row["chrom"], "start": row["start"], "end": row["end"],
                "strand": row["strand"],
                **{c: row[c] for c in a_cols},
                **{c: row[c] for c in region_cols_list},
                **{c: row[c] for c in query_cols_list},
                "c_has_hit": 0, "c_identity": 0.0, "c_coverage": 0.0,
                "c_n_exons": 0, "c_cds_length": 0, "c_composite_score": 0.0,
                "c_hit_score": 0.0, "c_rank": 99,
                "label_primary_f1": 0.0, "label_best_overlap_f1": 0.0,
                "label_primary_precision": 0.0, "label_primary_recall": 0.0,
                "label_best_gene": "", "label_n_pred_exons": 0, "label_pred_cds_bp": 0,
                "query_type": "human",
            })
        else:
            mouse_rows.append({
                "query_id": qid, "region_idx": ridx,
                "chrom": row["chrom"], "start": row["start"], "end": row["end"],
                "strand": row["strand"],
                **{c: row[c] for c in a_cols},
                **{c: row[c] for c in region_cols_list},
                "query_length_aa": r.get("mouse_protein_len_aa", 0),
                "query_length_log": np.log10(max(1, r.get("mouse_protein_len_aa", 1))),
                **{f"c_{k}": r.get(f"c_{k}", 0.0 if k != "rank" else 99)
                   for k in ["has_hit", "identity", "coverage", "n_exons",
                             "cds_length", "composite_score", "hit_score", "rank"]},
                **{f"label_{k}": r.get(f"label_{k}", 0.0 if "f1" in k else ("" if k == "best_gene" else 0))
                   for k in ["primary_f1", "best_overlap_f1", "primary_precision",
                             "primary_recall", "best_gene", "n_pred_exons", "pred_cds_bp"]},
                "query_type": "mouse",
            })

    mouse_df = pd.DataFrame(mouse_rows)
    combined = pd.concat([human_df, mouse_df], ignore_index=True)
    combined["query_type_num"] = (combined["query_type"] == "mouse").astype(int)

    out_path = sp_outdir / "training_data.tsv"
    combined.to_csv(out_path, sep="\t", index=False, float_format="%.6g")
    print(f"Combined training data: {len(combined)} rows → {out_path}")

    # ---- Step 3: Train dual routers and evaluate ----
    print("\n--- Step 3: Training dual routers ---")

    queries = combined["query_id"].unique()
    gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=RANDOM_SEED)
    train_idx, test_idx = next(gss.split(combined, groups=combined["query_id"]))
    train_df = combined.iloc[train_idx].copy()
    test_df = combined.iloc[test_idx].copy()

    print(f"Train: {train_df['query_id'].nunique()} queries, "
          f"Test: {test_df['query_id'].nunique()} queries")

    # Router-B: predict best_overlap F1
    X_train = train_df[FEATURE_COLS].fillna(0.0)
    model_b = XGBRegressor(
        n_estimators=200, max_depth=6, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        random_state=RANDOM_SEED, n_jobs=-1)
    model_b.fit(X_train, train_df["label_best_overlap_f1"].fillna(0.0))

    # Router-P: predict primary F1
    model_p = XGBRegressor(
        n_estimators=200, max_depth=6, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        random_state=RANDOM_SEED, n_jobs=-1)
    model_p.fit(X_train, train_df["label_primary_f1"].fillna(0.0))

    # Evaluate
    rb_best = evaluate_router(model_b, test_df, "label_best_overlap_f1",
                              "label_best_overlap_f1")
    rp_prim = evaluate_router(model_p, test_df, "label_primary_f1",
                              "label_primary_f1")
    rp_best = evaluate_router(model_p, test_df, "label_primary_f1",
                              "label_best_overlap_f1")

    base_primary = evaluate_baseline(test_df, "always_primary")
    base_oracle = evaluate_baseline(test_df, "oracle")

    # Feature importance
    imp_b = sorted(zip(FEATURE_COLS, model_b.feature_importances_),
                   key=lambda x: -x[1])
    print(f"\nRouter-B (best_overlap): F1={rb_best['best_overlap_mean']:.4f} "
          f"[{rb_best['best_overlap_ci_lower']:.4f}, {rb_best['best_overlap_ci_upper']:.4f}]")
    print(f"Top5 features: {[f'{n}={v:.3f}' for n,v in imp_b[:5]]}")

    print(f"Router-P (primary):       F1={rp_prim['primary_mean']:.4f} "
          f"[{rp_prim['primary_ci_lower']:.4f}, {rp_prim['primary_ci_upper']:.4f}]")
    print(f"  Best_overlap via P:     {rp_best['best_overlap_mean']:.4f}")

    print(f"Baseline always_primary:  best={base_primary['best_overlap_mean']:.4f}, "
          f"prim={base_primary['primary_mean']:.4f}")
    print(f"Oracle:                   best={base_oracle['best_overlap_mean']:.4f}, "
          f"prim={base_oracle['primary_mean']:.4f}")

    # Save
    summary = {
        "species": sp,
        "router_b": rb_best,
        "router_p_primary": rp_prim,
        "router_p_best_overlap": rp_best,
        "baseline_always_primary": base_primary,
        "baseline_oracle": base_oracle,
    }
    with open(sp_outdir / "evaluation_summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)

    # Paper Table II format line
    print(f"\nPaper Table II: {sp}")
    print(f"TPRS Phase 5: primary={rp_prim['primary_mean']:.4f} "
          f"({rp_prim['primary_ci_lower']:.3f},{rp_prim['primary_ci_upper']:.3f}), "
          f"best_overlap={rb_best['best_overlap_mean']:.4f} "
          f"({rb_best['best_overlap_ci_lower']:.3f},{rb_best['best_overlap_ci_upper']:.3f})")

    return summary


def main():
    if len(sys.argv) < 2:
        print("Usage: python scripts/tprs/phase5_run.py <species> [species2 ...]")
        print("  species: mouse, dog, chicken, zebrafish")
        return 1

    species_list = sys.argv[1:]
    all_results = {}

    for sp in species_list:
        if sp not in SPECIES_ORTHOLOG:
            print(f"Unknown species: {sp}, skipping")
            continue
        result = run_species(sp)
        all_results[sp] = result

    # Cross-species summary
    if len(all_results) > 1:
        print(f"\n{'='*70}")
        print("CROSS-SPECIES SUMMARY (Paper Table II format)")
        print(f"{'='*70}")
        print(f"{'Species':<12s} {'Primary F1':>25s} {'Best_overlap F1':>25s}")
        print(f"{'-'*62}")
        for sp, r in all_results.items():
            rp = r["router_p_primary"]
            rb = r["router_b"]
            print(f"{sp:<12s} {rp['primary_mean']:.4f} ({rp['primary_ci_lower']:.3f},{rp['primary_ci_upper']:.3f})  "
                  f"{rb['best_overlap_mean']:.4f} ({rb['best_overlap_ci_lower']:.3f},{rb['best_overlap_ci_upper']:.3f})")

    return 0


if __name__ == "__main__":
    sys.exit(main())
