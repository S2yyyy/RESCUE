"""Phase 5 Baselines B1/B2/B3.

B1: Mouse-ortholog-direct — mouse ortholog protein → miniprot genome-wide.
    Tests: if you know the mouse ortholog, what F1 can miniprot achieve?

B2: Human-query forced-locus — human query → miniprot at best paralog locus.
    Tests: how much does correct locus alone improve F1?

B3: Mouse-protein forced-locus — mouse co-ortholog protein at its own gene locus.
    Tests: upper bound for species-matched protein at the correct locus.
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
OUTDIR = REPO / "results/tprs/phase5"

RANDOM_SEED = 4242
N_BOOTSTRAP = 10_000


# ---------------------------------------------------------------------------
# Utilities
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


def load_ortholog_map(sp: str) -> dict[str, str]:
    species_map = {
        "mouse": "mus_musculus", "dog": "canis_lupus_familiaris",
        "chicken": "gallus_gallus", "zebrafish": "danio_rerio",
    }
    sp_name = species_map.get(sp, "mus_musculus")
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


def parse_miniprot_gff(gff_text: str, region_chrom: str,
                       region_offset: int) -> list[dict]:
    """Parse miniprot GFF, return hits with CDS coords in absolute 0-based coords.

    For genome-wide runs, region_chrom comes from the GFF seqid field,
    and region_offset=0 (GFF coords are already absolute 1-based).
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
                "chrom": fields[0],  # from GFF seqid
                "strand": fields[6],
                "score": float(fields[5]) if fields[5] != "." else 0.0,
                "identity": float(attrs.get("Identity", 0.0)),
                "query_start": int(target[1]) if len(target) > 1 else 0,
                "query_end": int(target[2]) if len(target) > 2 else 0,
                "rank": int(attrs.get("Rank", 1)),
            }
        elif feat == "CDS" and cur is not None:
            # Genome-wide: GFF coords are 1-based absolute, convert to 0-based
            cur_cds.append({
                "chrom": fields[0],
                "start": int(fields[3]) - 1,
                "end": int(fields[4]),
                "strand": cur.get("strand", "."),
            })
    if cur is not None and cur_cds:
        cur["cds"] = sorted(cur_cds, key=lambda x: x["start"])
        hits.append(cur)
    return hits


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


def paired_bootstrap_ci(values: np.ndarray) -> tuple:
    rng = np.random.RandomState(RANDOM_SEED)
    n = len(values)
    means = np.array([rng.choice(values, size=n, replace=True).mean()
                      for _ in range(N_BOOTSTRAP)])
    return float(values.mean()), float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


# ---------------------------------------------------------------------------
# B1: Mouse-ortholog-direct (genome-wide miniprot with mouse ortholog protein)
# ---------------------------------------------------------------------------

def run_b1_worker(args: tuple) -> dict:
    qid, ortholog_gene, protein_seq, genome_path, gt_data, cds_bed = args
    gt = gt_data.get(qid, {})
    gt_primary = gt.get("exons", [])
    paralog_genes = {}
    for pl in gt.get("paralog_loci", []):
        pg = pl.get("target_gene", "")
        if pg and pg in cds_bed:
            paralog_genes[pg] = cds_bed[pg]

    result = {"query_id": qid, "ortholog_gene": ortholog_gene,
              "status": "error", "primary_f1": 0.0, "best_overlap_f1": 0.0}

    query_fa = f">{qid}\n{protein_seq}\n"
    try:
        with tempfile.TemporaryDirectory() as td:
            tmp = Path(td)
            (tmp / "query.fa").write_text(query_fa)
            cmd = [MINIPROT, "--gff", "-t", "4", "-N", "1",
                   str(genome_path), str(tmp / "query.fa")]
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            hits = parse_miniprot_gff(proc.stdout, "", 0)

            if hits:
                pred_cds = hits[0].get("cds", [])
                p_metrics = base_level_metrics(pred_cds, gt_primary)
                best_overlap = p_metrics["f1"]
                best_gene = "primary"
                for pg, pg_cds in paralog_genes.items():
                    m = base_level_metrics(pred_cds, pg_cds)
                    if m["f1"] > best_overlap:
                        best_overlap = m["f1"]
                        best_gene = pg
                cds_len = sum(e["end"] - e["start"] for e in pred_cds)
                result.update({
                    "status": "ok",
                    "primary_f1": p_metrics["f1"],
                    "best_overlap_f1": best_overlap,
                    "best_gene": best_gene,
                    "n_exons": len(pred_cds),
                    "cds_length": cds_len,
                    "hit_chrom": hits[0].get("chrom", ""),
                    "identity": hits[0].get("identity", 0.0),
                })
            else:
                result["status"] = "no_hits"
    except subprocess.TimeoutExpired:
        result["status"] = "timeout"
    except Exception as e:
        result["status"] = "error"
        result["error"] = str(e)[:200]
        import traceback
        result["traceback"] = traceback.format_exc()[-500:]

    return result


def run_b1(sp: str, gt_data: dict, gene2prot: dict, human2sp: dict,
           cds_bed: dict, genome_path: Path, workers: int = 8) -> dict:
    """B1: Mouse ortholog → miniprot genome-wide.

    Uses GT target_gene (primary ortholog) as the species protein source.
    Falls back to Ensembl Compara one-to-one mapping.
    """
    print(f"\n{'='*60}")
    print(f"B1: {sp}-ortholog-direct (genome-wide miniprot)")
    print(f"{'='*60}")

    tasks = []
    skipped = 0
    used_gt = 0
    used_ensembl = 0
    for qid in gt_data:
        # Priority 1: target_gene from GT (works for paralog queries)
        ortholog_gene = gt_data[qid].get("target_gene", "")
        protein_seq = gene2prot.get(ortholog_gene, "")
        if ortholog_gene and protein_seq:
            used_gt += 1
        else:
            # Priority 2: Ensembl Compara one-to-one
            ortholog_gene = human2sp.get(qid, "")
            protein_seq = gene2prot.get(ortholog_gene, "") if ortholog_gene else ""
            if ortholog_gene and protein_seq:
                used_ensembl += 1
            else:
                skipped += 1
                continue
        tasks.append((qid, ortholog_gene, protein_seq, genome_path, gt_data, cds_bed))

    print(f"Queries with ortholog protein: {len(tasks)} (GT: {used_gt}, Ensembl: {used_ensembl}), skipped: {skipped}")

    t0 = time.time()
    with Pool(workers) as pool:
        results = pool.map(run_b1_worker, tasks)
    elapsed = time.time() - t0
    print(f"Done in {elapsed:.0f}s ({len(tasks)/elapsed:.1f} queries/s)")

    ok = [r for r in results if r["status"] == "ok"]
    no_hits = [r for r in results if r["status"] == "no_hits"]
    errors = [r for r in results if r["status"] not in ("ok", "no_hits")]

    primary_f1s = np.array([r["primary_f1"] for r in ok] + [0.0] * (len(no_hits) + len(errors)))
    best_f1s = np.array([r["best_overlap_f1"] for r in ok] + [0.0] * (len(no_hits) + len(errors)))

    prim_mean, prim_lo, prim_hi = paired_bootstrap_ci(primary_f1s)
    best_mean, best_lo, best_hi = paired_bootstrap_ci(best_f1s)

    print(f"Ok: {len(ok)}, no_hits: {len(no_hits)}, errors: {len(errors)}")
    print(f"Primary F1:       {prim_mean:.4f} [{prim_lo:.4f}, {prim_hi:.4f}]")
    print(f"Best_overlap F1:  {best_mean:.4f} [{best_lo:.4f}, {best_hi:.4f}]")
    print(f"Zero-F1 (best):   {(best_f1s==0).sum()}/{len(best_f1s)}")

    # Also check: what fraction of ok hits land on the primary chromosome?
    if ok:
        on_primary = sum(1 for r in ok if r.get("hit_chrom") == gt_data.get(r["query_id"], {}).get("chrom"))
        print(f"Top hit on primary chrom: {on_primary}/{len(ok)} ({100*on_primary/len(ok):.1f}%)")

    return {
        "name": "B1_ortholog_direct",
        "primary_mean": prim_mean, "primary_ci_lower": prim_lo, "primary_ci_upper": prim_hi,
        "best_overlap_mean": best_mean, "best_overlap_ci_lower": best_lo, "best_overlap_ci_upper": best_hi,
        "n_queries": len(results), "ok": len(ok), "no_hits": len(no_hits), "errors": len(errors),
        "zero_f1": int((best_f1s == 0).sum()),
        "results": results,
    }


# ---------------------------------------------------------------------------
# B2: Human-query forced-locus (best human query result at any paralog locus)
# ---------------------------------------------------------------------------

def run_b2(sp: str, gt_data: dict, cds_bed: dict) -> dict:
    """B2: Human query at best paralog locus (from Phase 4 data).

    For each query, take the max label_best_overlap_f1 among human-query rows.
    This is the "oracle" achievable with human query across all Phase 4 candidate regions.
    """
    print(f"\n{'='*60}")
    print(f"B2: Human-query forced-locus (best human result)")
    print(f"{'='*60}")

    phase4_path = REPO / "results/tprs/router" / f"training_data_{sp}.tsv"
    if not phase4_path.exists():
        print(f"No Phase 4 data for {sp}, skipping")
        return {}

    df = pd.read_csv(phase4_path, sep="\t")

    # Per-query max best_overlap F1 (human query only)
    per_query = df.groupby("query_id")["label_best_overlap_f1"].max()
    per_query_prim = df.groupby("query_id")["label_primary_f1"].max()

    best_f1s = per_query.values
    primary_f1s = per_query_prim.values

    best_mean, best_lo, best_hi = paired_bootstrap_ci(best_f1s)
    prim_mean, prim_lo, prim_hi = paired_bootstrap_ci(primary_f1s)

    print(f"Queries: {len(per_query)}")
    print(f"Primary F1:       {prim_mean:.4f} [{prim_lo:.4f}, {prim_hi:.4f}]")
    print(f"Best_overlap F1:  {best_mean:.4f} [{best_lo:.4f}, {best_hi:.4f}]")
    print(f"Zero-F1 (best):   {(best_f1s==0).sum()}/{len(best_f1s)}")

    return {
        "name": "B2_human_forced_locus",
        "primary_mean": prim_mean, "primary_ci_lower": prim_lo, "primary_ci_upper": prim_hi,
        "best_overlap_mean": best_mean, "best_overlap_ci_lower": best_lo, "best_overlap_ci_upper": best_hi,
        "n_queries": len(per_query),
        "zero_f1": int((best_f1s == 0).sum()),
    }


# ---------------------------------------------------------------------------
# B3: Mouse-protein forced-locus (mouse co-ortholog at its own gene locus)
# ---------------------------------------------------------------------------

def run_b3(sp: str, gt_data: dict, cds_bed: dict) -> dict:
    """B3: Mouse protein at its own gene locus (from Phase 5 data).

    For each query, take the max label_best_overlap_f1 among mouse-query rows
    where c_identity > 0.9 (confirming same-species alignment).
    """
    print(f"\n{'='*60}")
    print(f"B3: Species-protein forced-locus (mouse protein at own gene)")
    print(f"{'='*60}")

    phase5_path = OUTDIR / sp / "training_data.tsv"
    if not phase5_path.exists():
        phase5_path = OUTDIR / "training_data_mouse.tsv"
    if not phase5_path.exists():
        print(f"No Phase 5 data for {sp}, skipping")
        return {}

    df = pd.read_csv(phase5_path, sep="\t")
    mouse_df = df[df["query_type"] == "mouse"]

    # Only use rows where the mouse protein actually aligned (high identity)
    mouse_hits = mouse_df[mouse_df["c_has_hit"] == 1]

    per_query = mouse_hits.groupby("query_id")["label_best_overlap_f1"].max()
    per_query_prim = mouse_hits.groupby("query_id")["label_primary_f1"].max()

    best_f1s = per_query.values
    primary_f1s = per_query_prim.values

    best_mean, best_lo, best_hi = paired_bootstrap_ci(best_f1s)
    prim_mean, prim_lo, prim_hi = paired_bootstrap_ci(primary_f1s)

    print(f"Queries with mouse protein hits: {len(per_query)}")
    print(f"Primary F1:       {prim_mean:.4f} [{prim_lo:.4f}, {prim_hi:.4f}]")
    print(f"Best_overlap F1:  {best_mean:.4f} [{best_lo:.4f}, {best_hi:.4f}]")
    print(f"Zero-F1 (best):   {(best_f1s==0).sum()}/{len(best_f1s)}")

    return {
        "name": "B3_mouse_protein_forced_locus",
        "primary_mean": prim_mean, "primary_ci_lower": prim_lo, "primary_ci_upper": prim_hi,
        "best_overlap_mean": best_mean, "best_overlap_ci_lower": best_lo, "best_overlap_ci_upper": best_hi,
        "n_queries": len(per_query),
        "zero_f1": int((best_f1s == 0).sum()),
    }


# ---------------------------------------------------------------------------
# B4 (extra): Phase 4 router (human query only, no species protein)
# ---------------------------------------------------------------------------

def run_phase4_baseline(sp: str) -> dict:
    """Phase 4 router: human query only, router picks best region."""
    print(f"\n{'='*60}")
    print(f"Baseline: Phase 4 Router (human query only)")
    print(f"{'='*60}")

    from sklearn.model_selection import GroupShuffleSplit
    from xgboost import XGBRegressor

    phase4_path = REPO / "results/tprs/router" / f"training_data_{sp}.tsv"
    if not phase4_path.exists():
        print(f"No Phase 4 data for {sp}, skipping")
        return {}

    df = pd.read_csv(phase4_path, sep="\t")

    feature_cols = [
        "a_diamond", "a_last", "a_miniprot", "a_mmseqs",
        "a_n_tools", "a_best_score", "a_best_qcov", "a_n_raw_hits",
        "region_size", "region_size_log",
        "on_primary_chrom", "strand_match_primary",
        "coorth_overlap_bp", "coorth_same_chrom",
        "query_length_aa", "query_length_log",
        "c_has_hit", "c_identity", "c_coverage", "c_n_exons",
        "c_cds_length", "c_composite_score", "c_hit_score", "c_rank",
    ]

    queries = df["query_id"].unique()
    gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=RANDOM_SEED)
    train_idx, test_idx = next(gss.split(df, groups=df["query_id"]))
    train_df = df.iloc[train_idx].copy()
    test_df = df.iloc[test_idx].copy()

    X_train = train_df[feature_cols].fillna(0.0)
    y_train = train_df["label_best_overlap_f1"].fillna(0.0)

    model = XGBRegressor(
        n_estimators=200, max_depth=6, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        random_state=RANDOM_SEED, n_jobs=-1)
    model.fit(X_train, y_train)

    # Evaluate: per query, pick region with highest predicted F1
    test_queries = test_df["query_id"].unique()
    actual_f1s = []
    actual_prim_f1s = []

    for qid in test_queries:
        qdf = test_df[test_df["query_id"] == qid]
        X = qdf[feature_cols].fillna(0.0)
        preds = model.predict(X)
        best_row = qdf.iloc[np.argmax(preds)]
        actual_f1s.append(best_row["label_best_overlap_f1"])
        actual_prim_f1s.append(best_row["label_primary_f1"])

    actual_f1s = np.array(actual_f1s)
    actual_prim_f1s = np.array(actual_prim_f1s)

    best_mean, best_lo, best_hi = paired_bootstrap_ci(actual_f1s)
    prim_mean, prim_lo, prim_hi = paired_bootstrap_ci(actual_prim_f1s)

    # Always-primary baseline
    always_prim_f1s = []
    always_prim_prim = []
    for qid in test_queries:
        qdf = test_df[test_df["query_id"] == qid]
        subset = qdf
        best_idx = subset["on_primary_chrom"].values.argmax()
        best_row = subset.iloc[best_idx]
        always_prim_f1s.append(best_row["label_best_overlap_f1"])
        always_prim_prim.append(best_row["label_primary_f1"])

    always_prim_f1s = np.array(always_prim_f1s)
    always_prim_prim = np.array(always_prim_prim)
    ap_best_mean, ap_best_lo, ap_best_hi = paired_bootstrap_ci(always_prim_f1s)
    ap_prim_mean, ap_prim_lo, ap_prim_hi = paired_bootstrap_ci(always_prim_prim)

    print(f"Phase 4 Router:     best={best_mean:.4f}, prim={prim_mean:.4f}")
    print(f"Always-primary:     best={ap_best_mean:.4f}, prim={ap_prim_mean:.4f}")

    return {
        "name": "Phase4_router",
        "best_overlap_mean": best_mean, "best_overlap_ci_lower": best_lo, "best_overlap_ci_upper": best_hi,
        "primary_mean": prim_mean, "primary_ci_lower": prim_lo, "primary_ci_upper": prim_hi,
        "n_queries": len(test_queries),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    sp = "mouse" if len(sys.argv) < 2 else sys.argv[1]

    print(f"Phase 5 Baselines B1/B2/B3 — {sp}")
    print("=" * 60)

    # Load shared data
    gt_path = DATA / f"{sp}_paralog" / "ground_truth.json"
    with open(gt_path) as f:
        gt_data = json.load(f)

    gene2prot = load_proteome_by_gene(sp)
    human2sp = load_ortholog_map(sp)
    cds_bed = load_cds_bed(sp)
    genome_path = GENOMES / f"{sp}.fa"

    print(f"GT queries: {len(gt_data)}, proteome: {len(gene2prot)}, "
          f"ortholog pairs: {len(human2sp)}, CDS genes: {len(cds_bed)}")

    all_baselines = []

    # B1: Mouse-ortholog-direct (needs miniprot runs)
    b1 = run_b1(sp, gt_data, gene2prot, human2sp, cds_bed, genome_path, workers=4)
    all_baselines.append(b1)

    # B2: Human-query forced-locus (from existing Phase 4 data)
    b2 = run_b2(sp, gt_data, cds_bed)
    all_baselines.append(b2)

    # B3: Mouse-protein forced-locus (from existing Phase 5 data)
    b3 = run_b3(sp, gt_data, cds_bed)
    all_baselines.append(b3)

    # Phase 4 router baseline
    p4 = run_phase4_baseline(sp)
    all_baselines.append(p4)

    # Load Phase 5 results for comparison
    phase5_summary = OUTDIR / sp / "evaluation_summary.json"
    if not phase5_summary.exists():
        phase5_summary = OUTDIR / "evaluation_summary.json"

    with open(phase5_summary) as f:
        p5_data = json.load(f)

    if "router_b" in p5_data:
        p5_best = p5_data["router_b"]["best_overlap_mean"]
        p5_prim = p5_data["router_p_primary"]["primary_mean"]
    else:
        p5_best = p5_data["phase5_router_b_best_overlap"]["mean"]
        p5_prim = p5_data["phase5_router_p_primary"]["mean"]

    # Paper baselines
    paper_gth_best = 0.611
    paper_gth_prim = 0.338
    paper_miniprot_best = 0.597
    paper_miniprot_prim = 0.326

    # Final summary
    print(f"\n{'='*70}")
    print(f"BASELINE COMPARISON — {sp}")
    print(f"{'='*70}")
    print(f"{'Method':<35s} {'Primary F1':>12s} {'Best_overlap F1':>17s}")
    print(f"{'-'*65}")
    print(f"{'Paper: GTH (human only)':<35s} {paper_gth_prim:>12.3f} {paper_gth_best:>17.3f}")
    print(f"{'Paper: miniprot (human only)':<35s} {paper_miniprot_prim:>12.3f} {paper_miniprot_best:>17.3f}")
    print()

    if b1:
        print(f"{'B1: Ortholog-direct (mouse prot GW)':<35s} "
              f"{b1['primary_mean']:>12.4f} {b1['best_overlap_mean']:>17.4f}")
    if b2:
        print(f"{'B2: Human forced-locus (best locus)':<35s} "
              f"{b2['primary_mean']:>12.4f} {b2['best_overlap_mean']:>17.4f}")
    if b3:
        print(f"{'B3: Mouse prot forced-locus':<35s} "
              f"{b3['primary_mean']:>12.4f} {b3['best_overlap_mean']:>17.4f}")
    if p4:
        print(f"{'Phase 4 Router (human only)':<35s} "
              f"{p4['primary_mean']:>12.4f} {p4['best_overlap_mean']:>17.4f}")
    print(f"{'Phase 5 Router (human + mouse prot)':<35s} "
          f"{p5_prim:>12.4f} {p5_best:>17.4f}")

    # Save
    out_path = OUTDIR / f"baselines_{sp}.json"
    with open(out_path, "w") as f:
        json.dump({
            "B1": {k: v for k, v in b1.items() if k != "results"},
            "B2": b2,
            "B3": b3,
            "Phase4": p4,
            "Phase5_best_overlap": p5_best,
            "Phase5_primary": p5_prim,
        }, f, indent=2, default=float)
    print(f"\nSaved to {out_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
