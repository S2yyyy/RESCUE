"""TPRS Phase 2+3 — C-stage filtering + co-ortholog guided recovery.

Phase 2: Run miniprot C-stage on merged candidate regions from all A-tools,
score each region by alignment quality, output scored candidates.

Phase 3: For blind-miss proteins where a co-ortholog candidate produces a
good C-stage alignment, expand the search region (±2Mb) and re-run miniprot.
The primary ortholog is often on the same chromosome as the co-ortholog.

Outputs:
  results/tprs/cstage/<protein_id>.<region_idx>.json  — per-candidate scores
  results/tprs/cstage/manifest.tsv                       — all scores summary
  results/tprs/recovery/<protein_id>.gff                — Phase 3 expanded hits
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
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import GENOMES, MINIPROT, PROJECT_ROOT, V3_DATA  # noqa: E402

REPO = PROJECT_ROOT
BED_DIR = REPO / "results/tprs/stage0_bed"
EXP3_BED = REPO / "results/exp3/stage0_bed"
GT_DIR = REPO / "results/exp2/gt"
DATA = V3_DATA
RAW = REPO / "results/tprs/raw"
CSTAGE = REPO / "results/tprs/cstage"
RECOVERY = REPO / "results/tprs/recovery"

A_TOOLS = [
    "miniprot_relaxed", "miniprot", "genblasta", "mmseqs", "last",
    "diamond", "spaln3", "genblastg", "tblastn"
]
SPECIES = ["mouse", "dog", "chicken", "zebrafish"]

# Miniprot binary
MINIPROT_BIN = MINIPROT


def load_bed(path: Path) -> dict[str, list[tuple]]:
    out: dict[str, list[tuple]] = defaultdict(list)
    if not path.exists():
        return out
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.strip().split("\t")
            # chrom, start, end, name, score, strand, q_start, q_end, q_cov
            out[f[3]].append((f[0], int(f[1]), int(f[2]), f[5],
                             float(f[4]) if len(f) > 4 else 0.0))
    return out


def load_gt(path: Path) -> dict[str, list[tuple]]:
    out: dict[str, list[tuple]] = defaultdict(list)
    with path.open() as fh:
        header = fh.readline().strip().split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.strip().split("\t")
            out[f[idx["query_id"]]].append((
                f[idx["chrom"]], int(f[idx["start"]]),
                int(f[idx["end"]]), f[idx["strand"]],
            ))
    return out


def recip_overlap_50(a, b) -> bool:
    if a[0] != b[0] or a[3] != b[3]:
        return False
    inter = max(0, min(a[2], b[2]) - max(a[1], b[1]))
    if inter == 0:
        return False
    return inter / (a[2] - a[1]) >= 0.5 and inter / (b[2] - b[1]) >= 0.5


def merge_regions(candidates: list[tuple], slop: int = 0) -> list[tuple]:
    """Merge overlapping candidate regions per chromosome+strand."""
    if not candidates:
        return []
    by_key: dict[tuple, list[tuple]] = defaultdict(list)
    for c in candidates:
        by_key[(c[0], c[3])].append(c)
    merged = []
    for (chrom, strand), regs in by_key.items():
        regs_sorted = sorted(regs, key=lambda x: x[1])
        cur_chrom, cur_start, cur_end, cur_strand = regs_sorted[0][:4]
        for r in regs_sorted[1:]:
            if r[1] <= cur_end + slop:
                cur_end = max(cur_end, r[2])
            else:
                merged.append((cur_chrom, cur_start - slop, cur_end + slop, cur_strand))
                cur_chrom, cur_start, cur_end, cur_strand = r[:4]
        merged.append((cur_chrom, cur_start - slop, cur_end + slop, cur_strand))
    return [(c[0], max(0, c[1]), c[2], c[3]) for c in merged]


def extract_region_fasta(genome: Path, chrom: str, start: int, end: int,
                          out_fa: Path) -> str:
    """Extract genomic region via samtools faidx."""
    region = f"{chrom}:{start}-{end}"
    cmd = ["samtools", "faidx", str(genome), region]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"samtools faidx failed: {proc.stderr[:200]}")
    with out_fa.open("w") as fh:
        fh.write(proc.stdout)
    return proc.stdout


def parse_miniprot_gff(gff_text: str) -> list[dict]:
    """Parse miniprot GFF output, return list of transcript models."""
    transcripts = []
    cur = None
    cur_cds = []
    for line in gff_text.splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        feat = fields[2]
        if feat == "mRNA":
            if cur is not None and cur_cds:
                cur["cds"] = cur_cds
                transcripts.append(cur)
            cur_cds = []
            attrs = dict(p.split("=", 1) for p in fields[8].split(";") if "=" in p)
            target = attrs.get("Target", "").split()
            cur = {
                "strand": fields[6],
                "score": float(fields[5]) if fields[5] != "." else 0.0,
                "identity": float(attrs.get("Identity", 0.0)),
                "positive": float(attrs.get("Positive", 0.0)),
                "query_start": int(target[1]) if len(target) > 1 else 0,
                "query_end": int(target[2]) if len(target) > 2 else 0,
                "rank": int(attrs.get("Rank", 1)),
            }
        elif feat == "CDS" and cur is not None:
            cur_cds.append({
                "start": int(fields[3]),
                "end": int(fields[4]),
            })
    if cur is not None and cur_cds:
        cur["cds"] = cur_cds
        transcripts.append(cur)
    return transcripts


def score_alignment(hit: dict, query_len: int) -> dict:
    """Score a miniprot alignment by identity, coverage, and CDS completeness."""
    cds = hit.get("cds", [])
    if not cds:
        return {"score": 0.0, "identity": 0.0, "coverage": 0.0,
                "n_exons": 0, "cds_length": 0, "label": "no_cds"}

    identity = hit.get("identity", 0.0)
    cds_length = sum(e["end"] - e["start"] for e in cds)
    coverage = cds_length / (query_len * 3) if query_len > 0 else 0.0
    n_exons = len(cds)

    # Composite score: identity * coverage * exon_bonus
    exon_bonus = min(1.0, n_exons / 3.0)  # penalty for single-exon hits
    composite = identity * min(1.0, coverage) * exon_bonus

    label = "OK"
    if identity < 0.3:
        label = "low_identity"
    elif coverage < 0.3:
        label = "low_coverage"
    elif n_exons == 1 and coverage < 0.8:
        label = "single_exon_partial"

    return {
        "score": round(composite, 4),
        "identity": round(identity, 4),
        "coverage": round(coverage, 4),
        "n_exons": n_exons,
        "cds_length": cds_length,
        "label": label,
    }


def run_miniprot_on_region(genome: Path, query_fa: Path, region: tuple,
                           timeout: int = 60) -> tuple[str, float]:
    """Run miniprot on a genomic region, return GFF output and wall time."""
    chrom, start, end, strand = region
    with tempfile.TemporaryDirectory() as td:
        scratch = Path(td)
        region_fa = scratch / "region.fa"
        extract_region_fasta(genome, chrom, start, end, region_fa)
        cmd = [MINIPROT_BIN, "--gff", "-t", "1", "-N", "5",
               str(region_fa), str(query_fa)]
        t0 = time.time()
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        wall = time.time() - t0
        return proc.stdout, wall


def run_phase2_cstage(sp: str, query_ids: list[str] | None = None,
                      max_regions_per_query: int = 20,
                      timeout: int = 60) -> dict:
    """Phase 2: C-stage filter on merged candidate regions."""
    log = logging.getLogger("phase2")
    genome = GENOMES / f"{sp}.fa"
    query_fa = DATA / f"{sp}_paralog" / "proteins.fa"
    gt_broad = GT_DIR / f"{sp}.paralog_broad.tsv"

    # Load GT
    gts = load_gt(gt_broad)
    all_qids = sorted(set(gts))
    if query_ids:
        all_qids = [q for q in all_qids if q in query_ids]

    # Load query lengths
    qlen = {}
    with open(query_fa) as fh:
        cur_id = None
        for line in fh:
            if line.startswith(">"):
                cur_id = line[1:].strip().split()[0]
                qlen[cur_id] = 0
            elif cur_id:
                qlen[cur_id] += len(line.strip())

    # Load all tool candidates and merge
    all_cands: dict[str, list[tuple]] = defaultdict(list)
    for tool in A_TOOLS:
        if tool == "miniprot_relaxed":
            bed = BED_DIR / f"{tool}.{sp}.paralog.bed"
        else:
            bed = EXP3_BED / f"{tool}.{sp}.paralog.bed"
        if not bed.exists():
            continue
        for qid, cands in load_bed(bed).items():
            all_cands[qid].extend(cands)

    CSTAGE.mkdir(parents=True, exist_ok=True)
    manifest_rows = []

    for qid in all_qids:
        if qid not in qlen:
            continue
        ql = qlen[qid]

        # Merge overlapping candidates
        raw_cands = all_cands.get(qid, [])
        merged = merge_regions(raw_cands, slop=50000)
        # Sort by some priority (number of tools supporting, then size)
        # For now, just take top max_regions_per_query
        merged = merged[:max_regions_per_query]

        # Write query protein FASTA
        with tempfile.TemporaryDirectory() as td:
            scratch = Path(td)
            q_fa = scratch / "query.fa"
            # Extract single query from multi-fasta
            with open(query_fa) as fin, open(q_fa, "w") as fout:
                found = False
                for line in fin:
                    if line.startswith(">"):
                        found = (line[1:].strip().split()[0] == qid)
                    if found:
                        fout.write(line)

            for i, region in enumerate(merged):
                chrom, start, end, strand = region
                out_json = CSTAGE / f"{qid}.{i}.json"
                if out_json.exists():
                    continue

                try:
                    gff_out, wall = run_miniprot_on_region(
                        genome, q_fa, region, timeout=timeout)
                except subprocess.TimeoutExpired:
                    manifest_rows.append({
                        "query_id": qid, "region_idx": i, "chrom": chrom,
                        "start": start, "end": end, "strand": strand,
                        "status": "timeout", "wall_s": timeout,
                        "score": 0.0, "identity": 0.0, "coverage": 0.0,
                        "n_exons": 0, "cds_length": 0, "label": "timeout",
                    })
                    continue
                except Exception as e:
                    manifest_rows.append({
                        "query_id": qid, "region_idx": i, "chrom": chrom,
                        "start": start, "end": end, "strand": strand,
                        "status": "error", "wall_s": 0.0,
                        "score": 0.0, "identity": 0.0, "coverage": 0.0,
                        "n_exons": 0, "cds_length": 0,
                        "label": f"error:{type(e).__name__}",
                    })
                    continue

                hits = parse_miniprot_gff(gff_out)
                if not hits:
                    manifest_rows.append({
                        "query_id": qid, "region_idx": i, "chrom": chrom,
                        "start": start, "end": end, "strand": strand,
                        "status": "no_hits", "wall_s": round(wall, 1),
                        "score": 0.0, "identity": 0.0, "coverage": 0.0,
                        "n_exons": 0, "cds_length": 0, "label": "no_hits",
                    })
                    continue

                # Take best hit (rank 1)
                best = hits[0]
                scores = score_alignment(best, ql)

                rec = {
                    "query_id": qid, "region_idx": i, "chrom": chrom,
                    "start": start, "end": end, "strand": strand,
                    "status": "ok", "wall_s": round(wall, 1),
                    **scores,
                    "cds_json": json.dumps(best.get("cds", [])),
                }
                manifest_rows.append(rec)

                out_json.parent.mkdir(parents=True, exist_ok=True)
                with open(out_json, "w") as fh:
                    json.dump(rec, fh, indent=2)

        if len(all_qids) > 0 and all_qids.index(qid) % 50 == 0:
            log.info("phase2 %s: %d/%d queries done", sp, all_qids.index(qid), len(all_qids))

    # Write manifest
    manifest_path = CSTAGE / f"manifest_{sp}.tsv"
    cols = ["query_id", "region_idx", "chrom", "start", "end", "strand",
            "status", "wall_s", "score", "identity", "coverage",
            "n_exons", "cds_length", "label"]
    with open(manifest_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in manifest_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    log.info("phase2 %s: wrote %d rows to %s", sp, len(manifest_rows), manifest_path)
    return {"species": sp, "n_queries": len(all_qids), "n_regions": len(manifest_rows)}


def main() -> int:
    ap = argparse.ArgumentParser(description="TPRS Phase 2+3: C-stage filter + family recovery")
    ap.add_argument("--species", nargs="*", default=SPECIES)
    ap.add_argument("--max-regions", type=int, default=20)
    ap.add_argument("--timeout", type=int, default=60,
                    help="Timeout per miniprot run (seconds)")
    ap.add_argument("--phase3", action="store_true",
                    help="Run Phase 3 expanded recovery on blind-miss proteins")
    args = ap.parse_args()

    CSTAGE.mkdir(parents=True, exist_ok=True)
    log_path = CSTAGE / "phase2.log"
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] TPRS %(message)s",
        handlers=[logging.StreamHandler(sys.stdout),
                  logging.FileHandler(log_path, mode="a")],
    )
    log = logging.getLogger("phase2")

    if not MINIPROT_BIN:
        log.error("miniprot binary not found. Set MINIPROT_BIN env var.")
        return 1

    for sp in args.species:
        log.info("begin phase2 %s", sp)
        result = run_phase2_cstage(sp, max_regions_per_query=args.max_regions,
                                   timeout=args.timeout)
        log.info("done phase2 %s: %d queries, %d regions",
                 sp, result["n_queries"], result["n_regions"])

    if args.phase3:
        log.info("Phase 3 not yet implemented in this script")
        # Phase 3 will be implemented as a separate analysis + recovery step

    return 0


if __name__ == "__main__":
    sys.exit(main())
