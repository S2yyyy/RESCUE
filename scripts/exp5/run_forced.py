"""Exp 5 — Forced-locus driver.

For each (tool, species, regime, padding):
  - Load proteins and ground truth.
  - For each protein, extract its GT-region FASTA (padded).
  - Call the tool's runner on that region.
  - Evaluate the returned hit via v3_evaluate.evaluate_hit.
  - Write a per-protein TSV row.

Usage:
    python -m scripts.exp5.run_forced \
        --tool miniprot --species mouse --regime one2one --padding 50000 \
        --max-proteins 50 --workers 8
"""
from __future__ import annotations

import argparse
import importlib
import logging
import os
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from .common import (REPO, PER_PROTEIN_COLS,
                     load_proteins, load_ground_truth,
                     extract_region_fa, gt_region, chrom_lengths,
                     evaluate_hit_default)


OUT_DIR = REPO / "results/exp5/per_protein"

RUNNERS = {
    "miniprot": "scripts.exp5.runners.miniprot",
    "spaln3":   "scripts.exp5.runners.spaln3",
    "gth":      "scripts.exp5.runners.gth",
    "exonerate":"scripts.exp5.runners.exonerate",
    "genewise": "scripts.exp5.runners.genewise",
    "scipio":   "scripts.exp5.runners.scipio",
    "genblastg":"scripts.exp5.runners.genblastg",
    "minisplice_miniprot": "scripts.exp5.runners.minisplice",
}


def _process_one(args_tuple):
    (tool, species, regime, padding, pid, prot_seq, gt_entry,
     chrom_len_map, options, scratch_root) = args_tuple

    runner_mod = importlib.import_module(RUNNERS[tool])
    chrom, rstart, rend, strand, region_offset = gt_region(
        gt_entry, padding, chrom_len_map)

    with tempfile.TemporaryDirectory(dir=scratch_root, prefix=f"{tool}_{pid}_") as td:
        scratch = Path(td)
        region_fa = scratch / "region.fa"
        try:
            extract_region_fa(species, chrom, rstart, rend, region_fa)
        except Exception as e:
            return {
                "protein_id": pid, "tool": tool, "species": species,
                "regime": regime, "padding": padding,
                "failed": True, "wallclock_s": 0,
                **{k: None for k in PER_PROTEIN_COLS
                   if k not in {"protein_id","tool","species","regime","padding","failed","wallclock_s"}},
            }

        t0 = time.time()
        try:
            hit = runner_mod.run_one(
                pid, prot_seq, region_fa, chrom, rstart, strand,
                scratch, options)
        except Exception as e:
            hit = None
        wallclock = time.time() - t0

    metrics = evaluate_hit_default(hit, gt_entry, len(prot_seq))
    row = {
        "protein_id": pid, "tool": tool, "species": species,
        "regime": regime, "padding": padding,
        "failed": metrics["failed"],
        "wallclock_s": round(wallclock, 3),
    }
    for k in PER_PROTEIN_COLS:
        if k not in row:
            row[k] = metrics.get(k)
    return row


def run_cell(tool: str, species: str, regime: str, padding: int,
             max_proteins: int | None, workers: int,
             query_ids: list[str] | None,
             options: dict, log: logging.Logger,
             scratch_root: str) -> Path:

    proteins = load_proteins(species, "one2one" if regime == "one2one" else "paralog")
    gt = load_ground_truth(species, "one2one" if regime == "one2one" else "paralog")

    if query_ids:
        pids = [q for q in query_ids if q in proteins and q in gt]
    else:
        pids = sorted(set(proteins) & set(gt))
    if max_proteins:
        pids = pids[:max_proteins]
    log.info("%s %s %s pad=%d: %d proteins", tool, species, regime, padding, len(pids))

    chrom_len_map = chrom_lengths(species)

    out_path = OUT_DIR / f"{tool}.{species}.{regime}.{padding}.tsv"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    args_iter = [
        (tool, species, regime, padding, pid, proteins[pid], gt[pid],
         chrom_len_map, options, scratch_root)
        for pid in pids
    ]

    rows: list[dict] = []
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(_process_one, a): a[4] for a in args_iter}
            done = 0
            for fut in as_completed(futs):
                rows.append(fut.result())
                done += 1
                if done % 100 == 0 or done == len(pids):
                    n_failed = sum(1 for r in rows if r["failed"])
                    log.info("  %d/%d done (%d failed) for %s %s %s",
                             done, len(pids), n_failed, tool, species, regime)
    else:
        for i, a in enumerate(args_iter):
            rows.append(_process_one(a))
            if (i + 1) % 50 == 0 or i + 1 == len(pids):
                n_failed = sum(1 for r in rows if r["failed"])
                log.info("  %d/%d done (%d failed)", i + 1, len(pids), n_failed)

    rows.sort(key=lambda r: r["protein_id"])
    with out_path.open("w") as fh:
        fh.write("\t".join(PER_PROTEIN_COLS) + "\n")
        for r in rows:
            fh.write("\t".join(
                "" if r.get(c) is None else str(r.get(c))
                for c in PER_PROTEIN_COLS) + "\n")
    log.info("wrote %s (%d rows; %d failed)", out_path, len(rows),
             sum(1 for r in rows if r["failed"]))
    return out_path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True, choices=sorted(RUNNERS))
    ap.add_argument("--species", required=True,
                    choices=["mouse", "dog", "chicken", "zebrafish"])
    ap.add_argument("--regime", default="one2one",
                    choices=["one2one", "paralog_strict", "paralog_broad"])
    ap.add_argument("--padding", type=int, default=50000)
    ap.add_argument("--max-proteins", type=int, default=None)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--queries", default=None,
                    help="Path to file with one protein_id per line")
    ap.add_argument("--scratch-root", default="/tmp")
    ap.add_argument("--timeout", type=int, default=120)
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO,
                        format="[%(asctime)s] B1 %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])
    log = logging.getLogger("run_forced")

    query_ids = None
    if args.queries:
        query_ids = [l.strip() for l in open(args.queries) if l.strip()]

    options = {"timeout": args.timeout, "species": args.species}
    run_cell(args.tool, args.species, args.regime, args.padding,
             args.max_proteins, args.workers, query_ids, options, log,
             args.scratch_root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
