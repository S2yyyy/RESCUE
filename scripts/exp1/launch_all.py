#!/usr/bin/env python3
"""Exp 1 scheduler: run 7 tools x 4 species x 2 regimes = 56 batches.

LPT-ish ordering: per-query (slow) tools first, then batch tools. Global cap
on concurrent batches from configs/exp1_tools.yaml (default 4). Each batch runs
as `python run_batch.py --tool ... --species ... --regime ...`.

Resumable: batches whose out JSON exists and is non-trivial are skipped.
Per-batch log -> results/exp1/logs/{tool}.{species}.{regime}.log
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from rescue_paths import PROJECT_ROOT  # noqa: E402

ROOT = PROJECT_ROOT
sys.path.insert(0, str(HERE))

import yaml  # noqa: E402
from common import out_path, log_path  # noqa: E402

CFG = yaml.safe_load((ROOT / "configs/exp1_tools.yaml").read_text())
MAX_CONC = CFG["defaults"]["max_concurrent_batches"]
PY = sys.executable
RUNNER = str(HERE / "run_batch.py")

# order: slow (per_query) first to minimize tail latency
_PER_QUERY = [t for t, c in CFG["tools"].items() if c["mode"] == "per_query"]
_BATCH = [t for t, c in CFG["tools"].items() if c["mode"] == "batch"]
TOOL_ORDER = _PER_QUERY + _BATCH


def already_done(tool, species, regime):
    suffix = CFG["regimes"][regime]
    p = out_path(tool, species, suffix)
    return p.exists() and p.stat().st_size > 100


def _run_one(tool, species, regime):
    suffix = CFG["regimes"][regime]
    log = log_path(tool, species, suffix)
    log.parent.mkdir(parents=True, exist_ok=True)
    cmd = [PY, RUNNER, "--tool", tool, "--species", species, "--regime", regime]
    t0 = time.time()
    with log.open("w") as fh:
        fh.write(f"$ {' '.join(cmd)}\n")
        fh.flush()
        rc = subprocess.call(cmd, stdout=fh, stderr=subprocess.STDOUT)
    return (tool, species, regime, rc, time.time() - t0)


def main():
    import argparse as _ap
    ap = _ap.ArgumentParser()
    ap.add_argument("--tools", nargs="*", default=None,
                    help="Only run these tools (default: all)")
    ap.add_argument("--species", nargs="*", default=None)
    ap.add_argument("--regimes", nargs="*", default=None)
    ap.add_argument("--max-workers", type=int, default=None,
                    help="Override max concurrent batches (default: from config)")
    args = ap.parse_args()
    only_tools = set(args.tools) if args.tools else None
    only_sp = set(args.species) if args.species else None
    only_rg = set(args.regimes) if args.regimes else None
    jobs = []
    for tool in TOOL_ORDER:
        if only_tools and tool not in only_tools:
            continue
        for species in CFG["species"]:
            if only_sp and species not in only_sp:
                continue
            for regime in CFG["regimes"]:
                if only_rg and regime not in only_rg:
                    continue
                if already_done(tool, species, regime):
                    print(f"skip (done): {tool}.{species}.{regime}")
                    continue
                jobs.append((tool, species, regime))
    max_w = args.max_workers if args.max_workers else MAX_CONC
    print(f"To run: {len(jobs)} batches, max {max_w} concurrent")
    print(f"Order: {TOOL_ORDER}")
    if not jobs:
        return

    t_all = time.time()
    with ProcessPoolExecutor(max_workers=max_w) as pool:
        futures = {pool.submit(_run_one, *j): j for j in jobs}
        for f in as_completed(futures):
            tool, species, regime, rc, dt = f.result()
            tag = "OK" if rc == 0 else f"FAIL rc={rc}"
            print(f"  [{tag}] {tool}.{species}.{regime}  {dt:.0f}s")
    print(f"TOTAL: {(time.time() - t_all)/3600:.2f} h")


if __name__ == "__main__":
    main()
