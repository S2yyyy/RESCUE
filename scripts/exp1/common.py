"""Shared helpers for Exp 1: data loading, GFF parsing, output format."""
from __future__ import annotations

import json
import os
import sys
import subprocess
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from rescue_paths import (  # noqa: E402
    CONDA,
    GENOMES,
    PROJECT_ROOT,
    RESULTS_ROOT,
    TOOLS_ROOT,
    V3_DATA,
)

ROOT = PROJECT_ROOT
sys.path.insert(0, str(ROOT / "scripts"))

from v3_tool_common import load_proteins, load_ground_truth, load_meta  # noqa: E402
from v3_run_miniprot import parse_miniprot_gff  # noqa: E402

DATA = V3_DATA
OUT = RESULTS_ROOT / "exp1/raw"
LOG = RESULTS_ROOT / "exp1/logs"
WORK = RESULTS_ROOT / "exp1/work"
SCORES = RESULTS_ROOT / "exp1/scores"
TOOLS = TOOLS_ROOT


def species_dir(species: str, regime_suffix: str) -> Path:
    return DATA / f"{species}{regime_suffix}"


def out_path(tool: str, species: str, regime_suffix: str) -> Path:
    reg = "paralog" if regime_suffix else "one2one"
    return OUT / f"{tool}.{species}.{reg}.json"


def log_path(tool: str, species: str, regime_suffix: str) -> Path:
    reg = "paralog" if regime_suffix else "one2one"
    return LOG / f"{tool}.{species}.{reg}.log"


def workdir(tool: str, species: str, regime_suffix: str) -> Path:
    reg = "paralog" if regime_suffix else "one2one"
    d = WORK / f"{tool}.{species}.{reg}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def load_pairs_proteins(species: str, regime_suffix: str) -> dict[str, str]:
    """Load protein FASTA for the batch; key=protein_id, value=seq."""
    fa = species_dir(species, regime_suffix) / "proteins.fa"
    seqs = {}
    cur_id, cur = None, []
    with fa.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(cur)
                cur_id = line[1:].strip().split()[0]
                cur = []
            else:
                cur.append(line.strip())
    if cur_id:
        seqs[cur_id] = "".join(cur)
    return seqs


def write_fasta(path: Path, seqs: dict[str, str]) -> None:
    with path.open("w") as fh:
        for pid, seq in seqs.items():
            fh.write(f">{pid}\n{seq}\n")


def save_batch_result(tool, species, regime_suffix, results, wall_time, extra=None):
    """Save in V3-compatible format."""
    reg = "paralog" if regime_suffix else "one2one"
    n_total = len(results)
    n_with_hits = sum(1 for h in results.values() if h)
    payload = {
        "tool": tool,
        "species": species,
        "regime": reg,
        "n_proteins": n_total,
        "n_with_hits": n_with_hits,
        "wall_time_sec": round(wall_time, 1),
        "results": results,
    }
    if extra:
        payload.update(extra)
    op = out_path(tool, species, regime_suffix)
    op.parent.mkdir(parents=True, exist_ok=True)
    with op.open("w") as fh:
        json.dump(payload, fh, indent=2)
    return op


def run_cmd(cmd, env=None, timeout=None, cwd=None, conda_env=None,
            stdin=None, capture=True, check=False):
    """subprocess.run with conda_env wrapping."""
    if conda_env:
        env_name = os.environ.get("RESCUE_CONDA_ENV", conda_env)
        if env_name:
            cmd = [CONDA, "run", "-n", env_name, "--no-capture-output"] + list(cmd)
    full_env = os.environ.copy()
    if env:
        full_env.update(env)
    return subprocess.run(
        cmd,
        env=full_env,
        cwd=cwd,
        input=stdin,
        capture_output=capture,
        text=True,
        timeout=timeout,
        check=check,
    )


# ------- GFF parsers (per tool, returning V3 hit format) ---------------

def parse_spaln_gff(gff_text: str) -> dict[str, list[dict]]:
    """Parse Spaln GFF3 output. Spaln uses gene/mRNA/CDS; Target in mRNA attr."""
    results = defaultdict(list)
    cur, cds = None, []
    rank_counter = defaultdict(int)
    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 9:
            continue
        chrom, _, feat, s, e, score, strand, _, attrs = p[:9]
        s, e = int(s), int(e)
        if feat == "mRNA" or feat == "gene":
            if cur and cds:
                cur["cds"] = sorted(cds, key=lambda c: c["start"])
                results[cur.pop("_pid")].append(cur)
            ad = dict(kv.split("=", 1) for kv in attrs.split(";") if "=" in kv)
            tgt = ad.get("Target", "").split()
            pid = tgt[0] if tgt else ad.get("ID", "")
            q_start = int(tgt[1]) if len(tgt) > 1 else 1
            q_end = int(tgt[2]) if len(tgt) > 2 else 0
            rank_counter[pid] += 1
            sval = float(score) if score not in (".", "") else 0.0
            cur = {
                "rank": rank_counter[pid],
                "chrom": chrom, "strand": strand, "score": sval,
                "identity": float(ad.get("Identity", 0)),
                "query_start": q_start, "query_end": q_end,
                "_pid": pid,
            }
            cds = []
        elif feat in ("CDS", "cds") and cur is not None:
            cds.append({"chrom": chrom, "start": s, "end": e})
    if cur and cds:
        cur["cds"] = sorted(cds, key=lambda c: c["start"])
        results[cur.pop("_pid")].append(cur)
    return {k: sorted(v, key=lambda h: h["rank"]) for k, v in results.items()}


def parse_generic_gff3(gff_text: str, fallback_pid: str | None = None) -> list[dict]:
    """For per-query tools: return a list of hits (single protein) from GFF3.

    Groups CDS under the same mRNA/gene ID. Scores read from column 6.
    """
    transcripts = {}
    order = []
    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 9:
            continue
        chrom, _, feat, s, e, score, strand, _, attrs = p[:9]
        ad = dict(kv.split("=", 1) for kv in attrs.split(";") if "=" in kv)
        if feat in ("mRNA", "gene", "transcript", "match"):
            tid = ad.get("ID") or ad.get("Name") or fallback_pid or "tx"
            if tid not in transcripts:
                order.append(tid)
            sval = float(score) if score not in (".", "") else 0.0
            transcripts.setdefault(tid, {
                "chrom": chrom, "strand": strand, "score": sval,
                "identity": float(ad.get("Identity", 0)),
                "query_start": 0, "query_end": 0, "cds": [],
            })
            tgt = ad.get("Target", "").split()
            if len(tgt) > 2:
                transcripts[tid]["query_start"] = int(tgt[1])
                transcripts[tid]["query_end"] = int(tgt[2])
        elif feat in ("CDS", "exon", "match_part", "protein_match", "coding_exon"):
            parent = ad.get("Parent") or ad.get("ID") or fallback_pid or "tx"
            if parent not in transcripts:
                order.append(parent)
                sval = float(score) if score not in (".", "") else 0.0
                transcripts[parent] = {
                    "chrom": chrom, "strand": strand, "score": sval,
                    "identity": 0.0, "query_start": 0, "query_end": 0, "cds": [],
                }
            transcripts[parent]["cds"].append(
                {"chrom": chrom, "start": int(s), "end": int(e)}
            )
    hits = []
    for r, tid in enumerate(order, 1):
        t = transcripts[tid]
        if not t["cds"]:
            continue
        t["cds"].sort(key=lambda c: c["start"])
        t["rank"] = r
        hits.append(t)
    return sorted(hits, key=lambda h: h["rank"])
