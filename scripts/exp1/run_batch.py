#!/usr/bin/env python3
"""Run one (tool, species, regime) batch for Exp 1.

Produces results/exp1/raw/{tool}.{species}.{regime}.json in V3 hit format.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))

import yaml  # noqa: E402
from common import (  # noqa: E402
    CONDA, GENOMES, TOOLS, SCORES,
    species_dir, out_path, log_path, workdir,
    load_pairs_proteins, write_fasta,
    save_batch_result, run_cmd,
    parse_miniprot_gff, parse_spaln_gff, parse_generic_gff3,
)
from rescue_paths import TBLASTN  # noqa: E402

CFG = yaml.safe_load((Path(__file__).resolve().parents[2] / "configs/exp1_tools.yaml").read_text())
THREADS = CFG["defaults"]["threads"]
MAX_PROT = int(os.environ.get("EXP1_MAX_PROT", "0") or 0)
EXP1_CONDA_ENV = os.environ.get("RESCUE_CONDA_ENV", "")


def _cap(d: dict) -> dict:
    if MAX_PROT > 0:
        return dict(list(d.items())[:MAX_PROT])
    return d
BATCH_TIMEOUT = CFG["defaults"]["batch_timeout_sec"]
PQ_TIMEOUT = CFG["defaults"]["per_query_timeout_sec"]


# ---------- batch-mode tools ------------------------------------------------

def _miniprot_cmd(genome_mpi: str, prot_fa: str, extra: list[str] | None = None):
    cmd = ["miniprot", "--gff", "-t", str(THREADS), "-N", "1"]
    if extra:
        cmd += extra
    cmd += [genome_mpi, prot_fa]
    return cmd


def run_miniprot(species, regime_suffix, wdir):
    prots = _cap(load_pairs_proteins(species, regime_suffix))
    fa = wdir / "proteins.fa"; write_fasta(fa, prots)
    mpi = str(GENOMES / f"{species}.mpi")
    t0 = time.time()
    res = run_cmd(_miniprot_cmd(mpi, str(fa)), conda_env=EXP1_CONDA_ENV,
                  timeout=BATCH_TIMEOUT)
    (wdir / "miniprot.gff").write_text(res.stdout)
    (wdir / "miniprot.stderr").write_text(res.stderr or "")
    results = parse_miniprot_gff(res.stdout)
    for pid in prots:
        results.setdefault(pid, [])
    return results, time.time() - t0


def _ensure_minisplice_score(species: str) -> Path:
    """Lazily generate per-genome minisplice score.tsv (cached)."""
    SCORES.mkdir(parents=True, exist_ok=True)
    out = SCORES / f"{species}.score.tsv"
    if out.exists() and out.stat().st_size > 0:
        return out
    model = TOOLS / "minisplice/models/vi2-7k.kan"
    cali = TOOLS / "minisplice/models/vi2-7k.kan.cali"
    genome = GENOMES / f"{species}.fa"
    binp = TOOLS / "minisplice/minisplice"
    cmd = [str(binp), "predict", "-t", str(THREADS),
           "-c", str(cali), str(model), str(genome)]
    with out.open("w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE, text=True,
                           timeout=BATCH_TIMEOUT)
    if r.returncode != 0:
        (SCORES / f"{species}.score.err").write_text(r.stderr or "")
        out.unlink(missing_ok=True)
        raise RuntimeError(f"minisplice predict failed for {species}: {r.stderr[:500]}")
    return out


def run_minisplice_miniprot(species, regime_suffix, wdir):
    prots = _cap(load_pairs_proteins(species, regime_suffix))
    fa = wdir / "proteins.fa"; write_fasta(fa, prots)
    score = _ensure_minisplice_score(species)
    mpi = str(GENOMES / f"{species}.mpi")
    extra = ["-j", "2", f"--spsc={score}"]
    t0 = time.time()
    res = run_cmd(_miniprot_cmd(mpi, str(fa), extra), conda_env=EXP1_CONDA_ENV,
                  timeout=BATCH_TIMEOUT)
    (wdir / "miniprot.gff").write_text(res.stdout)
    (wdir / "miniprot.stderr").write_text(res.stderr or "")
    results = parse_miniprot_gff(res.stdout)
    for pid in prots:
        results.setdefault(pid, [])
    return results, time.time() - t0


def run_spaln3(species, regime_suffix, wdir):
    prots = _cap(load_pairs_proteins(species, regime_suffix))
    fa = wdir / "proteins.fa"; write_fasta(fa, prots)
    prefix = str(GENOMES / species)  # uses {prefix}.bkp/.seq/.idx/...
    cmd = ["spaln", "-Q7", "-O0", "-M1,2", f"-t{THREADS}",
           f"-d{prefix}", str(fa)]
    t0 = time.time()
    res = run_cmd(cmd, conda_env=EXP1_CONDA_ENV, timeout=BATCH_TIMEOUT)
    (wdir / "spaln.gff").write_text(res.stdout)
    (wdir / "spaln.stderr").write_text(res.stderr or "")
    results = parse_spaln_gff(res.stdout)
    for pid in prots:
        results.setdefault(pid, [])
    return results, time.time() - t0


def run_gth(species, regime_suffix, wdir):
    """GTH whole-genome batch; -gff3out; -paralogs off -intermediate."""
    prots = _cap(load_pairs_proteins(species, regime_suffix))
    fa = wdir / "proteins.fa"; write_fasta(fa, prots)
    genome = str(GENOMES / f"{species}.fa")
    sp_model = CFG["gth_species_map"].get(species, "human")
    gff = wdir / "gth.gff"
    cmd = ["gth", "-genomic", genome, "-protein", str(fa),
           "-species", sp_model, "-prseedlength", "4", "-prminmatchlen", "8",
           "-gff3out", "-skipalignmentout", "-o", str(gff), "-force"]
    t0 = time.time()
    res = run_cmd(cmd, conda_env=EXP1_CONDA_ENV, timeout=BATCH_TIMEOUT)
    (wdir / "gth.stdout").write_text(res.stdout or "")
    (wdir / "gth.stderr").write_text(res.stderr or "")
    gff_text = gff.read_text() if gff.exists() else ""
    # GTH gff3 format: group mRNAs under protein; use Target attr to route.
    results = _parse_gth_gff3_per_protein(gff_text)
    for pid in prots:
        results.setdefault(pid, [])
    return results, time.time() - t0


def _parse_gth_gff3_per_protein(gff_text: str):
    """GTH GFF3 emits one mRNA per alignment with Target=protein_id ...; we
    group by the Target protein id and keep rank order of appearance."""
    from collections import defaultdict
    res = defaultdict(list)
    cur, cds, pid = None, [], None
    rank = defaultdict(int)
    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 9:
            continue
        chrom, _, feat, s, e, score, strand, _, attrs = p[:9]
        ad = dict(kv.split("=", 1) for kv in attrs.split(";") if "=" in kv)
        if feat in ("mRNA", "gene", "prediction"):
            if cur and cds and pid:
                cur["cds"] = sorted(cds, key=lambda c: c["start"])
                res[pid].append(cur)
            tgt = ad.get("Target", "").split()
            pid = tgt[0] if tgt else ""
            rank[pid] += 1
            cur = {
                "rank": rank[pid],
                "chrom": chrom, "strand": strand,
                "score": float(score) if score not in (".", "") else 0.0,
                "identity": float(ad.get("Identity", 0)),
                "query_start": int(tgt[1]) if len(tgt) > 1 else 0,
                "query_end": int(tgt[2]) if len(tgt) > 2 else 0,
            }
            cds = []
        elif feat in ("CDS", "exon") and cur is not None:
            cds.append({"chrom": chrom, "start": int(s), "end": int(e)})
    if cur and cds and pid:
        cur["cds"] = sorted(cds, key=lambda c: c["start"])
        res[pid].append(cur)
    return {k: sorted(v, key=lambda h: h["rank"]) for k, v in res.items()}


# ---------- per-query tools ---------------------------------------------

def _write_single_fa(pid, seq, wdir):
    fa = wdir / f"{pid}.fa"
    fa.write_text(f">{pid}\n{seq}\n")
    return fa


# ---- batched per-query pipelines (one search per genome, not per protein) ----

LEGACY_BLAST = Path(os.environ.get("RESCUE_LEGACY_BLAST_BIN_DIR", "tools/blastlegacy/bin"))
LEGACY_DB = Path(os.environ.get("RESCUE_LEGACY_BLAST_DB_DIR", GENOMES / "legacy"))


def _run_batch_blat(genome_fa, batch_fa, out_psl, log_stub):
    """Run BLAT with multi-protein FASTA input. One invocation for all."""
    cmd = ["blat", "-t=dnax", "-q=prot", "-noHead",
           "-minIdentity=18", "-minScore=15",
           str(genome_fa), str(batch_fa), str(out_psl)]
    with open(f"{log_stub}.stderr", "w") as ferr:
        r = subprocess.run([CONDA, "run",
                            "-n", "scipio", "--no-capture-output"] + cmd,
                           stdout=ferr, stderr=subprocess.STDOUT,
                           timeout=BATCH_TIMEOUT)
    return r.returncode == 0 and out_psl.exists() and out_psl.stat().st_size > 0


def _read_fasta_dict(path):
    """Read FASTA file into {id: sequence} dict."""
    seqs = {}
    cur_id = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(buf)
                cur_id = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur_id:
        seqs[cur_id] = "".join(buf)
    return seqs


def _get_available_mem_gb():
    with open("/proc/meminfo") as f:
        for line in f:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / (1024 * 1024)
    return 4.0


def _run_parallel_blat(genome_fa, query_fa, out_psl, wdir, label="", max_chunks=16,
                       mem_budget_gb=None):
    """Split query proteins into chunks, run parallel BLAT, merge PSL.

    Auto-detects chunk count based on available memory and genome size.
    Resumes from existing PSL if present (skips already-done proteins).
    mem_budget_gb: max RAM this batch may use for BLAT (avoids races with concurrent batches).
    """
    genome_fa = Path(genome_fa)
    query_fa = Path(query_fa)
    out_psl = Path(out_psl)

    all_prots = _read_fasta_dict(query_fa)
    if not all_prots:
        out_psl.write_text("")
        return

    done_pids = set()
    existing_psl_lines = []
    if out_psl.exists() and out_psl.stat().st_size > 0:
        with open(out_psl) as f:
            for line in f:
                existing_psl_lines.append(line)
                parts = line.split("\t")
                if len(parts) >= 10:
                    done_pids.add(parts[9])

    remaining = {pid: seq for pid, seq in all_prots.items() if pid not in done_pids}
    if not remaining:
        print(f"  [{label}] BLAT already complete ({len(done_pids)}/{len(all_prots)} proteins)")
        sys.stdout.flush()
        return

    if done_pids:
        print(f"  [{label}] BLAT resume: {len(done_pids)} done, {len(remaining)} remaining")
    else:
        print(f"  [{label}] BLAT starting: {len(remaining)} proteins")
    sys.stdout.flush()

    genome_size_gb = genome_fa.stat().st_size / (1024 ** 3)
    blat_mem_gb = genome_size_gb * 2.7 + 0.3
    avail_gb = _get_available_mem_gb()
    if mem_budget_gb is None:
        safe_gb = max(avail_gb - 4, 2)
    else:
        safe_gb = mem_budget_gb
    n_chunks = max(1, min(max_chunks, int(safe_gb / blat_mem_gb)))
    n_chunks = min(n_chunks, len(remaining))

    pids = list(remaining.keys())
    chunk_size = max(1, -(-len(pids) // n_chunks))
    chunks = []
    for i in range(0, len(pids), chunk_size):
        chunk_pids = pids[i:i + chunk_size]
        idx = len(chunks)
        chunk_fa = wdir / f"blat_chunk_{idx}.fa"
        with open(chunk_fa, "w") as f:
            for pid in chunk_pids:
                f.write(f">{pid}\n{remaining[pid]}\n")
        chunk_psl = wdir / f"blat_chunk_{idx}.psl"
        chunks.append((idx, chunk_fa, chunk_psl, len(chunk_pids)))

    print(f"  [{label}] Parallel BLAT: {len(remaining)} proteins -> {len(chunks)} chunks "
          f"(~{blat_mem_gb:.1f}GB/proc, {avail_gb:.0f}GB avail)")
    sys.stdout.flush()

    ooc_path = genome_fa.with_suffix(".ooc")
    ooc_args = [f"-ooc={ooc_path}"] if ooc_path.exists() else []
    if ooc_args:
        print(f"  [{label}] using BLAT overused-tile file {ooc_path}")
        sys.stdout.flush()

    t_blat = time.time()
    procs = []
    for idx, chunk_fa, chunk_psl, n_prot in chunks:
        cmd = [CONDA, "run",
               "-n", "scipio", "--no-capture-output",
               "blat", "-t=dnax", "-q=prot", "-noHead",
               "-minIdentity=18", "-minScore=15",
               *ooc_args,
               str(genome_fa), str(chunk_fa), str(chunk_psl)]
        stderr_f = open(wdir / f"blat_chunk_{idx}.stderr", "w")
        p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=stderr_f)
        procs.append((p, stderr_f, idx, chunk_psl, n_prot))

    finished = set()
    while len(finished) < len(procs):
        for p, stderr_f, idx, chunk_psl, n_prot in procs:
            if idx in finished:
                continue
            if p.poll() is not None:
                stderr_f.close()
                finished.add(idx)
                elapsed = time.time() - t_blat
                sz = chunk_psl.stat().st_size if chunk_psl.exists() else 0
                print(f"  [{label}] BLAT chunk {len(finished)}/{len(chunks)} done "
                      f"(chunk {idx}: {n_prot} proteins, {sz} bytes, {elapsed:.0f}s elapsed)")
                sys.stdout.flush()
        if len(finished) < len(procs):
            time.sleep(30)

    with open(out_psl, "w") as out:
        for line in existing_psl_lines:
            out.write(line)
        for _, _, idx, chunk_psl, _ in procs:
            if chunk_psl.exists():
                with open(chunk_psl) as inp:
                    shutil.copyfileobj(inp, out)

    total_time = time.time() - t_blat
    total_queries = 0
    with open(out_psl) as f:
        seen = set()
        for line in f:
            parts = line.split("\t")
            if len(parts) >= 10:
                seen.add(parts[9])
        total_queries = len(seen)
    print(f"  [{label}] BLAT done: {total_queries}/{len(all_prots)} proteins, "
          f"{out_psl.stat().st_size} bytes, {total_time:.0f}s")
    sys.stdout.flush()


def _run_batch_blastall(genome_fa, batch_fa, out_txt, log_stub, threads=8):
    """BLAST+ tblastn, producing m8 tabular output compatible with procompart.

    Legacy blastall 2.2.26 from bioconda blast-legacy was catastrophically slow
    under system load (42/2000 queries in 2.5h on chicken genome). BLAST+ 2.17
    tblastn with -outfmt 6 produces the same 12-column tabular format and has
    far better multi-thread scaling.
    """
    species = genome_fa.stem
    legacy_db = LEGACY_DB / species  # v4 indexes work for BLAST+
    cmd = [TBLASTN, "-query", str(batch_fa),
           "-db", str(legacy_db), "-out", str(out_txt),
           "-num_threads", str(threads),
           "-outfmt", "6", "-evalue", "10", "-seg", "no"]
    with open(f"{log_stub}.stderr", "w") as ferr:
        r = subprocess.run(cmd, stdout=ferr, stderr=subprocess.STDOUT,
                           timeout=BATCH_TIMEOUT)
    return r.returncode == 0 and out_txt.exists() and out_txt.stat().st_size > 0


def _scipio_per_protein(args):
    pid, seq, wdir, genome_fa, shared_psl = args
    try:
        fa = wdir / f"{pid}.fa"
        fa.write_text(f">{pid}\n{seq}\n")
        # Extract this protein's PSL rows only (scipio will still work with subset)
        sub_psl = wdir / f"{pid}.psl"
        with open(shared_psl) as fh, sub_psl.open("w") as out:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 10 and parts[9] == pid:
                    out.write(line)
        if sub_psl.stat().st_size == 0:
            return pid, []
        yaml_out = wdir / f"{pid}.yaml"
        if not (yaml_out.exists() and yaml_out.stat().st_size > 0):
            cmd = ["perl", str(TOOLS / "scipio/scipio-1.4/scipio.1.4.1.pl"),
                   f"--blat_output={sub_psl}",
                   "--min_score=0.1", "--min_identity=20", "--min_coverage=20",
                   str(genome_fa), str(fa)]
            r = subprocess.run([CONDA, "run",
                                "-n", "scipio", "--no-capture-output"] + cmd,
                               capture_output=True, text=True, timeout=PQ_TIMEOUT)
            yaml_out.write_text(r.stdout or "")
            if not r.stdout.strip():
                return pid, []
        if yaml_out.stat().st_size == 0:
            return pid, []
        cmd2 = ["perl", str(TOOLS / "scipio/scipio-1.4/yaml2gff.1.4.pl"),
                str(yaml_out)]
        r2 = subprocess.run([CONDA, "run",
                             "-n", "scipio", "--no-capture-output"] + cmd2,
                            capture_output=True, text=True, timeout=120)
        return pid, parse_generic_gff3(r2.stdout or "", fallback_pid=pid)
    except Exception:
        return pid, []


BLAT_MEM_BUDGET = float(os.environ.get("BLAT_MEM_BUDGET_GB", "0") or 0) or None


def run_scipio(species, regime_suffix, wdir, workers=8):
    prots = _cap(load_pairs_proteins(species, regime_suffix))
    genome = GENOMES / f"{species}.fa"
    batch_fa = wdir / "proteins.fa"; write_fasta(batch_fa, prots)
    shared_psl = wdir / "shared.psl"
    label = f"scipio/{species}/{regime_suffix}"
    t0 = time.time()

    _run_parallel_blat(genome, batch_fa, shared_psl, wdir, label=label,
                       mem_budget_gb=BLAT_MEM_BUDGET)

    print(f"  [{label}] Per-protein Scipio phase: {len(prots)} proteins, {workers} workers")
    sys.stdout.flush()
    ctx = mp.get_context("spawn")
    args = [(pid, seq, wdir, str(genome), str(shared_psl))
            for pid, seq in prots.items()]
    results = {pid: [] for pid in prots}
    t_scipio = time.time()
    done_count = 0
    hit_count = 0
    with ctx.Pool(workers) as pool:
        for pid, hits in pool.imap_unordered(_scipio_per_protein, args, chunksize=10):
            results[pid] = hits
            done_count += 1
            if hits:
                hit_count += 1
            if done_count % 100 == 0 or done_count == len(prots):
                elapsed = time.time() - t_scipio
                rate = done_count / elapsed if elapsed > 0 else 0
                eta = (len(prots) - done_count) / rate if rate > 0 else 0
                print(f"  [{label}] scipio {done_count}/{len(prots)} "
                      f"({hit_count} hits, {elapsed:.0f}s, ETA {eta:.0f}s)")
                sys.stdout.flush()
    total = time.time() - t0
    blat_time = t_scipio - t0
    scipio_time = time.time() - t_scipio
    print(f"  [{label}] TOTAL: {total:.0f}s (BLAT={blat_time:.0f}s, scipio={scipio_time:.0f}s)")
    sys.stdout.flush()
    return results, total


def _genblastg_setup_chunk_dir(chunk_dir, tg_bundle):
    """Install blastall/formatdb/alignscore symlinks inside a genblastG work dir."""
    for tool_name, target in (("blastall", tg_bundle / "blastall_threaded.sh"),
                               ("formatdb", tg_bundle / "formatdb")):
        link = chunk_dir / tool_name
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(target)
    aln = chunk_dir / "alignscore.txt"
    if aln.is_symlink() or aln.exists():
        aln.unlink()
    aln.symlink_to(tg_bundle / "alignscore.txt")


def _parse_genblastg_gff_dir(gff_dir, results):
    """Parse genblastg_out*.gff in gff_dir and append hits to results dict."""
    for gff_file in gff_dir.glob("genblastg_out*.gff"):
        text = gff_file.read_text()
        buf_by_pid = {}
        for line in text.splitlines():
            if not line or line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 9:
                continue
            attrs = dict(kv.split("=", 1) for kv in p[8].split(";") if "=" in kv)
            name = attrs.get("ID", attrs.get("Parent", ""))
            pid = name.split("-R")[0] if "-R" in name else name.split(".")[0]
            buf_by_pid.setdefault(pid, []).append(line)
        for pid, lines in buf_by_pid.items():
            results.setdefault(pid, []).extend(
                parse_generic_gff3("\n".join(lines), fallback_pid=pid)
            )


def _genblastg_batch(species, regime_suffix, wdir):
    """genBlastG: split proteins into N chunks, run genblastG in parallel
    sub-workdirs, merge GFFs.

    Chain/alignment phase is single-threaded inside genblastG; running multiple
    instances with disjoint query sets parallelizes it wall-clock-wise. Blast
    phase uses the -a 8 wrapper per instance.

    Env var GENBLASTG_CHUNKS overrides default chunk count
    (4 for one2one, 2 for paralog).

    Uses the thecgs matched-version bundle (blastall 2.2.19 + genblastG
    v1.0.138 + alignscore.txt). bioconda's genblastg + blast-legacy 2.2.26
    triggers "wuBlast/Blast file error".
    """
    prots = _cap(load_pairs_proteins(species, regime_suffix))
    tg_bundle = TOOLS / "genblast/thecgs"
    legacy_db = LEGACY_DB / f"{species}.fa"

    reg = "paralog" if regime_suffix else "one2one"
    default = "4" if reg == "one2one" else "2"
    n_chunks = int(os.environ.get("GENBLASTG_CHUNKS", default))
    n_chunks = max(1, min(n_chunks, len(prots)))

    # Full FASTA for reference (parsers, debugging).
    batch_fa_all = wdir / "proteins.fa"; write_fasta(batch_fa_all, prots)

    pids = list(prots.keys())
    chunk_size = -(-len(pids) // n_chunks)
    chunk_infos = []  # (idx, chunk_dir, chunk_fa, n_prot)
    for i in range(0, len(pids), chunk_size):
        idx = len(chunk_infos)
        chunk_pids = pids[i:i + chunk_size]
        chunk_dir = wdir / f"chunk_{idx}"
        chunk_dir.mkdir(parents=True, exist_ok=True)
        chunk_fa = chunk_dir / "proteins.fa"
        write_fasta(chunk_fa, {pid: prots[pid] for pid in chunk_pids})
        chunk_infos.append((idx, chunk_dir, chunk_fa, len(chunk_pids)))

    label = f"genblastg/{species}/{reg}"
    print(f"  [{label}] parallel genblastG: {len(prots)} proteins -> "
          f"{len(chunk_infos)} chunks (~{chunk_size}/chunk)")
    sys.stdout.flush()

    t0 = time.time()
    procs = []  # (Popen | None, stdout_f, stderr_f, idx, chunk_dir, n_prot)
    for idx, chunk_dir, chunk_fa, n_prot in chunk_infos:
        # Resume: if chunk has nontrivial GFF already, don't relaunch.
        existing_gffs = list(chunk_dir.glob("genblastg_out*.gff"))
        if existing_gffs and any(g.stat().st_size > 1024 for g in existing_gffs):
            print(f"  [{label}] chunk {idx} already has GFF; skipping genblastG")
            sys.stdout.flush()
            procs.append((None, None, None, idx, chunk_dir, n_prot))
            continue

        _genblastg_setup_chunk_dir(chunk_dir, tg_bundle)
        out_pref = chunk_dir / "genblastg_out"
        cmd = [str(tg_bundle / "genblastG"),
               "-p", "genblastg", "-q", str(chunk_fa), "-t", str(legacy_db),
               "-o", str(out_pref), "-gff"]
        env = os.environ.copy()
        env["PATH"] = f"{chunk_dir}:{env.get('PATH','')}"
        fo = open(chunk_dir / "genblastg.stdout", "w")
        fe = open(chunk_dir / "genblastg.stderr", "w")
        p = subprocess.Popen(cmd, cwd=str(chunk_dir), stdout=fo, stderr=fe, env=env)
        procs.append((p, fo, fe, idx, chunk_dir, n_prot))
        print(f"  [{label}] chunk {idx} launched ({n_prot} proteins, PID {p.pid})")
        sys.stdout.flush()

    finished = set()
    # Pre-skipped chunks are already done.
    for p, _, _, idx, _, _ in procs:
        if p is None:
            finished.add(idx)

    last_heartbeat = time.time()
    deadline = t0 + BATCH_TIMEOUT
    while len(finished) < len(procs):
        for p, fo, fe, idx, chunk_dir, n_prot in procs:
            if idx in finished:
                continue
            if p.poll() is not None:
                fo.close(); fe.close()
                finished.add(idx)
                elapsed = time.time() - t0
                gffs = list(chunk_dir.glob("genblastg_out*.gff"))
                sz = sum(g.stat().st_size for g in gffs) if gffs else 0
                print(f"  [{label}] chunk {idx} done "
                      f"({len(finished)}/{len(procs)}): {n_prot} prot, "
                      f"{sz} B GFF, {elapsed:.0f}s elapsed, rc={p.returncode}")
                sys.stdout.flush()
        if len(finished) < len(procs):
            now = time.time()
            if now > deadline:
                for p, fo, fe, idx, _, _ in procs:
                    if idx not in finished and p is not None:
                        p.kill()
                        fo.close(); fe.close()
                raise subprocess.TimeoutExpired("genblastg_parallel", BATCH_TIMEOUT)
            if now - last_heartbeat > 600:
                sizes = []
                for p, _, _, idx, chunk_dir, _ in procs:
                    if idx in finished:
                        continue
                    gffs = list(chunk_dir.glob("genblastg_out*.gff"))
                    sz = sum(g.stat().st_size for g in gffs) if gffs else 0
                    sizes.append(f"c{idx}:{sz // 1_000_000}MB")
                print(f"  [{label}] heartbeat {now - t0:.0f}s elapsed, "
                      f"{len(procs) - len(finished)} chunks running "
                      f"({', '.join(sizes)})")
                sys.stdout.flush()
                last_heartbeat = now
            time.sleep(60)

    results = {pid: [] for pid in prots}
    for _, _, _, idx, chunk_dir, _ in procs:
        _parse_genblastg_gff_dir(chunk_dir, results)
    return results, time.time() - t0


def _prosplign_batch(species, regime_suffix, wdir, workers=8):
    """tblastn (blastall) batch -> procompart -> prosplign per compartment."""
    prots = _cap(load_pairs_proteins(species, regime_suffix))
    genome = GENOMES / f"{species}.fa"
    batch_fa = wdir / "proteins.fa"; write_fasta(batch_fa, prots)
    tblastn_out = wdir / "tblastn.txt"
    compart_out = wdir / "compart.txt"
    t0 = time.time()
    if not tblastn_out.exists() or tblastn_out.stat().st_size == 0:
        print(f"  [prosplign] blastall (batch) -> {tblastn_out}")
        ok = _run_batch_blastall(genome, batch_fa, tblastn_out,
                                 str(wdir / "blastall"), threads=8)
        if not ok:
            raise RuntimeError("blastall failed; see blastall.stderr")
        print(f"  [prosplign] blastall done ({tblastn_out.stat().st_size} bytes)")
    # procompart once on whole output
    pc = TOOLS / "prosplign/procompart"
    env = os.environ.copy()
    env["LD_LIBRARY_PATH"] = f"{TOOLS/'prosplign/lib'}:{env.get('LD_LIBRARY_PATH','')}"
    with open(tblastn_out) as fi, compart_out.open("w") as fo, \
         open(wdir / "procompart.stderr", "w") as fe:
        r = subprocess.run([str(pc), "-t"], stdin=fi, stdout=fo, stderr=fe,
                           env=env, timeout=BATCH_TIMEOUT)
    print(f"  [prosplign] procompart -> {compart_out.stat().st_size} bytes")

    # prosplign on merged compart file; outputs alignments for all queries
    prosplign_out = wdir / "prosplign.txt"
    ps = TOOLS / "prosplign/prosplign"
    with open(wdir / "prosplign.stderr", "w") as fe:
        r = subprocess.run([str(ps), "-nogenbank", "-two_stages", "-full",
                            "-i", str(compart_out), "-nfa", str(genome),
                            "-pfa", str(batch_fa), "-o", str(prosplign_out)],
                           env=env, stderr=fe, stdout=subprocess.DEVNULL,
                           timeout=BATCH_TIMEOUT)
    # Parse prosplign text output per protein (simple exon blocks).
    import re
    text = prosplign_out.read_text() if prosplign_out.exists() else ""
    results = {pid: [] for pid in prots}
    cur_pid = None; cur_exons = []; cur_chrom = None; cur_strand = "+"
    def _flush():
        if cur_pid and cur_exons:
            hits = [{
                "rank": 1,
                "chrom": cur_chrom, "strand": cur_strand, "score": 0.0,
                "identity": 0.0,
                "query_start": 0, "query_end": 0,
                "cds": sorted([{"chrom": c, "start": s, "end": e}
                               for c,s,e in cur_exons], key=lambda x: x["start"]),
            }]
            results.setdefault(cur_pid, []).extend(hits)
    for line in text.splitlines():
        mh = re.match(r">\s*(\S+)", line)
        if mh:
            _flush()
            cur_pid = mh.group(1); cur_exons = []; cur_chrom = None
            continue
        me = re.match(r"\s*exon\s+(\S+)\s+(\d+)\s+(\d+)\s+([+-])", line)
        if me:
            c, s, e, st = me.group(1), int(me.group(2)), int(me.group(3)), me.group(4)
            cur_chrom = c; cur_strand = st
            cur_exons.append((c, min(s, e), max(s, e)))
    _flush()
    return results, time.time() - t0


PER_QUERY_FUNCS = {
    "scipio": run_scipio,
    "procompart_prosplign": _prosplign_batch,
    "genblastg": _genblastg_batch,
}


def run_per_query(tool, species, regime_suffix, wdir, workers=8):
    func = PER_QUERY_FUNCS[tool]
    # scipio uses pool internally; others are single-process batched
    try:
        return func(species, regime_suffix, wdir, workers=workers)
    except TypeError:
        return func(species, regime_suffix, wdir)


# ---------- dispatch --------------------------------------------------------

BATCH_FUNCS = {
    "miniprot": run_miniprot,
    "spaln3": run_spaln3,
    "minisplice_miniprot": run_minisplice_miniprot,
    "gth": run_gth,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", required=True, choices=list(CFG["tools"].keys()))
    ap.add_argument("--species", required=True, choices=CFG["species"])
    ap.add_argument("--regime", required=True, choices=list(CFG["regimes"].keys()))
    ap.add_argument("--workers", type=int, default=THREADS)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    suffix = CFG["regimes"][args.regime]
    op = out_path(args.tool, args.species, suffix)
    if op.exists() and not args.force and op.stat().st_size > 100:
        print(f"SKIP (exists): {op}")
        return
    wd = workdir(args.tool, args.species, suffix)
    print(f"[exp1] {args.tool} / {args.species} / {args.regime} -> {op}")

    tool_cfg = CFG["tools"][args.tool]
    if tool_cfg["mode"] == "batch":
        results, wall = BATCH_FUNCS[args.tool](args.species, suffix, wd)
    else:
        results, wall = run_per_query(args.tool, args.species, suffix, wd,
                                      workers=args.workers)
    save_batch_result(args.tool, args.species, suffix, results, wall)
    n_hits = sum(1 for h in results.values() if h)
    print(f"  done: {n_hits}/{len(results)} proteins with >=1 hit; {wall:.1f}s")


if __name__ == "__main__":
    main()
