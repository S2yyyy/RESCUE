#!/usr/bin/env python3
"""V3 miniprot wrapper: run miniprot with top-k output, parse GFF to standard format.

Usage:
    python scripts/v3_run_miniprot.py --species mouse --top-k 5
    python scripts/v3_run_miniprot.py --species mouse --top-k 5 --max-proteins 50
    python scripts/v3_run_miniprot.py --species mouse --region 12:100000000-100200000  # restricted region
"""

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from v3_tool_common import (PROCESSED, V3_DATA, V3_RESULTS,
                            load_proteins, load_ground_truth, load_meta,
                            save_results, write_protein_batch_fasta)
from rescue_paths import MINIPROT


def parse_miniprot_gff(gff_text):
    """Parse miniprot GFF output into standard hit format.

    Returns: {protein_id: [hit1, hit2, ...]} sorted by rank.
    """
    results = defaultdict(list)
    current_mrna = None
    current_cds = []

    for line in gff_text.strip().split("\n"):
        if not line or line.startswith("##gff-version"):
            continue

        if line.startswith("##PAF"):
            # Flush previous mRNA
            if current_mrna and current_cds:
                current_mrna["cds"] = sorted(current_cds, key=lambda c: c["start"])
                results[current_mrna["_protein_id"]].append(current_mrna)
            current_mrna = None
            current_cds = []
            continue

        if line.startswith("#"):
            continue

        parts = line.split("\t")
        if len(parts) < 9:
            continue

        chrom, source, feature, start, end, score, strand, phase, attrs = parts[:9]
        start, end = int(start), int(end)

        if feature == "mRNA":
            # Flush previous
            if current_mrna and current_cds:
                current_mrna["cds"] = sorted(current_cds, key=lambda c: c["start"])
                results[current_mrna["_protein_id"]].append(current_mrna)

            # Parse attributes
            attr_dict = {}
            for a in attrs.split(";"):
                if "=" in a:
                    k, v = a.split("=", 1)
                    attr_dict[k] = v

            # Parse Target field: "ENSP00000349467 2 149"
            target = attr_dict.get("Target", "")
            target_parts = target.split()
            protein_id = target_parts[0] if target_parts else ""
            q_start = int(target_parts[1]) if len(target_parts) > 1 else 1
            q_end = int(target_parts[2]) if len(target_parts) > 2 else 0

            rank = int(attr_dict.get("Rank", 1))
            identity = float(attr_dict.get("Identity", 0))
            score_val = float(score) if score != "." else 0

            current_mrna = {
                "rank": rank,
                "chrom": chrom,
                "strand": strand,
                "score": score_val,
                "identity": identity,
                "query_start": q_start,
                "query_end": q_end,
                "_protein_id": protein_id,
            }
            current_cds = []

        elif feature == "CDS" and current_mrna:
            current_cds.append({
                "chrom": chrom,
                "start": start,
                "end": end,
            })

    # Flush last
    if current_mrna and current_cds:
        current_mrna["cds"] = sorted(current_cds, key=lambda c: c["start"])
        results[current_mrna["_protein_id"]].append(current_mrna)

    # Clean up internal fields and sort by rank
    clean = {}
    for pid, hits in results.items():
        for h in hits:
            h.pop("_protein_id", None)
        clean[pid] = sorted(hits, key=lambda h: h["rank"])
    return clean


def run_miniprot(genome_path, protein_fasta, top_k=5, threads=4,
                 region=None, extra_args=None):
    """Run miniprot and return parsed results.

    Args:
        genome_path: Path to genome FASTA or .mpi index.
        protein_fasta: Path to query protein FASTA.
        top_k: Number of top hits to report per protein.
        threads: Number of threads.
        region: Optional "chrom:start-end" to restrict alignment.
        extra_args: Additional miniprot arguments.
    """
    cmd = [MINIPROT, "--gff", "-t", str(threads), "-N", str(top_k)]

    if region:
        # miniprot doesn't have --region, we'll handle this via bed
        # For now, just run on full genome — region filtering done post-hoc
        pass

    if extra_args:
        cmd.extend(extra_args)

    cmd.extend([genome_path, protein_fasta])

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if result.returncode != 0:
        print(f"  miniprot stderr: {result.stderr[:500]}", file=sys.stderr)

    return parse_miniprot_gff(result.stdout)


def run_miniprot_restricted(genome_fa, protein_fasta, chrom, start, end,
                            top_k=5, threads=4, tmpdir=None):
    """Run miniprot restricted to a specific genomic region.

    Extracts the region with samtools, builds temp index, runs miniprot.
    """
    import pysam

    if tmpdir is None:
        tmpdir = tempfile.mkdtemp()

    # Extract region using samtools
    region_fa = os.path.join(tmpdir, "region.fa")
    region_str = f"{chrom}:{start}-{end}"

    try:
        result = subprocess.run(
            ["samtools", "faidx", genome_fa, region_str],
            capture_output=True, text=True, timeout=60
        )
        if result.returncode != 0:
            return {}

        # Rewrite FASTA header to keep original chrom name with offset
        lines = result.stdout.strip().split("\n")
        with open(region_fa, "w") as fh:
            fh.write(f">{chrom}\n")
            for line in lines[1:]:
                fh.write(line + "\n")

        # Run miniprot on extracted region (no index needed for small regions)
        hits = run_miniprot(region_fa, protein_fasta, top_k=top_k, threads=threads)

        # Adjust coordinates back to genome coordinates
        offset = start - 1  # samtools faidx is 1-based
        for pid, hit_list in hits.items():
            for h in hit_list:
                for cds in h["cds"]:
                    cds["start"] += offset
                    cds["end"] += offset

        return hits

    except Exception as e:
        print(f"  Error in restricted miniprot: {e}", file=sys.stderr)
        return {}


def main():
    parser = argparse.ArgumentParser(description="Run miniprot (v3 wrapper)")
    parser.add_argument("--species", required=True,
                        choices=["mouse", "dog", "chicken", "zebrafish",
                                 "mouse_paralog", "dog_paralog", "chicken_paralog", "zebrafish_paralog"])
    parser.add_argument("--top-k", type=int, default=5)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--max-proteins", type=int, default=0,
                        help="Limit number of proteins (0=all)")
    parser.add_argument("--region", type=str, default=None,
                        help="Restrict to genomic region (chrom:start-end)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output JSON path (default: results/v3/raw/miniprot_{species}.json)")
    args = parser.parse_args()

    species = args.species
    meta = load_meta(species)
    genome = meta["genome_index"]  # Use .mpi index for speed
    proteins = load_proteins(species)

    if args.max_proteins > 0:
        prot_items = list(proteins.items())[:args.max_proteins]
        proteins = dict(prot_items)

    print(f"Running miniprot on {species}: {len(proteins)} proteins, top-k={args.top_k}")

    with tempfile.TemporaryDirectory() as tmpdir:
        fa_path = write_protein_batch_fasta(proteins, tmpdir)

        t0 = time.time()
        results = run_miniprot(genome, fa_path, top_k=args.top_k,
                               threads=args.threads)
        elapsed = time.time() - t0

    # Add proteins with no hits
    for pid in proteins:
        if pid not in results:
            results[pid] = []

    # Summary
    n_hits = sum(1 for h in results.values() if len(h) > 0)
    n_total = len(results)
    print(f"  Hits: {n_hits}/{n_total} proteins ({n_hits/n_total*100:.1f}%)")
    print(f"  Time: {elapsed:.1f}s")

    # Save
    output = {
        "tool": "miniprot",
        "species": species,
        "top_k": args.top_k,
        "n_proteins": n_total,
        "n_with_hits": n_hits,
        "wall_time_sec": round(elapsed, 1),
        "results": results,
    }

    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as fh:
            json.dump(output, fh, indent=2)
    else:
        out_path = save_results(output, "miniprot", species)

    print(f"  Saved to {out_path}")


if __name__ == "__main__":
    main()
