#!/usr/bin/env python3
"""V3 Common utilities for tool wrappers.

Standard output format for all tools:
  {
    protein_id: {
      "hits": [
        {
          "rank": 1,
          "chrom": "12",
          "strand": "+",
          "score": 614,
          "identity": 1.0,
          "query_start": 1,   # 1-based protein coords
          "query_end": 149,
          "query_len": 149,
          "cds": [
            {"chrom": "12", "start": 100168652, "end": 100168682},
            ...
          ]
        },
        ...  # up to top_k hits
      ]
    }
  }

Coordinates: 1-based, inclusive (GFF convention).
"""

import json
import os
import subprocess
import tempfile
from pathlib import Path

from rescue_paths import PROCESSED_DATA, V3_DATA, V3_RESULTS

PROCESSED = PROCESSED_DATA


def load_proteins(species):
    """Load protein sequences from v3 data."""
    fa_path = V3_DATA / species / "proteins.fa"
    seqs = {}
    cur_id, cur_seq = "", []
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_id:
                    seqs[cur_id] = "".join(cur_seq)
                cur_id = line[1:].strip().split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
    if cur_id:
        seqs[cur_id] = "".join(cur_seq)
    return seqs


def load_ground_truth(species):
    """Load ground truth from v3 data."""
    gt_path = V3_DATA / species / "ground_truth.json"
    with open(gt_path) as fh:
        return json.load(fh)


def load_meta(species):
    """Load species metadata."""
    meta_path = V3_DATA / species / "meta.json"
    with open(meta_path) as fh:
        return json.load(fh)


def save_results(results, tool, species):
    """Save tool results in standard format."""
    out_dir = V3_RESULTS / "raw"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{tool}_{species}.json"
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)
    return out_path


def write_single_protein_fasta(protein_id, sequence, tmpdir):
    """Write a single protein to a temp FASTA file."""
    fa_path = os.path.join(tmpdir, f"{protein_id}.fa")
    with open(fa_path, "w") as fh:
        fh.write(f">{protein_id}\n{sequence}\n")
    return fa_path


def write_protein_batch_fasta(proteins, tmpdir, batch_name="batch"):
    """Write multiple proteins to a FASTA file."""
    fa_path = os.path.join(tmpdir, f"{batch_name}.fa")
    with open(fa_path, "w") as fh:
        for pid, seq in proteins.items():
            fh.write(f">{pid}\n{seq}\n")
    return fa_path


def compute_query_coverage(hit, protein_len):
    """Compute fraction of protein covered by alignment."""
    if protein_len == 0:
        return 0.0
    covered = hit.get("query_end", 0) - hit.get("query_start", 1) + 1
    return covered / protein_len
