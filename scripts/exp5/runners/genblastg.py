"""genBlastG forced-locus runner.

Calls genBlastG with -p genblastg -gff. The binary requires being run from
its own thecgs/ directory because it shells out to a relative `./formatdb`.
We invoke it via cwd=thecgs and absolute paths to query/target/out.

Output filename has a parameter-suffix (e.g. `<out>_1.1c_2.3_s1_0_16_1.gff`),
so we glob for the matching .gff.
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

from ..common import run_subproc
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from rescue_paths import TOOLS_ROOT  # noqa: E402

GENBLAST_DIR = Path(os.environ.get("GENBLAST_DIR", TOOLS_ROOT / "genblast/thecgs"))
GENBLAST_BIN = "./genblastG"


def run_one(protein_id, protein_seq, region_fa, region_chrom,
            region_offset, region_strand, scratch, options) -> dict | None:
    # genBlastG runs formatdb on the target FASTA — point it at the scratch
    # copy so artefacts stay out of the parent extraction dir.
    region_local = scratch / "region.fa"
    if region_fa.resolve() != region_local.resolve():
        shutil.copy(region_fa, region_local)
    prot_fa = scratch / f"{protein_id}.prot.fa"
    prot_fa.write_text(f">{protein_id}\n{protein_seq}\n")
    out_prefix = scratch / f"{protein_id}.gbg"

    cmd = [GENBLAST_BIN, "-p", "genblastg",
           "-q", str(prot_fa), "-t", str(region_local),
           "-gff", "-o", str(out_prefix)]
    rc, _, _ = run_subproc(cmd, timeout=options.get("timeout", 300),
                            cwd=str(GENBLAST_DIR))
    if rc != 0:
        return None
    gff_files = list(scratch.glob(f"{protein_id}.gbg*.gff"))
    if not gff_files:
        return None
    text = gff_files[0].read_text()
    return _parse(text, region_chrom, region_offset)


def _parse(gff_text: str, region_chrom: str, region_offset: int) -> dict | None:
    transcripts: list[dict] = []
    cur: dict | None = None
    cur_cds: list[dict] = []

    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        feat = f[2]
        if feat == "transcript":
            if cur and cur_cds:
                cur["cds"] = sorted(cur_cds, key=lambda c: c["start"])
                transcripts.append(cur)
            cur_cds = []
            try:
                sc = float(f[5]) if f[5] != "." else 0.0
            except ValueError:
                sc = 0.0
            cur = {
                "rank": 1, "chrom": region_chrom,
                "strand": f[6], "score": sc, "identity": 0,
                "query_start": 1, "query_end": 0,
            }
        elif feat == "coding_exon" and cur is not None:
            try:
                s = int(f[3]); e = int(f[4])
            except ValueError:
                continue
            cur_cds.append({
                "chrom": region_chrom,
                "start": min(s, e) + region_offset,
                "end":   max(s, e) + region_offset,
            })
    if cur and cur_cds:
        cur["cds"] = sorted(cur_cds, key=lambda c: c["start"])
        transcripts.append(cur)
    if not transcripts:
        return None
    transcripts.sort(key=lambda h: -h["score"])
    return transcripts[0]
