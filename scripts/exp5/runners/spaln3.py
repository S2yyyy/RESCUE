"""Spaln3 forced-locus runner. Calls `spaln -O0 -Q7 -M1`."""
from __future__ import annotations

from pathlib import Path

from ..common import run_subproc

SPALN = "spaln"


def run_one(protein_id, protein_seq, region_fa, region_chrom,
            region_offset, region_strand, scratch, options) -> dict | None:
    prot_fa = scratch / f"{protein_id}.prot.fa"
    prot_fa.write_text(f">{protein_id}\n{protein_seq}\n")
    cmd = [SPALN, "-O0", "-Q0", "-M1", str(region_fa), str(prot_fa)]
    rc, stdout, _ = run_subproc(cmd, timeout=options.get("timeout", 120))
    prot_fa.unlink(missing_ok=True)
    if rc != 0 or not stdout.strip():
        return None
    return _parse(stdout, region_chrom, region_offset)


def _parse(gff_text: str, region_chrom: str, region_offset: int) -> dict | None:
    cur = None
    cur_cds: list[dict] = []
    best = None
    best_score = -1e18

    def _flush():
        nonlocal best, best_score
        if cur is None or not cur_cds:
            return
        sc = cur.get("score", 0)
        if sc > best_score:
            best_score = sc
            cur["cds"] = [
                {"chrom": region_chrom,
                 "start": int(c["start"]) + region_offset,
                 "end": int(c["end"]) + region_offset}
                for c in sorted(cur_cds, key=lambda c: int(c["start"]))
            ]
            cur["chrom"] = region_chrom
            best = dict(cur)

    for line in gff_text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 8:
            continue
        feat = f[2].lower()
        if feat in ("gene", "mrna"):
            _flush()
            cur_cds = []
            attr = dict(a.split("=", 1) for a in f[8].split(";") if "=" in a)
            tgt = attr.get("Target", attr.get("Name", "")).split()
            qs = int(tgt[1]) if len(tgt) > 1 else 1
            qe = int(tgt[2]) if len(tgt) > 2 else 0
            try:
                sc = float(f[5]) if f[5] != "." else 0.0
            except ValueError:
                sc = 0.0
            cur = {
                "rank": 1, "strand": f[6], "score": sc,
                "identity": 0,
                "query_start": qs, "query_end": qe,
            }
        elif feat == "cds" and cur is not None:
            cur_cds.append({"start": f[3], "end": f[4]})
            # Spaln embeds Target with q_start/q_end on CDS lines too
            if len(f) >= 9:
                attr = dict(a.split("=", 1) for a in f[8].split(";") if "=" in a)
                tgt = attr.get("Target", "").split()
                if len(tgt) >= 3:
                    try:
                        qs, qe = int(tgt[1]), int(tgt[2])
                        if cur.get("query_end", 0) == 0:
                            cur["query_start"], cur["query_end"] = qs, qe
                        else:
                            cur["query_start"] = min(cur["query_start"], qs)
                            cur["query_end"] = max(cur["query_end"], qe)
                    except ValueError:
                        pass
    _flush()
    return best
