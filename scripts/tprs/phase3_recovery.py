"""TPRS Phase 3 — Co-ortholog guided recovery for blind-miss proteins.

For proteins where NO A-tool finds the primary locus, but a co-ortholog
IS found by A-tools:
1. Run miniprot C-stage on the co-ortholog candidate region → get CDS
2. Extract CDS from genome, translate to protein sequence
3. Use MMseqs2 easy-search with this protein against the genome
4. The primary ortholog should be a top hit (higher identity)

Output:
  results/tprs/recovery/<protein_id>.gff         — Phase 3 recovery hits
  results/tprs/recovery/<protein_id>.cds.fa     — extracted CDS protein
  results/tprs/recovery/manifest_<species>.tsv  — summary
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
from rescue_paths import GENOMES, MINIPROT, PROJECT_ROOT, PROTEOMES, V3_DATA  # noqa: E402

REPO = PROJECT_ROOT
BED_DIR = REPO / "results/tprs/stage0_bed"
EXP3_BED = REPO / "results/exp3/stage0_bed"
GT_DIR = REPO / "results/exp2/gt"
DATA = V3_DATA
DBROOT = REPO / "results/exp3/db"
RECOVERY = REPO / "results/tprs/recovery"

A_TOOLS = ["miniprot_relaxed", "miniprot", "genblasta", "mmseqs", "last",
           "diamond", "spaln3", "genblastg", "tblastn"]
SPECIES = ["mouse", "dog", "chicken", "zebrafish"]
SPECIES_PROTEOME = {
    "mouse": "mouse_canonical.fa",
    "dog": "dog_canonical.fa",
    "chicken": "chicken_canonical.fa",
    "zebrafish": "zebrafish_canonical.fa",
}



def load_bed(path: Path) -> dict[str, list[tuple]]:
    out: dict[str, list[tuple]] = defaultdict(list)
    if not path.exists(): return out
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"): continue
            f = line.strip().split("\t")
            out[f[3]].append((f[0], int(f[1]), int(f[2]), f[5],
                             float(f[4]) if len(f) > 4 else 0.0))
    return out


def load_gt_json(sp: str) -> dict:
    path = DATA / f"{sp}_paralog" / "ground_truth.json"
    with open(path) as f: return json.load(f)


def load_gt_tsv(sp: str) -> dict[str, list[tuple]]:
    path = GT_DIR / f"{sp}.paralog_broad.tsv"
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
    if a[0] != b[0] or a[3] != b[3]: return False
    inter = max(0, min(a[2], b[2]) - max(a[1], b[1]))
    if inter == 0: return False
    return inter / (a[2] - a[1]) >= 0.5 and inter / (b[2] - b[1]) >= 0.5


def extract_region(genome: Path, chrom: str, start: int, end: int) -> str:
    region = f"{chrom}:{max(0, start)}-{end}"
    proc = subprocess.run(["samtools", "faidx", str(genome), region],
                         capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"samtools faidx failed: {proc.stderr[:200]}")
    return proc.stdout


def parse_miniprot_gff(gff_text: str, region_chrom: str, region_offset: int) -> list[dict]:
    """Parse miniprot GFF, return list of hits with CDS coordinates in genome space."""
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
        if len(fields) < 9: continue
        feat = fields[2]
        if feat == "mRNA":
            if cur is not None and cur_cds:
                cur["cds"] = sorted(cur_cds, key=lambda x: x["start"])
                hits.append(cur)
            cur_cds = []
            attrs = dict(p.split("=", 1) for p in fields[8].split(";") if "=" in p)
            target = attrs.get("Target", "").split()
            cur = {
                "chrom": region_chrom,
                "strand": fields[6],
                "score": float(fields[5]) if fields[5] != "." else 0.0,
                "identity": float(attrs.get("Identity", 0.0)),
                "query_start": int(target[1]) if len(target) > 1 else 0,
                "query_end": int(target[2]) if len(target) > 2 else 0,
                "rank": int(attrs.get("Rank", 1)),
            }
        elif feat == "CDS" and cur is not None:
            cur_cds.append({
                "chrom": region_chrom,
                "start": int(fields[3]) - 1 + region_offset,
                "end": int(fields[4]) + region_offset,
                "strand": cur.get("strand", "."),
                "phase": int(fields[7]) if fields[7] != "." else 0,
            })
    if cur is not None and cur_cds:
        cur["cds"] = sorted(cur_cds, key=lambda x: x["start"])
        hits.append(cur)
    return hits


def extract_cds_sequence(genome: Path, hit: dict) -> str:
    """Extract CDS DNA sequence from genome coordinates, handling GFF phase."""
    cds = hit.get("cds", [])
    if not cds:
        return ""
    strand = hit.get("strand", "+")
    # Sort by start coordinate (genomic order)
    cds_sorted = sorted(cds, key=lambda x: x["start"])

    # Extract plus-strand sequence for each exon
    exon_seqs = []
    for exon in cds_sorted:
        seq = extract_region(genome, exon["chrom"], exon["start"], exon["end"])
        seq = "".join(seq.splitlines()[1:]).upper()
        phase = exon.get("phase", 0)
        exon_seqs.append((seq, phase, exon))

    if strand == "-":
        # Minus strand: transcription goes high→low.
        # Process exons in reverse genomic order, revcomp each, trim phase.
        trimmed = []
        for seq, phase, exon in reversed(exon_seqs):
            # Reverse complement
            trans = str.maketrans("ACGT", "TGCA")
            rc = seq.translate(trans)[::-1]
            # Trim phase bases from 5' end (beginning of rc)
            if phase > 0 and len(rc) > phase:
                rc = rc[phase:]
            trimmed.append(rc)
        return "".join(trimmed)
    else:
        # Plus strand: trim phase from 5' end of each exon
        trimmed = []
        for seq, phase, exon in exon_seqs:
            if phase > 0 and len(seq) > phase:
                seq = seq[phase:]
            trimmed.append(seq)
        return "".join(trimmed)


def translate(cds_dna: str) -> str:
    """Translate DNA to protein (standard code)."""
    gencode = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein = []
    for i in range(0, len(cds_dna) - 2, 3):
        codon = cds_dna[i:i+3]
        aa = gencode.get(codon, 'X')
        if aa == '*': break
        protein.append(aa)
    return "".join(protein)


def search_with_mmseqs(protein_fasta: str, sp: str, tmp_dir: Path,
                        timeout: int = 600) -> str:
    """Run MMseqs2 easy-search with the recovered protein against genome DB."""
    db = DBROOT / "mmseqs2" / sp / "genomedb"
    out = tmp_dir / "result.tsv"
    cmd = [
        "mmseqs", "easy-search",
        str(tmp_dir / "query.fa"),
        str(db),
        str(out),
        str(tmp_dir / "mmseqs_tmp"),
        "-s", "1.0", "--search-type", "2",
        "--threads", "4", "--max-seqs", "200", "--format-mode", "0",
    ]
    # Write query
    (tmp_dir / "query.fa").write_text(protein_fasta)
    (tmp_dir / "mmseqs_tmp").mkdir(exist_ok=True)

    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        raise RuntimeError(f"mmseqs failed (rc={proc.returncode}): "
                          f"{proc.stderr[:200]}")
    log = logging.getLogger("phase3")
    n_hits = proc.stdout.count("\n") if proc.stdout else 0
    log.debug("mmseqs %s: %d hits in %.0fs", sp, n_hits, 0)
    if out.exists():
        return out.read_text()
    return ""


def mmseqs_hits_to_bed(mmseqs_output: str, min_identity: float = 0.3) -> list[tuple]:
    """Parse MMseqs2 output to candidate loci.  MMseqs2 fident is 0–1 fraction."""
    candidates = []
    for line in mmseqs_output.strip().splitlines():
        if not line.strip(): continue
        f = line.split("\t")
        if len(f) < 12: continue
        # MMseqs2 format: query target fident alnlen mismatch gapopen qstart qend tstart tend evalue bitscore
        pident = float(f[2])
        if pident < min_identity: continue
        tstart = int(f[8]); tend = int(f[9])
        candidates.append((
            f[1],  # target (chromosome)
            min(tstart, tend) - 1,
            max(tstart, tend),
            "+" if tstart < tend else "-",
        ))
    return candidates


def _any_overlap(a, b, lenient_strand: bool = True) -> bool:
    """Check if intervals overlap at all (1bp sufficient), with optional strand leniency."""
    if a[0] != b[0]:
        return False
    if a[3] == b[3] or lenient_strand:
        inter = max(0, min(a[2], b[2]) - max(a[1], b[1]))
        return inter > 0
    return False


def find_primary_locus(candidates: list[tuple], gt_entries: list[tuple],
                       lenient_strand: bool = True) -> bool:
    """Check if any candidate overlaps the ground-truth primary locus.
    Uses any-overlap (not reciprocal 50%) since MMseqs2 tblastn hits are
    often short (100-200bp) compared to GT loci (7-32kb)."""
    for c in candidates:
        for g in gt_entries:
            if _any_overlap(c, g, lenient_strand):
                return True
    return False


def load_proteome_by_gene(sp: str) -> dict[str, str]:
    """Load species proteome, indexed by Ensembl gene ID from FASTA header."""
    fa_path = PROTEOMES / SPECIES_PROTEOME.get(sp, "")
    if not fa_path.exists():
        return {}
    gene_to_protein: dict[str, str] = {}
    cur_gene = None
    cur_seq: list[str] = []
    with open(fa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_gene and cur_seq:
                    gene_to_protein[cur_gene] = "".join(cur_seq)
                cur_seq = []
                # Extract gene:ENSxxx from header
                for token in line[1:].strip().split():
                    if token.startswith("gene:"):
                        cur_gene = token.split(":")[1].split(".")[0]
                        break
                else:
                    cur_gene = None
            elif cur_gene:
                cur_seq.append(line.strip())
        if cur_gene and cur_seq:
            gene_to_protein[cur_gene] = "".join(cur_seq)
    return gene_to_protein


def run_phase3_recovery(sp: str, max_proteins: int | None = None,
                        timeout: int = 120, query_ids: list[str] | None = None) -> dict:
    """Run Phase 3 recovery for blind-miss proteins with found co-orthologs.

    Uses the canonical proteome (indexed by Ensembl gene ID) to get the
    co-ortholog protein directly, avoiding error-prone CDS extraction from
    miniprot alignments.
    """
    log = logging.getLogger("phase3")
    genome = GENOMES / f"{sp}.fa"
    gt_json = load_gt_json(sp)
    gt_tsv = load_gt_tsv(sp)
    qids_all = sorted(set(gt_tsv))
    if query_ids:
        qids_all = [q for q in qids_all if q in set(query_ids)]

    # Load proteome indexed by gene ID
    gene2prot = load_proteome_by_gene(sp)
    log.info("phase3 %s: loaded %d proteins from proteome", sp, len(gene2prot))

    # Load all tool candidates
    all_cands: dict[str, list[tuple]] = {}
    for tool in A_TOOLS:
        if tool == "miniprot_relaxed":
            bed = BED_DIR / f"{tool}.{sp}.paralog.bed"
        else:
            bed = EXP3_BED / f"{tool}.{sp}.paralog.bed"
        if bed.exists(): all_cands[tool] = load_bed(bed)

    # Identify blind-miss proteins at K=200 (union of all tools)
    blind_miss = []
    for q in qids_all:
        primary_gt = gt_tsv.get(q, [None])[0]
        if primary_gt is None:
            blind_miss.append(q)
            continue
        found = False
        seen = set()
        for tool, cands in all_cands.items():
            for c in cands.get(q, [])[:200]:
                key = (c[0], c[1], c[2], c[3])
                if key in seen: continue
                seen.add(key)
                if recip_overlap_50(c, primary_gt):
                    found = True; break
            if found: break
        if not found:
            blind_miss.append(q)

    log.info("phase3 %s: %d blind-miss at K=200", sp, len(blind_miss))

    # For each blind-miss, find co-ortholog candidates
    recoverable = []
    for q in blind_miss:
        gt_entry = gt_json.get(q, {})
        paralogs = gt_entry.get("paralog_loci", [])
        if not paralogs:
            continue

        # Find any candidate on a paralog chromosome+strand
        coorth_candidates = []
        for pl in paralogs:
            for tool, cands in all_cands.items():
                for c in cands.get(q, [])[:200]:
                    if c[0] == pl["chrom"] and c[3] == pl["strand"]:
                        coorth_candidates.append((tool, c, pl))
                        break
                if coorth_candidates and coorth_candidates[-1][1][0] == pl["chrom"]:
                    break

        if coorth_candidates:
            recoverable.append((q, gt_entry, coorth_candidates))

    log.info("phase3 %s: %d recoverable (co-ortholog found)", sp, len(recoverable))
    if max_proteins:
        recoverable = recoverable[:max_proteins]

    RECOVERY.mkdir(parents=True, exist_ok=True)

    # Build query list: collect protein sequences for all recoverable queries
    query_list: list[dict] = []  # {query_id, protein_seq, chrom, start, end, gt_entries, gene, identity}
    for q, gt_entry, coorth_cands in recoverable:
        primary_gt = gt_tsv.get(q, [None])[0]
        if primary_gt is None:
            continue

        protein_seq = None
        best_pl = None
        best_cand = None
        for tool, coorth_cand, pl in coorth_cands:
            gene_id = pl.get("target_gene", "")
            prot = gene2prot.get(gene_id, "")
            if prot and len(prot) >= 10:
                protein_seq = prot
                best_pl = pl
                best_cand = coorth_cand
                break

        if protein_seq is None:
            query_list.append({"query_id": q, "protein_seq": None, "status": "no_proteome_protein"})
            continue

        chrom, start, end, strand = best_cand[:4]
        query_list.append({
            "query_id": q, "protein_seq": protein_seq,
            "coorth_chrom": chrom, "coorth_start": start, "coorth_end": end,
            "gt_entries": gt_tsv.get(q, []),
            "coorth_gene": best_pl.get("target_gene", ""),
            "coorth_identity": best_pl.get("identity", 0.0),
            "primary_chrom": primary_gt[0],
        })

    valid_queries = [x for x in query_list if x["protein_seq"] is not None]
    no_prot = [x for x in query_list if x["protein_seq"] is None]
    log.info("phase3 %s: %d valid proteins, %d no-proteome", sp,
             len(valid_queries), len(no_prot))

    if valid_queries:
        # Single batched MMseqs2 run with all co-ortholog proteins
        with tempfile.TemporaryDirectory() as td:
            scratch = Path(td)
            fasta_parts = []
            for x in valid_queries:
                fasta_parts.append(f">{x['query_id']}\n{x['protein_seq']}\n")
            multi_fasta = "".join(fasta_parts)
            log.info("phase3 %s: running batched MMseqs2 with %d queries (%d aa total)",
                     sp, len(valid_queries), sum(len(x["protein_seq"]) for x in valid_queries))

            try:
                mmseqs_out = search_with_mmseqs(multi_fasta, sp, scratch)
            except Exception as e:
                log.error("phase3 %s: batched MMseqs2 failed: %s", sp, e)
                mmseqs_out = ""

            # Group hits by query ID
            hits_by_query: dict[str, list[tuple]] = defaultdict(list)
            for line in mmseqs_out.strip().splitlines():
                if not line.strip():
                    continue
                f = line.strip().split("\t")
                if len(f) < 12:
                    continue
                pident = float(f[2])
                if pident < 0.3:
                    continue
                tstart = int(f[8])
                tend = int(f[9])
                hits_by_query[f[0]].append((
                    f[1],  # chrom
                    min(tstart, tend) - 1,
                    max(tstart, tend),
                    "+" if tstart < tend else "-",
                ))

            # Check each query for primary locus overlap
            for x in valid_queries:
                candidates = hits_by_query.get(x["query_id"], [])
                found = find_primary_locus(candidates, x["gt_entries"])

                out_cds = RECOVERY / f"{x['query_id']}.cds.fa"
                out_cds.write_text(f">{x['query_id']}\n{x['protein_seq']}\n")

                x["status"] = "recovered" if found else "not_found"
                x["mmseqs_hits"] = len(candidates)
                x["found"] = found

    # Build result list
    results = []
    for x in query_list:
        r = {"query_id": x["query_id"], "status": x.get("status", "no_proteome_protein")}
        if x["protein_seq"] is not None:
            r.update({
                "coorth_chrom": x["coorth_chrom"],
                "coorth_start": x["coorth_start"],
                "coorth_end": x["coorth_end"],
                "primary_chrom": x["primary_chrom"],
                "protein_length": len(x["protein_seq"]),
                "mmseqs_hits": x.get("mmseqs_hits", 0),
                "found": x.get("found", False),
                "coorth_gene": x["coorth_gene"],
                "coorth_identity": x["coorth_identity"],
            })
        results.append(r)

    # Summary
    n_found = sum(1 for r in results if r.get("found"))
    n_attempted = len(valid_queries)
    log.info("phase3 %s DONE: %d/%d recovered (%.1f%%)",
             sp, n_found, n_attempted,
             100 * n_found / max(n_attempted, 1))

    # Write manifest
    manifest_path = RECOVERY / f"manifest_{sp}.tsv"
    cols = ["query_id", "status", "coorth_chrom", "coorth_start", "coorth_end",
            "primary_chrom", "protein_length", "mmseqs_hits",
            "found", "coorth_gene", "coorth_identity"]
    with open(manifest_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in results:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    return {"species": sp, "n_blind": len(blind_miss),
            "n_recoverable": len(recoverable), "n_attempted": n_attempted,
            "n_recovered": n_found}


def main() -> int:
    ap = argparse.ArgumentParser(description="TPRS Phase 3: Co-ortholog guided recovery")
    ap.add_argument("--species", nargs="*", default=SPECIES)
    ap.add_argument("--max-proteins", type=int, default=None,
                    help="Max proteins to attempt recovery on (for quick tests)")
    ap.add_argument("--timeout", type=int, default=120)
    ap.add_argument("--query-ids", nargs="*", default=None,
                    help="Only process these specific query IDs")
    args = ap.parse_args()

    RECOVERY.mkdir(parents=True, exist_ok=True)
    log_path = RECOVERY / "phase3.log"
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] TPRS %(message)s",
        handlers=[logging.StreamHandler(sys.stdout),
                  logging.FileHandler(log_path, mode="a")],
    )
    log = logging.getLogger("phase3")

    if not MINIPROT:
        log.error("miniprot not found")
        return 1

    totals = {"n_blind": 0, "n_recoverable": 0, "n_attempted": 0, "n_recovered": 0}
    for sp in args.species:
        log.info("begin phase3 %s", sp)
        result = run_phase3_recovery(sp, max_proteins=args.max_proteins,
                                     timeout=args.timeout,
                                     query_ids=args.query_ids)
        for k in totals:
            totals[k] += result[k]
        log.info("done phase3 %s: %d/%d recovered", sp,
                 result["n_recovered"], result["n_attempted"])

    log.info("PHASE3 ALL DONE: %d/%d recovered (%.1f%%)",
             totals["n_recovered"], totals["n_attempted"],
             100 * totals["n_recovered"] / max(totals["n_attempted"], 1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
