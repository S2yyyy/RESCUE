"""Phase 5 Step 1: Build candidate region → mouse protein mapping.

For each candidate region from Phase 4 training data:
1. Intersect with CDS BED to find overlapping genes
2. Look up the gene's protein in the mouse proteome
3. For primary-locus candidates: also try Ensembl Compara ortholog
4. Fall back to human query if no protein found

Output: candidate_protein_map.json + coverage report.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import (  # noqa: E402
    GENE_STRUCTURES,
    ORTHOLOGS,
    PROJECT_ROOT,
    PROTEOMES,
    V3_DATA,
)

REPO = PROJECT_ROOT
DATA = V3_DATA
GENESTRUCT = GENE_STRUCTURES
OUTDIR = REPO / "results/tprs/phase5"
PHASE4_DATA = REPO / "results/tprs/router/training_data_mouse.tsv"

SPECIES_PROTEOME = {
    "mouse": "mouse_canonical.fa",
    "dog": "dog_canonical.fa",
    "chicken": "chicken_canonical.fa",
    "zebrafish": "zebrafish_canonical.fa",
}


def load_proteome_by_gene(sp: str) -> dict[str, str]:
    """Load species proteome indexed by Ensembl gene ID."""
    path = PROTEOMES / SPECIES_PROTEOME.get(sp, "")
    if not path.exists():
        return {}
    gene2prot = {}
    with open(path) as fh:
        cur_gene, cur_seq = None, []
        for line in fh:
            if line.startswith(">"):
                if cur_gene and cur_seq:
                    gene2prot[cur_gene] = "".join(cur_seq)
                cur_seq = []
                header = line[1:].strip()
                gene_part = [p for p in header.split() if p.startswith("gene:")]
                cur_gene = gene_part[0].split(":")[1].split(".")[0] if gene_part else None
            elif cur_gene:
                cur_seq.append(line.strip())
        if cur_gene and cur_seq:
            gene2prot[cur_gene] = "".join(cur_seq)
    return gene2prot


def load_cds_bed(sp: str) -> dict[str, dict]:
    """Load CDS BED → {gene_id: {chrom, exons, gene_range}}."""
    cds = {}
    bed_path = GENESTRUCT / f"{sp}_CDS.bed"
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip():
                continue
            f = line.strip().split("\t")
            gene_id = f[3].rsplit("_", 1)[0] if "_" in f[3] else f[3]
            if gene_id not in cds:
                cds[gene_id] = {"chrom": f[0], "exons": [], "strand": f[5]}
            cds[gene_id]["exons"].append({"start": int(f[1]), "end": int(f[2])})
    # Add gene_range for quick overlap checks
    for gene_id, info in cds.items():
        starts = [e["start"] for e in info["exons"]]
        ends = [e["end"] for e in info["exons"]]
        info["gene_range"] = (min(starts), max(ends))
    return cds


def load_ortholog_map(sp: str) -> dict[str, str]:
    """Load human ENSP → target-species gene one2one ortholog mapping."""
    species_map = {
        "mouse": "mus_musculus",
        "dog": "canis_lupus_familiaris",
        "chicken": "gallus_gallus",
        "zebrafish": "danio_rerio",
    }
    sp_name = species_map.get(sp, "mus_musculus")
    path = ORTHOLOGS / f"human_vs_{sp_name}_one2one.tsv"
    if not path.exists():
        return {}
    human2gene = {}
    with open(path) as fh:
        header = next(fh).strip().split("\t")
        gene_idx = header.index("gene_stable_id")
        prot_idx = header.index("homology_protein_stable_id")
        for line in fh:
            if not line.strip():
                continue
            f = line.strip().split("\t")
            human_prot = f[prot_idx]
            target_gene = f[gene_idx]
            if human_prot not in human2gene:
                human2gene[human_prot] = target_gene
    return human2gene


def build_gene_index(cds_bed: dict) -> dict:
    """Build a simple interval index: (chrom, start, end) → list of gene_ids.

    Uses gene_range from CDS BED for coarse overlap detection.
    """
    index = defaultdict(list)
    for gene_id, info in cds_bed.items():
        chrom = info["chrom"]
        gstart, gend = info["gene_range"]
        # Use a coarse binning: 1Mb bins
        for bin_start in range(gstart // 1_000_000 * 1_000_000,
                               gend + 1_000_000, 1_000_000):
            index[(chrom, bin_start)].append((gene_id, gstart, gend))
    return index


def find_best_gene(region_chrom: str, region_start: int, region_end: int,
                   cds_bed: dict, gene_index: dict,
                   primary_gene: str = "") -> tuple[str, str]:
    """Find the best gene overlapping a candidate region.

    Returns (gene_id, protein_source) where protein_source is
    "target_gene_overlap", "primary_ortholog", or "" (none found).

    Priority: if region overlaps primary_gene → use primary_gene.
    Otherwise: use the gene with most overlap (largest CDS in region).
    """
    best_gene = ""
    best_overlap = 0

    # Check genes in nearby bins
    for bin_start in range(region_start // 1_000_000 * 1_000_000,
                           region_end + 1_000_000, 1_000_000):
        for gene_id, gstart, gend in gene_index.get((region_chrom, bin_start), []):
            overlap = max(0, min(region_end, gend) - max(region_start, gstart))
            if overlap > best_overlap:
                best_overlap = overlap
                best_gene = gene_id
            # If gene overlaps primary_gene, prioritize it
            if primary_gene and gene_id == primary_gene and overlap > 0:
                return gene_id, "primary_ortholog"

    if best_gene:
        return best_gene, "target_gene_overlap"
    return "", ""


def main():
    sp = "mouse"

    print("Loading data...")
    cds_bed = load_cds_bed(sp)
    gene2prot = load_proteome_by_gene(sp)
    human2mouse = load_ortholog_map(sp)
    df = pd.read_csv(PHASE4_DATA, sep="\t")

    print(f"CDS genes: {len(cds_bed)}, proteome proteins: {len(gene2prot)}, "
          f"ortholog pairs: {len(human2mouse)}")
    print(f"Training data: {len(df)} rows, {df['query_id'].nunique()} queries")

    # Load GT for primary gene mapping
    gt_path = DATA / f"{sp}_paralog" / "ground_truth.json"
    with open(gt_path) as f:
        gt_data = json.load(f)

    # Build gene index
    gene_index = build_gene_index(cds_bed)

    # Process each unique region
    # Group by (query_id, chrom, start, end) to get unique regions
    region_cols = ["query_id", "chrom", "start", "end"]
    regions = df[region_cols].drop_duplicates()
    print(f"Unique regions: {len(regions)}")

    mapping = {}
    stats = {
        "total": len(regions),
        "gene_found": 0,
        "protein_found": 0,
        "ortholog_used": 0,
        "fallback_human": 0,
    }

    for _, row in regions.iterrows():
        qid = row["query_id"]
        chrom = str(row["chrom"])
        rstart = int(row["start"])
        rend = int(row["end"])
        region_key = f"{qid}|{chrom}|{rstart}|{rend}"

        # Get primary gene from GT
        gt = gt_data.get(qid, {})
        primary_gene = gt.get("target_gene", "")

        # Find best overlapping gene
        best_gene, source = find_best_gene(
            chrom, rstart, rend, cds_bed, gene_index, primary_gene)

        # If no gene from CDS overlap, try Ensembl ortholog (for primary locus)
        if not best_gene:
            ortholog_gene = human2mouse.get(qid, "")
            if ortholog_gene and ortholog_gene in cds_bed:
                og_info = cds_bed[ortholog_gene]
                og_start, og_end = og_info["gene_range"]
                # Check if this region overlaps the ortholog gene
                if max(0, min(rend, og_end) - max(rstart, og_start)) > 0:
                    best_gene = ortholog_gene
                    source = "primary_ortholog"
                    stats["ortholog_used"] += 1

        # Look up protein
        protein_seq = ""
        if best_gene:
            stats["gene_found"] += 1
            protein_seq = gene2prot.get(best_gene, "")
            if protein_seq:
                stats["protein_found"] += 1
            else:
                # Try ortholog as fallback
                ortholog_gene = human2mouse.get(qid, "")
                if ortholog_gene and ortholog_gene != best_gene:
                    protein_seq = gene2prot.get(ortholog_gene, "")
                    if protein_seq:
                        best_gene = ortholog_gene
                        source = "primary_ortholog_fallback"
                        stats["protein_found"] += 1

        if not protein_seq:
            stats["fallback_human"] += 1

        mapping[region_key] = {
            "query_id": qid,
            "chrom": chrom,
            "start": rstart,
            "end": rend,
            "best_gene": best_gene,
            "protein_source": source,
            "has_protein": len(protein_seq) > 0,
            "protein_length": len(protein_seq),
        }

    # Coverage report
    print(f"\n{'='*60}")
    print("COVERAGE REPORT")
    print(f"{'='*60}")
    print(f"Total unique regions:      {stats['total']}")
    print(f"Regions with gene overlap: {stats['gene_found']} "
          f"({100*stats['gene_found']/stats['total']:.1f}%)")
    print(f"Regions with protein:      {stats['protein_found']} "
          f"({100*stats['protein_found']/stats['total']:.1f}%)")
    print(f"Ortholog used:             {stats['ortholog_used']}")
    print(f"Fallback to human query:   {stats['fallback_human']} "
          f"({100*stats['fallback_human']/stats['total']:.1f}%)")

    # Per-query stats
    query_coverage = defaultdict(lambda: {"total": 0, "with_protein": 0})
    for key, info in mapping.items():
        qid = info["query_id"]
        query_coverage[qid]["total"] += 1
        if info["has_protein"]:
            query_coverage[qid]["with_protein"] += 1

    queries_with_coverage = sum(
        1 for qc in query_coverage.values()
        if qc["with_protein"] > 0)
    queries_full_coverage = sum(
        1 for qc in query_coverage.values()
        if qc["with_protein"] == qc["total"])
    print(f"\nQueries with ≥1 protein-mapped region: {queries_with_coverage}/{len(query_coverage)} "
          f"({100*queries_with_coverage/len(query_coverage):.1f}%)")
    print(f"Queries with ALL regions mapped: {queries_full_coverage}/{len(query_coverage)} "
          f"({100*queries_full_coverage/len(query_coverage):.1f}%)")

    # Save mapping
    OUTDIR.mkdir(parents=True, exist_ok=True)
    map_path = OUTDIR / "candidate_protein_map.json"
    with open(map_path, "w") as f:
        json.dump({"stats": stats, "mapping": mapping}, f, indent=2)
    print(f"\nMapping saved to {map_path}")

    # Also save a compact TSV for phase5_species_query.py
    tsv_path = OUTDIR / "candidate_protein_map.tsv"
    with open(tsv_path, "w") as f:
        f.write("region_key\tquery_id\tchrom\tstart\tend\tbest_gene\t"
                "protein_source\thas_protein\tprotein_length\n")
        for key, info in sorted(mapping.items()):
            f.write(f"{key}\t{info['query_id']}\t{info['chrom']}\t"
                    f"{info['start']}\t{info['end']}\t{info['best_gene']}\t"
                    f"{info['protein_source']}\t{info['has_protein']}\t"
                    f"{info['protein_length']}\n")
    print(f"TSV saved to {tsv_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
