#!/usr/bin/env python3
"""Build ground-truth tables for 4 species x 3 regimes (12 TSVs).

Produces (default: $RESCUE_DATA_ROOT/gt/):
  <species>.one2one.tsv        - primary locus per one-to-one query
  <species>.paralog_strict.tsv - primary (highest-identity) locus per paralog query
  <species>.paralog_broad.tsv  - primary locus plus every acceptable
                                 co-ortholog locus from paralog_info.json

Schema: query_id  chrom  start  end  strand  gene_id
Coordinates are 0-based half-open; gene body span = min/max of CDS exons.

Inputs:
  data/v3/<species>/ground_truth.json
  data/v3/<species>_paralog/{ground_truth.json, paralog_info.json, pairs.tsv}
  data/processed/gene_structures/<species>_CDS.bed
  data/processed/orthologs/human_vs_<species>_paralog_rich.tsv

Standard library only.
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import GENE_STRUCTURES, ORTHOLOGS, V3_DATA, DATA_ROOT  # noqa: E402

SPECIES_TO_ENSEMBL = {
    "mouse": "mus_musculus",
    "dog": "canis_lupus_familiaris",
    "chicken": "gallus_gallus",
    "zebrafish": "danio_rerio",
}

EXPECTED_QID_COUNTS = {
    "one2one": {"mouse": 2000, "dog": 2000, "chicken": 2000, "zebrafish": 2000},
    "paralog_strict": {"mouse": 500, "dog": 500, "chicken": 488, "zebrafish": 500},
    "paralog_broad": {"mouse": 500, "dog": 500, "chicken": 488, "zebrafish": 500},
}

GT_HEADER = ["query_id", "chrom", "start", "end", "strand", "gene_id"]


@dataclass(frozen=True)
class PrimaryLocus:
    chrom: str
    start: int
    end: int
    strand: str
    gene_id: str


@dataclass(frozen=True)
class CDSBlock:
    chrom: str
    start: int
    end: int
    gene_id: str
    transcript_id: str
    strand: str


@dataclass(frozen=True)
class CanonicalLocus:
    chrom: str
    start: int
    end: int
    strand: str
    transcript_id: str
    cds_total: int


def load_primary_gt_json(path: Path) -> dict[str, PrimaryLocus]:
    with open(path) as fh:
        data = json.load(fh)
    out: dict[str, PrimaryLocus] = {}
    for protein_id, entry in data.items():
        exons = entry.get("exons", [])
        if not exons:
            continue
        starts = [int(e["start"]) for e in exons]
        ends = [int(e["end"]) for e in exons]
        out[protein_id] = PrimaryLocus(
            chrom=str(entry["chrom"]),
            start=min(starts),
            end=max(ends),
            strand=entry["strand"],
            gene_id=entry["target_gene"],
        )
    return out


def load_cds_bed(path: Path) -> dict[str, dict[str, list[CDSBlock]]]:
    out: dict[str, dict[str, list[CDSBlock]]] = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            out[f[3]][f[4]].append(CDSBlock(
                chrom=f[0], start=int(f[1]), end=int(f[2]),
                gene_id=f[3], transcript_id=f[4], strand=f[5],
            ))
    return out


def canonical_locus_for_gene(gene_to_tx: dict[str, list[CDSBlock]]) -> CanonicalLocus | None:
    """Longest total CDS; ties broken by lexicographically smallest tx ID."""
    best: CanonicalLocus | None = None
    for tx_id in sorted(gene_to_tx.keys()):
        blocks = gene_to_tx[tx_id]
        if not blocks:
            continue
        total = sum(b.end - b.start for b in blocks)
        cand = CanonicalLocus(
            chrom=blocks[0].chrom,
            start=min(b.start for b in blocks),
            end=max(b.end for b in blocks),
            strand=blocks[0].strand,
            transcript_id=tx_id,
            cds_total=total,
        )
        if best is None or cand.cds_total > best.cds_total:
            best = cand
    return best


def build_gene_canonical_map(cds: dict[str, dict[str, list[CDSBlock]]]) -> dict[str, CanonicalLocus]:
    return {g: c for g, tx in cds.items() if (c := canonical_locus_for_gene(tx)) is not None}


def load_pairs_query_ids(path: Path) -> set[str]:
    out: set[str] = set()
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            qid = row.get("human_protein") or row.get("protein_id")
            if qid:
                out.add(qid)
    return out


def load_homology_table(path: Path) -> dict[tuple[str, str], str]:
    """{(human_protein, target_gene) -> homology_type}."""
    out: dict[tuple[str, str], str] = {}
    if not path.exists():
        return out
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            human_p = row.get("protein_stable_id", "")
            target_g = row.get("homology_gene_stable_id", "")
            ht = row.get("homology_type", "")
            if human_p and target_g and ht:
                out.setdefault((human_p, target_g), ht)
    return out


def write_gt_tsv(path: Path, rows: list[tuple[str, str, int, int, str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(GT_HEADER)
        writer.writerows(rows)
    tmp.rename(path)
    return len(rows)


def build_primary_rows(gt: dict[str, PrimaryLocus]) -> list[tuple[str, str, int, int, str, str]]:
    rows = [
        (qid, locus.chrom, locus.start, locus.end, locus.strand, locus.gene_id)
        for qid, locus in gt.items()
    ]
    rows.sort(key=lambda r: (r[0], r[5], r[1], r[2], r[3]))
    return rows


def build_paralog_broad_rows(
    species: str,
    paralog_info: dict[str, dict],
    primary_strict: dict[str, PrimaryLocus],
    canonical_by_gene: dict[str, CanonicalLocus],
    homology_lookup: dict[tuple[str, str], str],
    audit_rows: list[dict],
    log: logging.Logger,
) -> list[tuple[str, str, int, int, str, str]]:
    rows: list[tuple[str, str, int, int, str, str]] = []
    seen: set[tuple[str, str]] = set()
    n_missing_coords = 0
    n_homology_warnings = 0

    for qid in sorted(paralog_info.keys()):
        info = paralog_info[qid]
        for t in info.get("all_targets", []):
            target_gene = t.get("target_gene")
            if not target_gene:
                continue
            key = (qid, target_gene)
            if key in seen:
                continue
            seen.add(key)
            is_best = bool(t.get("is_best", False))

            homology_type = homology_lookup.get((qid, target_gene), "")
            if homology_type and homology_type not in (
                "ortholog_one2many", "ortholog_many2many"
            ):
                log.error("broad: unexpected homology_type=%s for %s -> %s",
                          homology_type, qid, target_gene)
                n_homology_warnings += 1

            if is_best and qid in primary_strict:
                p = primary_strict[qid]
                rows.append((qid, p.chrom, p.start, p.end, p.strand, target_gene))
                coord_source = "paralog_ground_truth_json"
                chrom, start, end, strand = p.chrom, p.start, p.end, p.strand
            else:
                canon = canonical_by_gene.get(target_gene)
                if canon is None:
                    n_missing_coords += 1
                    audit_rows.append({
                        "species": species, "query_id": qid,
                        "human_gene": info.get("human_gene", ""),
                        "target_gene": target_gene, "is_best": is_best,
                        "identity": t.get("identity", ""),
                        "homology_type": homology_type,
                        "coordinate_source": "missing_in_cds_bed",
                        "chrom": t.get("chrom") or "", "start": "", "end": "",
                        "strand": "", "issue": "missing_coords_in_cds_bed",
                    })
                    continue
                rows.append((qid, canon.chrom, canon.start, canon.end,
                             canon.strand, target_gene))
                coord_source = "cds_bed_canonical_transcript"
                chrom, start, end, strand = canon.chrom, canon.start, canon.end, canon.strand

            audit_rows.append({
                "species": species, "query_id": qid,
                "human_gene": info.get("human_gene", ""),
                "target_gene": target_gene, "is_best": is_best,
                "identity": t.get("identity", ""),
                "homology_type": homology_type,
                "coordinate_source": coord_source,
                "chrom": chrom, "start": start, "end": end,
                "strand": strand, "issue": "",
            })

    rows.sort(key=lambda r: (r[0], r[5], r[1], r[2], r[3]))
    log.info("broad species=%s rows=%d missing_coords=%d homology_warnings=%d",
             species, len(rows), n_missing_coords, n_homology_warnings)
    return rows


def build_species(species: str, out_dir: Path, audit_rows: list[dict],
                  log: logging.Logger) -> list[dict]:
    log.info("=== species=%s ===", species)
    v3_dir = V3_DATA / species
    paralog_dir = V3_DATA / f"{species}_paralog"

    have_one2one = (v3_dir / "ground_truth.json").exists()
    have_paralog = (paralog_dir / "ground_truth.json").exists() and \
        (paralog_dir / "paralog_info.json").exists()
    if not have_one2one:
        log.warning("%s: no one2one query set at %s, regime skipped",
                    species, v3_dir)
    if not have_paralog:
        log.warning("%s: no paralog query set at %s, regimes skipped",
                    species, paralog_dir)
    if not have_one2one and not have_paralog:
        return []

    primary_one2one = (load_primary_gt_json(v3_dir / "ground_truth.json")
                       if have_one2one else {})
    primary_paralog = (load_primary_gt_json(paralog_dir / "ground_truth.json")
                       if have_paralog else {})
    paralog_info = (json.loads((paralog_dir / "paralog_info.json").read_text())
                    if have_paralog else {})
    pairs_qids = (load_pairs_query_ids(paralog_dir / "pairs.tsv")
                  if (paralog_dir / "pairs.tsv").exists() else set())

    cds = load_cds_bed(GENE_STRUCTURES / f"{species}_CDS.bed")
    canonical_by_gene = build_gene_canonical_map(cds)
    homology_lookup = load_homology_table(
        ORTHOLOGS / f"human_vs_{SPECIES_TO_ENSEMBL[species]}_paralog_rich.tsv")

    log.info("one2one queries=%d, paralog queries=%d, CDS genes=%d",
             len(primary_one2one), len(primary_paralog), len(cds))

    rows_one2one = build_primary_rows(primary_one2one)
    rows_strict = build_primary_rows(primary_paralog)
    rows_broad = (build_paralog_broad_rows(
        species, paralog_info, primary_paralog, canonical_by_gene,
        homology_lookup, audit_rows, log) if have_paralog else [])

    out_dir.mkdir(parents=True, exist_ok=True)
    n1 = n2 = n3 = 0
    if have_one2one:
        n1 = write_gt_tsv(out_dir / f"{species}.one2one.tsv", rows_one2one)
    if have_paralog:
        n2 = write_gt_tsv(out_dir / f"{species}.paralog_strict.tsv", rows_strict)
        n3 = write_gt_tsv(out_dir / f"{species}.paralog_broad.tsv", rows_broad)
    log.info("wrote %s.one2one rows=%d, strict rows=%d, broad rows=%d",
             species, n1, n2, n3)

    summaries = []
    regimes = []
    if have_one2one:
        regimes.append(("one2one", rows_one2one))
    if have_paralog:
        regimes += [("paralog_strict", rows_strict),
                    ("paralog_broad", rows_broad)]
    for regime, rows in regimes:
        n_input = len(primary_one2one) if regime == "one2one" else len(primary_paralog)
        n_qids = len({r[0] for r in rows})
        seen_pairs: set[tuple[str, str]] = set()
        n_dup = sum(1 for r in rows
                    if (r[0], r[5]) in seen_pairs or seen_pairs.add((r[0], r[5])))
        summaries.append({
            "species": species,
            "regime": regime,
            "path": str(out_dir / f"{species}.{regime}.tsv"),
            "n_query_ids": n_qids,
            "n_rows": len(rows),
            "expected_query_ids": EXPECTED_QID_COUNTS[regime][species],
            "n_missing_query_ids": (
                len(pairs_qids - {r[0] for r in rows}) if regime != "one2one" else 0),
            "n_duplicate_query_gene_rows": n_dup,
            # OK when every input query is represented exactly once; the
            # expected_query_ids column documents the full-benchmark counts
            # (smaller subsets, e.g. the test dataset, are fine).
            "status": "OK" if n_qids == n_input and n_dup == 0 else "CHECK",
        })
    return summaries


def write_manifest(rows: list[dict], path: Path) -> None:
    cols = ["species", "regime", "path", "n_query_ids", "n_rows",
            "expected_query_ids", "n_missing_query_ids",
            "n_duplicate_query_gene_rows", "status"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_audit(rows: list[dict], path: Path) -> None:
    cols = ["species", "query_id", "human_gene", "target_gene", "is_best",
            "identity", "homology_type", "coordinate_source", "chrom",
            "start", "end", "strand", "issue"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    p = argparse.ArgumentParser(description="Build ground-truth tables")
    p.add_argument("--out-dir", default=str(DATA_ROOT / "gt"))
    p.add_argument("--species", nargs="+",
                   default=["mouse", "dog", "chicken", "zebrafish"])
    p.add_argument("--log-level", default="INFO")
    args = p.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    handlers: list[logging.Handler] = [
        logging.StreamHandler(sys.stderr),
        logging.FileHandler(out_dir / "build_gt_tables.log", mode="w"),
    ]
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=handlers, force=True,
    )
    log = logging.getLogger("build_gt_tables")

    manifest_rows: list[dict] = []
    audit_rows: list[dict] = []
    for species in args.species:
        try:
            manifest_rows.extend(build_species(species, out_dir, audit_rows, log))
        except Exception as exc:
            log.exception("species %s failed: %s", species, exc)
            for regime in ("one2one", "paralog_strict", "paralog_broad"):
                manifest_rows.append({"species": species, "regime": regime,
                                      "status": f"FAIL_EXCEPTION: {exc}"})

    write_manifest(manifest_rows, out_dir / "manifest.tsv")
    write_audit(audit_rows, out_dir / "broad_target_audit.tsv")
    bad = [r for r in manifest_rows if r.get("status") != "OK"]
    if bad:
        log.warning("%d manifest rows not OK", len(bad))
        for r in bad:
            log.warning("  %s/%s status=%s", r["species"], r["regime"], r["status"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
