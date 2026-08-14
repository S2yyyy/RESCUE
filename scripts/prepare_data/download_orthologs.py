#!/usr/bin/env python3
"""Download Ensembl Compara homology tables (release 113) and build the
per-species ortholog TSVs used by the RESCUE benchmark.

Source: the pre-computed Compara homology flat dumps at
  https://ftp.ensembl.org/pub/release-113/tsv/ensembl-compara/homologies/homo_sapiens/

Outputs (under $RESCUE_ORTHOLOG_DIR, default data/processed/orthologs/):
  human_vs_<species>_one2one.tsv      - homology_type == ortholog_one2one
                                        (both directions kept)
  human_vs_<species>_paralog_rich.tsv - homology_type in
                                        {ortholog_one2many, ortholog_many2many},
                                        human side only, with an added
                                        n_target_loci column (distinct target
                                        genes per human gene)

Standard library only.
"""
from __future__ import annotations

import csv
import gzip
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from rescue_paths import DATA_ROOT, ORTHOLOGS  # noqa: E402

RELEASE = 113
FTP_DIR = (
    f"https://ftp.ensembl.org/pub/release-{RELEASE}/tsv/"
    "ensembl-compara/homologies/homo_sapiens/"
)
RAW_DIR = DATA_ROOT / "raw" / f"release-{RELEASE}"

SPECIES = {
    "mouse": "mus_musculus",
    "dog": "canis_lupus_familiaris",
    "chicken": "gallus_gallus",
    "zebrafish": "danio_rerio",
}

COLUMNS = [
    "gene_stable_id", "protein_stable_id", "species", "identity",
    "homology_type", "homology_gene_stable_id", "homology_protein_stable_id",
    "homology_species", "homology_identity", "dn", "ds", "goc_score",
    "wga_coverage", "is_high_confidence", "homology_id",
]

ONE2MANY_TYPES = {"ortholog_one2many", "ortholog_many2many"}


def find_dump_url() -> str:
    with urllib.request.urlopen(FTP_DIR, timeout=60) as resp:
        html = resp.read().decode("utf-8", "replace")
    files = re.findall(r'href="([^"]*homologies\.tsv\.gz)"', html)
    if not files:
        files = re.findall(r'href="([^"]+\.tsv\.gz)"', html)
    protein = [f for f in files if "protein_default" in f]
    chosen = (protein or files)[0]
    return FTP_DIR + chosen


def download(url: str) -> Path:
    dest = RAW_DIR / url.rsplit("/", 1)[1]
    if dest.exists() and dest.stat().st_size > 0:
        print(f"skip (exists): {dest}")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    print(f"download: {url}")
    with urllib.request.urlopen(url, timeout=600) as resp, open(tmp, "wb") as fh:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)
    tmp.rename(dest)
    return dest


def build_tables(dump: Path) -> None:
    ORTHOLOGS.mkdir(parents=True, exist_ok=True)
    one2one_rows = {sp: [] for sp in SPECIES.values()}
    paralog_rows = {sp: [] for sp in SPECIES.values()}

    with gzip.open(dump, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pair = {row["species"], row["homology_species"]}
            if "homo_sapiens" not in pair:
                continue
            target = (pair - {"homo_sapiens"})
            if not target:
                continue
            target = target.pop()
            if target not in one2one_rows:
                continue
            ht = row["homology_type"]
            if ht == "ortholog_one2one":
                one2one_rows[target].append(row)
            elif ht in ONE2MANY_TYPES and row["species"] == "homo_sapiens":
                paralog_rows[target].append(row)

    for short, ens in SPECIES.items():
        out = ORTHOLOGS / f"human_vs_{ens}_one2one.tsv"
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                               lineterminator="\n", extrasaction="ignore")
            w.writeheader()
            w.writerows(one2one_rows[ens])
        print(f"wrote {out} rows={len(one2one_rows[ens])}")

        rows = paralog_rows[ens]
        counts = defaultdict(set)
        for r in rows:
            counts[r["gene_stable_id"]].add(r["homology_gene_stable_id"])
        out = ORTHOLOGS / f"human_vs_{ens}_paralog_rich.tsv"
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLUMNS + ["n_target_loci"],
                               delimiter="\t", lineterminator="\n",
                               extrasaction="ignore")
            w.writeheader()
            for r in rows:
                r = dict(r)
                r["n_target_loci"] = len(counts[r["gene_stable_id"]])
                w.writerow(r)
        print(f"wrote {out} rows={len(rows)}")


def main() -> int:
    url = find_dump_url()
    dump = download(url)
    build_tables(dump)
    return 0


if __name__ == "__main__":
    sys.exit(main())
