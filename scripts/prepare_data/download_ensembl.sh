#!/usr/bin/env bash
# Download all primary data for the RESCUE benchmark from the Ensembl
# release 113 FTP archive. Files land in $RESCUE_DATA_ROOT/raw/release-113/
# (default: ./data/raw/release-113/).
#
# Usage: bash scripts/prepare_data/download_ensembl.sh
set -euo pipefail

DATA_ROOT="${RESCUE_DATA_ROOT:-$PWD/data}"
RAW="$DATA_ROOT/raw/release-113"
mkdir -p "$RAW"

FTP="https://ftp.ensembl.org/pub/release-113"

# short_name  ensembl_name  FileName  assembly
SPECIES="
mouse       mus_musculus           Mus_musculus           GRCm39
dog         canis_lupus_familiaris Canis_lupus_familiaris ROS_Cfam_1.0
chicken     gallus_gallus          Gallus_gallus          bGalGal1.mat.broiler.GRCg7b
zebrafish   danio_rerio            Danio_rerio            GRCz11
human       homo_sapiens           Homo_sapiens           GRCh38
"

echo "$SPECIES" | while read -r short sp spfile asm; do
    [ -z "$short" ] && continue
    for kind in dna pep gtf; do
        case "$kind" in
            dna) url="$FTP/fasta/$sp/dna/$spfile.$asm.dna.toplevel.fa.gz" ;;
            pep) url="$FTP/fasta/$sp/pep/$spfile.$asm.pep.all.fa.gz" ;;
            gtf) url="$FTP/gtf/$sp/$spfile.$asm.113.gtf.gz" ;;
        esac
        dest="$RAW/$(basename "$url")"
        if [ -s "$dest" ]; then
            echo "skip (exists): $dest"
        else
            echo "download: $url"
            curl -fSL --retry 3 -o "$dest.tmp" "$url"
            mv "$dest.tmp" "$dest"
        fi
    done
done

echo "Done. Raw files in $RAW"
