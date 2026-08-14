# Data preparation pipeline

Rebuilds the complete RESCUE benchmark data tree from Ensembl release 113
primary sources. All steps are deterministic (no random sampling). See
`docs/DATA.md` for the full data documentation.

## Steps

Run from the repository root. All scripts honor the `RESCUE_*` environment
variables documented in the main README (default layout: `data/`).

```bash
# 1. Download Ensembl release 113 genomes, proteomes, and GTFs
#    (~10 GB, lands in data/raw/release-113/)
bash scripts/prepare_data/download_ensembl.sh

# 2. Download Ensembl Compara homology dump and build the per-species
#    ortholog tables (one2one + paralog_rich TSVs)
python scripts/prepare_data/download_orthologs.py

# 3. Unpack genome FASTAs into data/processed/genomes/ (full commands below)

# 4. Build gene-structure BEDs (CDS + exons) from the GTFs
python scripts/prepare_data/build_gene_structures.py

# 5. Build canonical proteomes (longest-CDS transcript per gene)
python scripts/prepare_data/build_canonical_proteome.py

# 6. Build the benchmark query sets (one-to-one and paralog)
python scripts/prepare_data/prepare_one2one.py
python scripts/prepare_data/prepare_paralog.py

# 7. Build the ground-truth tables (one2one / paralog_strict / paralog_broad)
python scripts/prepare_data/build_gt_tables.py

# 8. Verify against the checksums of the data used in the paper
bash scripts/prepare_data/verify_checksums.sh
```

Step 3 in full — unpack each species' toplevel FASTA to its short name:

```bash
gunzip -c data/raw/release-113/Mus_musculus.GRCm39.dna.toplevel.fa.gz \
  > data/processed/genomes/mouse.fa
gunzip -c data/raw/release-113/Canis_lupus_familiaris.ROS_Cfam_1.0.dna.toplevel.fa.gz \
  > data/processed/genomes/dog.fa
gunzip -c data/raw/release-113/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz \
  > data/processed/genomes/chicken.fa
gunzip -c data/raw/release-113/Danio_rerio.GRCz11.dna.toplevel.fa.gz \
  > data/processed/genomes/zebrafish.fa
```

External aligners additionally need indexed genomes (e.g. `samtools faidx`
and `miniprot -d`); those indexes are built by the aligners' own tooling and
are not part of this pipeline.

## Notes

- BED row order is not canonicalized; if `verify_checksums.sh` reports a BED
  mismatch, compare sorted content before concluding a real difference.
- The canonical proteome rule (longest total CDS per gene) and the exact
  query-set selection logic are described in `docs/DATA.md` §2.
