# RESCUE test dataset (mouse, paralog regime)

A minimal, self-contained dataset for testing RESCUE installation and usage
end-to-end in a few minutes: 40 human query proteins whose mouse orthologs
have 2–6 paralogous loci, plus a 26 Mb mini-genome containing exactly those
loci (±250 kb flanks). All coordinates in the ground truth and annotation
files are remapped to the mini-genome.

Derived deterministically from the full benchmark data (Ensembl release 113);
see `docs/DATA.md` for provenance and `scripts/prepare_data/make_test_dataset.py`
for the exact construction.

## Layout

```
data/
  processed/
    genomes/mouse.fa                        mini genome (40 regions, 26 Mb)
    gene_structures/mouse_CDS.bed           CDS of the 634 genes in the regions
    proteins/mouse_canonical.fa             canonical proteins of those genes
    orthologs/human_vs_mus_musculus_one2one.tsv   ortholog rows for those genes
  v3/mouse_paralog/
    proteins.fa                             40 human query proteins
    ground_truth.json                       primary + alternative GT loci (remapped)
    paralog_info.json                       per-query target locus metadata
    pairs.tsv                               query ↔ target gene pairs
    meta.json
```

## Usage

See the "Quick Start with the Test Dataset" section of the repository
README. Expected runtime: ~2–5 minutes on 8 cores. Expected outputs: per-tool
locus recall (`results/tprs/recall/phase1_recall.tsv`), C-stage manifest
(`results/tprs/cstage/manifest_mouse.tsv`), and the co-ortholog recovery
summary (phase 3 log). Router training and final evaluation (phases 4–5) are
part of the full reproduction path; see the repository README.

Note: metrics computed on 40 proteins are a smoke test of the installation,
not a reproduction of the paper numbers.
