# Stage 0 — Intermediate Output Protocol

Standalone post-processing layer. Reads each tool's native output and
produces a unified candidate BED. **No tool binary or source is modified.**

See `EXPERIMENT_PLAN_V4_PER_TOOL_ABC.md` §Stage 0 for the full contract.

## Layout

```
scripts/stage0/
├── common.py        # HSP, Candidate, overlap math, query-length I/O
├── parsers.py       # Per-tool read-only parsers
├── aggregate.py     # 200kb HSP merge + reciprocal-50% dedup + top-K
├── eval_locus.py    # Search-OK / Select-OK evaluator
├── normalize.py     # CLI driver
└── test_stage0.py   # Self-tests
```

## Usage

Normalize a tool's output:

```bash
python -m scripts.stage0.normalize \
    --tool tblastn \
    --input results/tblastn/human_vs_mouse.tsv \
    --query-lengths data/v3/human/proteins.fa \
    --top-k 50 \
    -o results/stage0/tblastn.human_vs_mouse.bed
```

Evaluate Search-OK / Select-OK:

```bash
python -m scripts.stage0.eval_locus \
    --bed results/stage0/tblastn.human_vs_mouse.bed \
    --gt data/v3/mouse/gt_loci.tsv \
    -o results/stage0/tblastn.human_vs_mouse.verdict.tsv
```

## Supported tools

| Tool | Class | Parser |
|---|---|---|
| tblastn, diamond | HSP | BLAST `-outfmt 6` |
| last | HSP | LAST BlastTab |
| blat | HSP | BLAT PSL |
| miniprot (+variants) | native | PAF with `ms:i:` / `AS:i:` tags |
| spaln, gth | native | GFF3 with `Target=` attribute |
| scipio | native | YAML (line-based parser, no PyYAML dep) |
| procompart | native | compartment TSV |
| genblasta | native | `.gblast` pipe-separated |

HSP tools route through `aggregate_hsps` (200 kb same-strand merge, score =
Σ bitscores). Native tools pass through and only undergo dedup + top-K.

## GT table format

TSV with header:

```
query_id  chrom  start  end  strand  gene_id
```

Coordinates are 0-based half-open. Gene body span = `(min(exon_start), max(exon_end))`.

For 1:many broad/strict evaluation, emit one row per acceptable GT locus
(broad = all orthologs from Compara `homology_type ∈ {one2one, one2many,
many2many}`; strict = primary one2one only).

## Tests

```bash
python -m scripts.stage0.test_stage0
```
