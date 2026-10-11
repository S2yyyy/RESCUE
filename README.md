# RESCUE: Protein-to-Genome Alignment Experiments

This repository contains the release version of the code and minimal summary
artifacts for the RESCUE paper experiments. Some directory and script names
use the development name `TPRS` (e.g. `scripts/tprs/`); these are internal
identifiers only — the method is called RESCUE throughout.

The repository intentionally does not include large reference genomes,
proteomes, external aligner installations, or machine-specific run logs. Full
reruns require the external data layout described below. The checked-in
summary artifacts under `results/` document the quantitative results
reported in the paper. The paralog-rich benchmark sets evaluated in the
paper (1,988 proteins across four species) are included under
`data/benchmark/`. A small test dataset (40 proteins, 26 Mb mini-genome)
is included under `data/` for the quick start below.

## Repository Layout

- `scripts/stage0/`: normalization of heterogeneous tool outputs into a common
  candidate-region BED format, plus locus-recall utilities.
- `scripts/exp1/`: end-to-end benchmark runners and metrics for existing
  protein-to-genome alignment tools.
- `scripts/exp5/`: forced-locus diagnosis used to separate identity ceiling
  from locus-search failure.
- `scripts/exp6/`: fixed hybrid and oracle analyses for the candidate-selection
  gap.
- `scripts/tprs/`: RESCUE pipeline Phases 1-5 (Stage I candidate generation,
  Stage II learned joint selection), ablations, heuristic baselines, and
  feature-importance extraction.
- `scripts/prepare_data/`: deterministic rebuild of the full benchmark data
  tree from Ensembl release 113 (see `docs/DATA.md`).
- `results/`: minimal checked-in summary artifacts.
- `tests/`: release integrity checks.
- `docs/DATA.md`: data sources, versions, benchmark construction, and
  checksums for all datasets used in the paper.
- `data/benchmark/`: the paralog-rich benchmark sets evaluated in the paper
  (see "Benchmark Data" below).

## Benchmark Data

The paralog-rich benchmark is provided under
`data/benchmark/{mouse,dog,chicken,zebrafish}_paralog/` (488–500 proteins per
species, 1,988 cases in total). Each case is a human query protein whose gene
has at least two annotated orthologous loci in the target species (Ensembl
Compara release-113 `ortholog_one2many` / `ortholog_many2many`). Each
directory contains the query proteins (`proteins.fa`), the reference loci
(`ground_truth.json`), the ortholog-pair metadata (`pairs.tsv`), and the full
paralog annotation (`paralog_info.json`). Construction details and data
sources are documented in `docs/DATA.md`, and the sets can also be rebuilt
deterministically from Ensembl release 113 with `scripts/prepare_data/`. The
40-protein dataset under `data/` remains the quick-start smoke test described
below and is unrelated to this benchmark data.

## Experiment Map

| Paper claim or output | Code | Included summary artifact |
| --- | --- | --- |
| End-to-end benchmark of existing aligners | `scripts/exp1/` | `results/exp1/summary.tsv` |
| Forced-locus diagnosis of identity ceiling and search gap | `scripts/exp5/` | `results/exp5/cell/per_cell.tsv` |
| Best fixed recipe and oracle candidate-selection gap | `scripts/exp6/` | — |
| RESCUE candidate generation, alignment scoring, recovery, and learned router | `scripts/tprs/phase1_recall.py`, `phase2_cstage.py`, `phase3_recovery.py`, `phase4_router.py`, `phase5_run.py` | `results/tprs/phase5/evaluation_summary.json`, `results/tprs/phase5/ablation_summary_all_species.json` |
| Component ablation E2-E6 | `scripts/tprs/phase5_ablation.py` | `results/tprs/phase5/ablation_summary_all_species.json` |
| Non-learned router heuristic baselines | `scripts/tprs/phase5_heuristic_baselines.py` | `results/tprs/phase5/heuristic_baselines.json` |
| Router feature importance | `scripts/tprs/extract_phase5_feature_importance.py` | `results/tprs/phase5/feature_importance_phase5.json` |

## Installation

RESCUE requires Python 3.10 or later (developed and tested with Python 3.11).

Install the Python package and its four core dependencies
(numpy, pandas, scikit-learn, xgboost):

```bash
python -m pip install .
```

This provides the `rescue` command-line entry point plus the per-phase
commands `rescue-recall`, `rescue-cstage`, `rescue-recovery`, and
`rescue-router`. Alternatively, run without installing by using
`python -m pip install -r requirements.txt` and invoking the scripts under
`scripts/tprs/` directly.

RESCUE calls five external bioinformatics binaries at runtime: DIAMOND,
MMseqs2, LAST, miniprot, and samtools. All five are available from Bioconda:

```bash
conda install -c bioconda diamond mmseqs2 last miniprot samtools
```

The paper experiments used DIAMOND 2.1.24, MMseqs2 18.8cc5c, LAST 1651,
miniprot 0.18, and samtools 1.5. Binary paths can be overridden with the
environment variables `DIAMOND_BIN`, `MMSEQS_BIN`, `LASTAL_BIN`,
`MINIPROT_BIN`, and `SAMTOOLS_BIN`.

## Quick Start with the Test Dataset

A small self-contained test dataset (40 mouse paralog proteins, 26 Mb
mini-genome) ships with the repository under `data/` (see `data/README.md`).
It runs the RESCUE pipeline end-to-end in a few minutes on 8 cores.

```bash
# 1. Index the mini-genome and build the search databases
samtools faidx data/processed/genomes/mouse.fa
miniprot -d data/processed/genomes/mouse.mpi data/processed/genomes/mouse.fa
mkdir -p results/exp3/db/diamond_protdb results/exp3/db/last_protdb \
         results/exp3/db/mmseqs2/mouse
diamond makedb --in data/v3/mouse_paralog/proteins.fa \
    -d results/exp3/db/diamond_protdb/mouse.paralog --quiet
lastdb -p -c -R01 results/exp3/db/last_protdb/mouse.paralog \
    data/v3/mouse_paralog/proteins.fa
mmseqs createdb data/processed/genomes/mouse.fa \
    results/exp3/db/mmseqs2/mouse/genomedb --dbtype 2

# 2. Build the ground-truth tables
python scripts/prepare_data/build_gt_tables.py \
    --out-dir results/exp2/gt --species mouse

# 3. Stage I: candidate generation (relaxed search, 4 tools)
bash scripts/tprs/run_relaxed_search.sh all
bash scripts/tprs/run_relaxed_search.sh miniprot-relaxed

# 4. Normalize outputs and evaluate locus recall
mkdir -p results/tprs/stage0_bed
python -m scripts.stage0.normalize --tool miniprot \
    --input results/tprs/raw/miniprot_relaxed.mouse.paralog.gff \
    --query-lengths data/v3/mouse_paralog/proteins.fa --top-k 50 \
    -o results/tprs/stage0_bed/miniprot_relaxed.mouse.paralog.bed
python scripts/tprs/phase1_recall.py \
    --tools miniprot_ultra diamond_relaxed mmseqs_relaxed last_relaxed

# 5. Stage II preparation: candidate alignment + co-ortholog recovery
python scripts/tprs/phase2_cstage.py --species mouse
python scripts/tprs/phase3_recovery.py --species mouse
```

**Inputs** (all provided by the test dataset): human query proteins
(`data/v3/mouse_paralog/proteins.fa`, FASTA), the mouse mini-genome
(`data/processed/genomes/mouse.fa`, FASTA), and ground-truth loci
(`ground_truth.json`, 0-based half-open CDS exon coordinates).

**Outputs**: per-tool locus recall (`results/tprs/recall/phase1_recall.tsv`),
per-region alignment results (`results/tprs/cstage/manifest_mouse.tsv`), and
the co-ortholog recovery summary (phase 3 log, e.g. "16/17 recovered").

Router training and final evaluation (phases 4–5) involve XGBoost model
training and are part of the full reproduction path rather than this smoke
test; on the test dataset they complete in a few minutes:

```bash
python scripts/tprs/phase4_router.py --species mouse
python scripts/tprs/phase5_run.py mouse
```

Note: with 40 proteins these metrics are a smoke test of the installation,
not a reproduction of the paper results. The full benchmark data can be
regenerated with `scripts/prepare_data/` (see `docs/DATA.md`).

## Minimal Reproduction

Create a Python environment and install the Python dependencies:

```bash
python -m pip install -r requirements.txt
```

Some benchmark-harness scripts additionally require `pyyaml` and `pysam`;
the figure generators under `scripts/exp5/figures/` require `matplotlib`.
None of these are needed to run RESCUE itself.

Run the release smoke checks:

```bash
python -m scripts.stage0.test_stage0
python tests/test_release_integrity.py
python -m compileall -q scripts
```

These commands use the checked-in summary artifacts and do not require raw
genomes or aligner binaries.

## Full Rerun Inputs

Full reruns require external reference data and command-line aligners. The
full data tree can be rebuilt deterministically from Ensembl release 113 with
the pipeline in `scripts/prepare_data/` (see `scripts/prepare_data/README.md`
and `docs/DATA.md`). The code uses repository-relative defaults and can be
pointed to external resources with:

```bash
export RESCUE_REPO="$PWD"
export RESCUE_DATA_ROOT="$PWD/data"
export RESCUE_GENOME_DIR="$RESCUE_DATA_ROOT/processed/genomes"
export RESCUE_GENE_STRUCTURE_DIR="$RESCUE_DATA_ROOT/processed/gene_structures"
export RESCUE_PROTEOME_DIR="$RESCUE_DATA_ROOT/processed/proteins"
export RESCUE_ORTHOLOG_DIR="$RESCUE_DATA_ROOT/processed/orthologs"
export RESCUE_V3_DATA="$RESCUE_DATA_ROOT/v3"
```

Expected external-data layout:

```text
data/
  v3/
    mouse/
    mouse_paralog/
    dog/
    dog_paralog/
    chicken/
    chicken_paralog/
    zebrafish/
    zebrafish_paralog/
  processed/
    genomes/
    gene_structures/
    proteins/
    orthologs/
```

## Benchmark Tools (optional)

The following tools were used to benchmark existing aligners against RESCUE
in the paper. They are **not** required to install or run RESCUE itself; they
are only needed to reproduce the benchmark comparisons. Versions below are
the ones used in the paper; most are available from Bioconda.

| Tool | Version | Source |
| --- | --- | --- |
| SPALN | 3.0.8 | https://github.com/ogotoh/spaln |
| GenomeThreader | 1.7.1 | http://genomethreader.org/ |
| Scipio (+ Perl) | 1.4.1 | https://www.webscipio.org/ |
| genBlastG / genBlastA | 1.38 / 1.0.4 | http://genome.sfu.ca/genblast/ |
| Exonerate | 2.4.0 | https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate |
| GeneWise (WISE2) | 2.4.1 | https://www.ebi.ac.uk/~birney/wise2/ |
| NCBI BLAST+ (tblastn) | 2.17.0 | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/ |
| BLAT | 35 | https://github.com/ucscGenomeBrowser/kent |
| miniSplice | 0.4 | https://github.com/lh3/minisplice |

Tool paths can be overridden with environment variables such as
`MINIPROT_BIN`, `SAMTOOLS_BIN`, `DIAMOND_BIN`, `MMSEQS_BIN`, `LASTAL_BIN`,
`TBLASTN_BIN`, `GENBLAST_DIR`, `SCIPIO_ROOT`, and `WISECONFIGDIR`.

## Representative Full Commands

Stage 0 normalization and recall:

```bash
bash scripts/tprs/run_relaxed_search.sh all
python scripts/tprs/phase1_recall.py --top-k 50
```

RESCUE learned-selection experiments:

```bash
python scripts/tprs/phase4_router.py --all
python scripts/tprs/phase5_run.py mouse
python scripts/tprs/phase5_ablation.py --all
python scripts/tprs/phase5_heuristic_baselines.py
python scripts/tprs/extract_phase5_feature_importance.py
```

Baseline and forced-locus experiments:

```bash
python scripts/exp1/launch_all.py
python scripts/exp1/summarize.py
python scripts/exp5/run_forced.py --help
python scripts/exp5/aggregate.py
```

## Release Privacy Notes

This release was sanitized to avoid machine-specific paths, server addresses,
local user names, shell history, local environment files, generated caches, and
run logs. The release integrity test includes a privacy-pattern scan over all
checked-in files.
