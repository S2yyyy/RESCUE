# RESCUE: Protein-to-Genome Alignment Experiments

This repository contains the release version of the code and minimal summary
artifacts for the RESCUE paper experiments. Some historical scripts use the
development name `TPRS`; in this release, `TPRS` refers to the implementation
of the RESCUE candidate-generation and learned-selection pipeline.

The repository intentionally does not include large reference genomes,
proteomes, external aligner installations, or machine-specific run logs. Full
reruns require the external data layout described below. The included summary
files are sufficient to regenerate the paper tables and plotted result panels.

## Repository Layout

- `scripts/stage0/`: normalization of heterogeneous tool outputs into a common
  candidate-region BED format, plus locus-recall utilities.
- `scripts/exp1/`: end-to-end benchmark runners and metrics for existing
  protein-to-genome alignment tools.
- `scripts/exp5/`: forced-locus diagnosis used to separate identity ceiling
  from locus-search failure.
- `scripts/exp6/`: fixed hybrid and oracle analyses for the candidate-selection
  gap.
- `scripts/tprs/`: RESCUE/TPRS Phases 1-5, ablations, heuristic baselines, and
  feature-importance extraction.
- `paper/v4/`: manuscript source fragments and figure/table generation scripts.
- `results/`: minimal checked-in summary artifacts used by `paper/v4`.
- `tests/`: release integrity checks.

## Experiment Map

| Paper claim or output | Code | Included summary artifact |
| --- | --- | --- |
| Baseline end-to-end benchmark of existing aligners | `scripts/exp1/` | `results/exp1/summary.tsv`, `paper/v2/figures/out/table_exp1_benchmark.csv` |
| Forced-locus diagnosis of identity ceiling and search gap | `scripts/exp5/` | `results/exp5/cell/per_cell.tsv` |
| Fixed hybrid/oracle candidate-selection gap | `scripts/exp6/` | `paper/v4/tables/table_s_exp6_hybrid.tex` |
| RESCUE candidate pool, C-stage scoring, recovery, and router | `scripts/tprs/phase1_recall.py`, `phase2_cstage.py`, `phase3_recovery.py`, `phase4_router.py`, `phase5_run.py` | `results/tprs/phase5/evaluation_summary.json`, `results/tprs/phase5/ablation_summary_all_species.json` |
| Component ablation E2-E6 | `scripts/tprs/phase5_ablation.py` | `results/tprs/phase5/ablation_summary_all_species.json` |
| Non-learned router heuristic baselines | `scripts/tprs/phase5_heuristic_baselines.py` | `results/tprs/phase5/heuristic_baselines.json` |
| Router feature importance | `scripts/tprs/extract_phase5_feature_importance.py` | `results/tprs/phase5/feature_importance_phase5.json` |
| Paper tables and figures | `paper/v4/tables/`, `paper/v4/figures/` | `paper/v4/tables/*.tex` and checked-in summary artifacts |

## Minimal Reproduction

Create a Python environment and install the Python dependencies:

```bash
python -m pip install -r requirements.txt
```

Run the release smoke checks:

```bash
python -m scripts.stage0.test_stage0
python tests/test_release_integrity.py
python -m compileall -q scripts paper/v4/figures paper/v4/tables
```

Regenerate summary-driven paper outputs:

```bash
python paper/v4/tables/gen_new_tables.py
python paper/v4/figures/gen_fig1_paralog_gap.py
python paper/v4/figures/gen_fig2_pipeline.py
python paper/v4/figures/gen_fig_forced_locus.py
python paper/v4/figures/gen_fig4_ablation_forest.py
python paper/v4/figures/gen_fig_feature_importance.py
```

These commands use the checked-in summary artifacts and do not require raw
genomes or aligner binaries.

## Full Rerun Inputs

Full reruns require external reference data and command-line aligners. The code
uses repository-relative defaults and can be pointed to external resources with:

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

External tools used by full reruns include `miniprot`, `samtools`, `diamond`,
`mmseqs`, `lastal`, `spaln`, `gth`, `blat`, `tblastn`, `genblastG`, `genewise`,
and Scipio. Tool paths can be overridden with environment variables such as
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
