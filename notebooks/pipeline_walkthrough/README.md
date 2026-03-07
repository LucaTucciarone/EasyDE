# EasyDE Pipeline Walkthrough — Interactive Notebooks

This folder contains an interactive, notebook-based version of the full
EasyDE pipeline. Instead of Snakemake orchestration, each step is a
standalone R Markdown notebook you can run top-to-bottom in RStudio.

## Notebook Sequence

| Notebook | Pipeline Step | Description |
|----------|---------------|-------------|
| `00_setup.Rmd` | — | Config, libraries, define contrasts and strata loops |
| `01_validate.Rmd` | Step 01 | Validate config and contrasts |
| `02_prepare_coldata.Rmd` | Step 02 | Filter samples, build coldata, preflight QC |
| `03_deseq_base.Rmd` | Step 03 | Base DESeq2 analysis |
| `04_ruvseq.Rmd` | Step 04 | RUVseq normalization (skipped in paired mode) |
| `05_deseq_final.Rmd` | Step 05 | Final DESeq2 with W factors (skipped in paired mode) |
| `06_fgsea.Rmd` | Step 06 | Gene set enrichment analysis |
| `07_aggregate.Rmd` | Step 07 | Aggregate results + status tile plot |

## How to Use

1. Open each notebook in RStudio **in order** (00 → 07)
2. Set your working directory to the pipeline root in `00_setup.Rmd`
3. Run all cells top-to-bottom (`Ctrl+Alt+R` or "Run All")
4. Each notebook loops over all contrasts × strata using simple `for` loops

## Prerequisites

- R environment with DESeq2, RUVSeq, fgsea installed (use the EasyDE conda env)
- Working directory set to the pipeline root (where `Snakefile` lives)
- Data files in `data/counts/` and `data/sample_metadata.csv`
- Config file (e.g. `config/config_traits.yaml`)

## Paired vs Unpaired

In paired mode (`paired_analysis: true` in config), notebooks 04 and 05
will automatically skip — matching the pipeline behavior.
