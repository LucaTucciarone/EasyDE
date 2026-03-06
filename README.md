# EasyDE

A modular R pipeline for pseudobulk differential expression analysis with
RUVseq correction and gene set enrichment. Designed for scRNA-seq data
([PanKbase](workflow/scripts/fetch/PANKBASE.md) and similar resources)
but works with any pseudobulk count matrices.

Orchestrated by Snakemake for local or SLURM cluster execution.

---

## Table of Contents

1.  [Overview](#overview)
2.  [Installation](#installation)
3.  [Quick Start](#quick-start)
4.  [Project Layout](#project-layout)
5.  [Configuration](#configuration)
6.  [Running with Snakemake](#running-with-snakemake)
7.  [Running on SLURM](#running-on-slurm)
8.  [Running Step by Step](#running-step-by-step)
9.  [Methods](#methods)
10. [Output Files](#output-files)
11. [Exploratory Notebook](#exploratory-notebook)
12. [Resource Files](#resource-files)
13. [Troubleshooting](#troubleshooting)

---

## Overview

The pipeline runs in seven sequential steps. Steps 02-06 are parallelized
across every **contrast x stratum** (cell type) combination, while steps 01
and 07 run once per contrast:

```
01  Validate config + contrasts          (once per run — soft gate)
02  Filter samples, build coldata        (per contrast x stratum)
03  Base DESeq2                          (per contrast x stratum)
04  RUVseq normalization + k selection   (per contrast x stratum)
05  Final DESeq2 with W factors          (per contrast x stratum)
06  fGSEA pathway enrichment             (per contrast x stratum)
07  Aggregate results                    (once per contrast)
```

Step 01 validates all contrasts and flags any that have issues (bad column
references, missing groups, etc.) without blocking the rest. Steps 02-06
check preflight status and **skip** flagged strata gracefully, writing
sentinel files so Snakemake always sees the expected outputs.

**All scripts must be run from the pipeline root directory.** Results can be
written anywhere (set `outputs.results_dir` in your config), but code
expects the `workflow/`, `config/`, and `contrasts/` paths relative to `pwd`.

---

## Installation

### 1. Install micromamba

```bash
# macOS (Homebrew)
brew install micromamba

# Linux
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
```

### 2. Create the environment

```bash
cd /path/to/EasyDE
micromamba env create -f installation/EasyDE_install.yml
```

This installs R, DESeq2, RUVSeq, fgsea, Snakemake, and all dependencies.
Takes ~5-10 minutes on first run.

### 3. Activate

```bash
micromamba activate EasyDE
```

Activate this environment every time before running any pipeline command.

### 4. Verify

```bash
Rscript -e "cat('R ok\n')"               # should print "R ok" and exit 0
snakemake --version                        # should print >=7.x
```

> **Using conda/mamba?** Replace `micromamba` with `conda` or `mamba` — the
> syntax is identical.

<details>
<summary><strong>Known issue: Conda R exit code 255</strong></summary>

Some conda-installed R builds return exit code 255 on completion regardless
of success. Test with `Rscript -e "cat('hello\n')"; echo $?` — it should
return 0. If not, either reinstall R (`mamba install -c conda-forge r-base
--force-reinstall`) or point the Snakefile to your system R.

</details>

---

## Quick Start

```bash
# 1. Activate environment
micromamba activate EasyDE

# 2. Fetch data from PanKbase (optional — skip if you have your own data)
Rscript workflow/scripts/fetch/pankbase_helpers.R \
    --config config/config_traits.yaml \
    --pankbase-config workflow/scripts/fetch/pankbase_config.yaml \
    --strata "Beta|Alpha|Delta"

# 3. Run the full pipeline (all contrasts, all strata)
snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml

# 4. Check results
cat results/Diabetic_vs_ND/summary/contrast_summary.csv | column -t -s,
```

---

## Project Layout

```
EasyDE/
├── Snakefile                             <- Snakemake workflow (do not edit)
│
├── config/                               <- EDIT THESE
│   ├── config_traits.yaml                <- unpaired analysis config
│   └── config_treatments.yaml            <- paired analysis config
│
├── contrasts/                            <- EDIT THESE
│   ├── contrasts.csv                     <- your contrast definitions
│   ├── contrasts_treatments.csv          <- treatment contrast definitions
│   └── contrasts_template.csv            <- documented template to copy from
│
├── profiles/                             <- Snakemake execution profiles
│   ├── local/config.yaml                 <- local (laptop/workstation)
│   └── slurm/config.yaml                 <- SLURM cluster
│
├── data/                                 <- YOUR DATA (or fetched from PanKbase)
│   ├── counts/                           <- one CSV per cell type
│   │   ├── Beta.counts.csv
│   │   ├── Alpha.counts.csv
│   │   └── ...
│   ├── sample_metadata.csv               <- one row per biosample
│   └── cache/                            <- PanKbase download cache
│
├── resources/                            <- shipped with the repo
│   └── gsea_files/
│       ├── reactome_kegg.gmt.txt         <- merged pathway set (default)
│       ├── rpl_genes.csv                 <- ribosomal large subunit exclusion
│       ├── rps_genes.csv                 <- ribosomal small subunit exclusion
│       └── mtr_genes.csv                 <- mitochondrial gene exclusion
│
├── workflow/scripts/                     <- pipeline code (do not edit)
│   ├── 01_prepare_inputs.R               <- validation (soft gate)
│   ├── 02_prepare_coldata.R              <- sample filtering + coldata
│   ├── 03_run_deseq_base.R              <- base DESeq2
│   ├── 04_run_ruvseq.R                  <- RUVseq k selection
│   ├── 05_run_deseq_final.R             <- final DESeq2 + W factors
│   ├── 06_run_fgsea.R                   <- gene set enrichment
│   ├── 07_aggregate_results.R           <- per-contrast summary
│   ├── fetch/                            <- PanKbase data adapter
│   │   ├── pankbase_helpers.R
│   │   ├── pankbase_config.yaml
│   │   └── PANKBASE.md                   <- PanKbase documentation
│   └── utils/                            <- shared R utilities
│       ├── io_utils.R
│       ├── filter_utils.R
│       ├── logging_utils.R
│       ├── plot_utils.R
│       ├── stats_utils.R
│       └── validation_utils.R
│
├── notebooks/                            <- post-run exploration
│   └── explore_results.ipynb             <- R/Jupyter results dashboard
│
├── installation/
│   └── EasyDE_install.yml                <- conda/micromamba environment
│
├── results/                              <- created on first run (see Output Files)
└── logs/                                 <- created on first run
```

### What to edit

| What | Where | When |
|------|-------|------|
| Analysis parameters | `config/config_traits.yaml` | Always — set paths, thresholds, toggles |
| Contrast definitions | `contrasts/contrasts.csv` | Always — define your comparisons |
| Execution settings | `profiles/local/config.yaml` | Optional — adjust core count, latency |
| Data files | `data/counts/*.counts.csv` + `data/sample_metadata.csv` | Provide your own or fetch from PanKbase |

### What NOT to edit

The `workflow/scripts/` directory, `Snakefile`, and `resources/` are part of
the pipeline distribution. Edits there will break reproducibility.

---

## Configuration

### Config files

EasyDE separates **what** to analyze from **how** to execute:

| File | Purpose |
|------|---------|
| `config/config_traits.yaml` | Unpaired analysis: thresholds, paths, step toggles, DESeq2/RUVseq/fGSEA parameters |
| `config/config_treatments.yaml` | Paired analysis: same structure, with `paired_analysis: true` and treatment-specific contrasts |
| `profiles/local/config.yaml` | Local execution: core count, keep-going, latency wait |
| `profiles/slurm/config.yaml` | SLURM execution: job count, sbatch template, account/partition |

The pipeline config (traits or treatments) controls the R scripts. The
profile controls Snakemake's execution engine. They are independent.

### Contrast definitions

Copy the template and edit:

```bash
cp contrasts/contrasts_template.csv contrasts/contrasts.csv
```

Each row defines one contrast. Key columns:

| Column | Description | Example |
|--------|-------------|---------|
| `contrast_id` | Unique name | `Diabetic_vs_ND` |
| `contrast_var` | Metadata column to test | `Derived diabetes status` |
| `biosample_id_col` | Sample ID column | `sample_accession` |
| `control_grp` | Reference level (leave blank for continuous) | `Normal` |
| `experimental_grp` | Test level (leave blank for continuous) | `Diabetes` |
| `strata` | Cell types, pipe-separated | `Beta\|Alpha\|Delta` |
| `covariates` | Confounders, pipe-separated | `Age (years)\|Gender\|BMI` |
| `latent_corr_vars` | Variables to check for W-factor correlation | `Age (years)\|Gender\|BMI` |
| `filter_col_N` / `filter_val_N` | Subset filters (up to 3 pairs) | `treatment` / `no_treatment` |

See `contrasts/contrasts_template.csv` for full column documentation with
examples for binary, continuous, and filtered contrasts.

### Count matrix format

One CSV per stratum at `data/counts/{stratum_safe_name}.counts.csv`. The safe
name replaces spaces and special characters with underscores
(`"CD4 T cell"` becomes `CD4_T_cell.counts.csv`).

```
gene,SAMN001,SAMN002,SAMN003,...
INS,142,0,88,...
GCG,0,201,0,...
```

Genes as rows, samples as columns. First column is the gene symbol.

### Sample metadata format

`data/sample_metadata.csv` — one row per biosample, must contain:
- The `biosample_id_col` from your contrasts file (e.g. `sample_accession`)
- The contrast variable column (e.g. `Derived diabetes status`)
- All covariates listed in `covariates`
- The donor ID column set in `filtering.donor_col` (e.g. `donor_accession`)

---

## Running with Snakemake

Snakemake is the recommended way to run EasyDE. It handles the full DAG of
dependencies, parallelizes steps 02-06 across all contrast x stratum
combinations, and automatically retries failed jobs.

### Run all contrasts

```bash
micromamba activate EasyDE

snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml
```

This runs **all** contrasts defined in your contrasts CSV across **all**
their strata. For 20 contrasts x 11 cell types, that is ~1,100 parallel jobs.

### Run one contrast

```bash
snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml \
    results/Diabetic_vs_ND/summary/contrast_summary.csv
```

### Dry run (see what would execute)

```bash
snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml -n
```

### Run treatment (paired) analysis

```bash
snakemake --profile profiles/local \
    --config pipeline_config=config/config_treatments.yaml
```

### Clean and re-run

```bash
# Remove results for one contrast
rm -rf results/Diabetic_vs_ND/

# Remove all results
rm -rf results/ logs/

# Then re-run
snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml
```

---

## Running on SLURM

For cluster execution, use the SLURM profile. Edit
`profiles/slurm/config.yaml` to set your account, partition, and QOS:

```yaml
# profiles/slurm/config.yaml (edit these)
default-resources:
  slurm_account: your_account
  slurm_partition: your_partition
  slurm_qos: your_qos
  mem_mb: 8000
  time: 240      # minutes
  cpus: 1
```

Then run:

```bash
snakemake --profile $PWD/profiles/slurm \
    --config pipeline_config=config/config_traits.yaml
```

> **Note:** Use `$PWD/profiles/slurm` (absolute path) for the SLURM profile.
> SLURM jobs inherit the working directory, so all paths resolve correctly.

### SLURM with a cluster config

For production runs with absolute output paths (e.g. writing results to a
shared project directory), create a cluster-specific config:

```bash
snakemake --profile $PWD/profiles/slurm \
    --config pipeline_config=config/TSCC_config_Run_452026_allTraits.yaml
```

---

## Running Step by Step

For debugging or educational purposes, you can run each script manually.
All commands assume you are in the pipeline root directory.

```bash
micromamba activate EasyDE

CONFIG="config/config_traits.yaml"
CONTRAST="Diabetic_vs_ND"
STRATUM="Beta"

# Step 01 — Validate (once per run)
Rscript workflow/scripts/01_prepare_inputs.R --config "$CONFIG"

# Step 02 — Prepare coldata
Rscript workflow/scripts/02_prepare_coldata.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 03 — Base DESeq2
Rscript workflow/scripts/03_run_deseq_base.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 04 — RUVseq
Rscript workflow/scripts/04_run_ruvseq.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 05 — Final DESeq2 with W factors
Rscript workflow/scripts/05_run_deseq_final.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 06 — fGSEA
Rscript workflow/scripts/06_run_fgsea.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

# Step 07 — Aggregate (once per contrast, no --stratum needed)
Rscript workflow/scripts/07_aggregate_results.R \
    --config "$CONFIG" --contrast "$CONTRAST"
```

---

## Methods

### DESeq2 (steps 03 and 05)

Standard DESeq2 pseudobulk analysis. Step 03 fits the base model with
user-specified covariates. If RUVseq is enabled, step 05 re-fits the model
with additional W latent factors to control for unwanted technical variation.
Both steps produce log2 fold changes (MLE and apeglm-shrunk) and BH-adjusted
p-values.

For **continuous** contrasts (e.g. Age, BMI, HbA1c), the pipeline detects
that `control_grp` and `experimental_grp` are blank and treats the contrast
variable as numeric. Log2 fold changes represent the effect per unit increase.

For **paired** analysis (`paired_analysis: true` in config), the pipeline
matches samples by donor across treatment conditions and includes donor as a
blocking factor in the DESeq2 formula.

### RUVseq (step 04)

RUVseq removes unwanted technical variation using empirical negative control
genes — genes that are not differentially expressed (base DESeq2 p-value
above `ruvseq.empirical_pval_threshold`, default 0.5). These genes capture
technical batch effects without biological signal.

**Automatic k selection:** The pipeline runs `RUVg` for k = 1 through
`ruvseq.max_latent_factors` (default 10). For each k, it computes the
Relative Log Expression (RLE) matrix and measures between-sample variance via
a one-way ANOVA F-statistic. The optimal k is chosen at the "elbow" of the
F-statistic curve — the point where the second derivative is maximized,
indicating the transition from meaningful variance reduction to diminishing
returns. This is analogous to the scree-plot elbow method used in PCA.

After selecting the best k, the pipeline tests each W factor for correlation
with the contrast variable (linear regression, p < 0.05). W factors that
correlate with the variable of interest are pruned from the final formula to
avoid regressing out real biological signal.

### fGSEA (step 06)

Fast Gene Set Enrichment Analysis using the Wald test statistic as ranking
metric. Ribosomal (RPL/RPS) and mitochondrial (MT-) genes are excluded
before ranking to prevent these highly variable housekeeping genes from
dominating pathway results.

Runs on base results, RUV-corrected results, or both (set `fgsea.run_on`
in config). Default gene set: merged Reactome + KEGG pathways from MSigDB.

### Preflight validation (steps 01-02)

Step 01 performs structural validation: checks that all files exist, required
columns are present, contrasts are well-formed, and groups have matching
samples. Per-contrast issues are **flagged but not fatal** — only structural
problems (missing CSV columns, duplicate contrast IDs) halt the pipeline.

Step 02 performs data-level validation for each stratum: checks count matrix
values (no negatives, no all-zero rows), verifies design matrix rank, and
detects near-singular covariates. Results are written to
`intermediates/preflight.csv` and propagated downstream — steps 03-06 read
the preflight status and skip failed strata gracefully.

---

## Output Files

Results can be written to any directory by setting `outputs.results_dir` in
your config file. The structure within that directory is:

```
results/
└── {contrast_id}/
    ├── {stratum}/
    │   ├── intermediates/                       working files
    │   │   ├── coldata.csv                      step 02 — filtered sample metadata
    │   │   ├── counts_filtered.csv              step 02 — filtered count matrix
    │   │   ├── name_map.csv                     step 02 — original<->safe name map
    │   │   ├── preflight.csv                    step 02 — validation results
    │   │   ├── dds_base.rds                     step 03 — DESeq2 object
    │   │   ├── model_info_base.csv              step 03 — model formula + stats
    │   │   ├── ruvseq_best_coldata.csv          step 04 — coldata + W factors
    │   │   ├── ruvseq_summary.csv               step 04 — k selection + pruning
    │   │   ├── ruvseq_anova.csv                 step 04 — F-stat per k
    │   │   ├── dds_ruv.rds                      step 05 — final DESeq2 object
    │   │   └── model_info_ruv.csv               step 05 — RUV model formula + stats
    │   │
    │   ├── finals/                              clean result tables
    │   │   ├── results_base.tsv                 step 03 — base DE results (all genes)
    │   │   ├── results_ruv.tsv                  step 05 — RUV DE results (all genes)
    │   │   ├── fgsea_base.tsv                   step 06 — base pathway enrichment
    │   │   └── fgsea_ruv.tsv                    step 06 — RUV pathway enrichment
    │   │
    │   └── plots/
    │       ├── volcano_base.pdf                 step 03
    │       ├── volcano_ruv.pdf                  step 05
    │       ├── ruvseq_diagnostics.pdf           step 04 — correlation + elbow plot
    │       ├── fgsea_base.pdf                   step 06
    │       └── fgsea_ruv.pdf                    step 06
    │
    └── summary/                                 per-contrast aggregation (step 07)
        ├── contrast_summary.csv                 one row per stratum, all counts
        ├── preflight_summary.csv                preflight results across strata
        ├── {contrast}_results_base.tsv          merged base DE (all strata)
        ├── {contrast}_results_ruv.tsv           merged RUV DE (all strata)
        ├── {contrast}_fgsea_base.tsv            merged base fGSEA (all strata)
        ├── {contrast}_fgsea_ruv.tsv             merged RUV fGSEA (all strata)
        ├── errors.log                           errors only (if any)
        ├── warnings.log                         warnings only (if any)
        └── log_summary.tsv                      all flagged log lines
```

### DE result columns (`results_base.tsv` / `results_ruv.tsv`)

| Column | Description |
|--------|-------------|
| `gene` | Gene symbol |
| `baseMean` | Mean normalized count across all samples |
| `log2FoldChange` | DESeq2 MLE estimate |
| `log2FoldChange_shrunk` | apeglm-shrunk estimate (NA if shrinkage disabled/failed) |
| `lfcSE` | Standard error of log2FoldChange |
| `stat` | Wald statistic (default fGSEA ranking metric) |
| `pvalue` | Raw p-value |
| `padj` | BH-adjusted p-value |

### Pipeline status values (`contrast_summary.csv`)

| Status | Meaning |
|--------|---------|
| `success` | All steps completed, DE results produced |
| `partial` | Coldata built but DESeq2 did not produce results (e.g. too few samples after filtering, model fitting failed) |
| `failed` | Stratum could not be processed at all |
| `preflight_fail` | Validation detected a blocking issue in step 02 |

---

## Exploratory Notebook

After running the pipeline, use the Jupyter notebook to visualize results
across all contrasts and cell types:

```bash
cd notebooks/
jupyter notebook explore_results.ipynb
```

The notebook produces:
- **DE features boxplot** — RUV-corrected DEG counts per contrast, colored by cell type
- **FGSEA heatmaps** — top significant pathways with significance stars
- **Diagonal pathway heatmap** — pathways ordered by their dominant contrast x cell type
- **Cell-type-divergent pathways** — pathways enriched in opposite directions across cell types
- **Pipeline status tile plot** — contrast x cell type completion overview

---

## Resource Files

All resources are included in the repository — no additional downloads needed.

### Gene exclusion lists

Used by step 06 to remove ribosomal and mitochondrial genes before fGSEA ranking:

```
resources/gsea_files/rpl_genes.csv     ribosomal proteins, large subunit
resources/gsea_files/rps_genes.csv     ribosomal proteins, small subunit
resources/gsea_files/mtr_genes.csv     mitochondrially encoded genes
```

### Pathway GMT files

Three MSigDB gene set collections, pre-merged into a single default file:

```
resources/gsea_files/c2.cp.reactome.v*.symbols.gmt.txt    Reactome
resources/gsea_files/c2.cp.kegg_medicus.v*.symbols.gmt.txt KEGG Medicus
resources/gsea_files/c2.cp.kegg_legacy.v*.symbols.gmt.txt  KEGG Legacy
resources/gsea_files/reactome_kegg.gmt.txt                 merged (default)
```

To use a different collection, set `fgsea.gene_sets_file` in your config.

**Updating pathways:** Download new GMT files from
[MSigDB](https://www.gsea-msigdb.org/gsea/downloads.jsp) and re-merge:

```bash
cat c2.cp.reactome.*.gmt c2.cp.kegg_medicus.*.gmt c2.cp.kegg_legacy.*.gmt \
    > resources/gsea_files/reactome_kegg.gmt.txt
```

---

## Troubleshooting

**`every gene contains at least one zero`**
DESeq2 cannot estimate size factors. Usually caused by too few samples or
extreme sparsity. Reduce `filtering.min_gene_counts` or check that you have
enough samples in both groups.

**`lfcShrink failed`**
apeglm requires at least one sample per group with non-zero counts. This is
a warning — `log2FoldChange_shrunk` will be NA but `log2FoldChange` is valid.

**`no count matrix found for stratum 'X'`**
The count file at `data/counts/X.counts.csv` (after safe-name conversion)
does not exist. Check that `strata` values in contrasts.csv match filenames
in `data/counts/`. Spaces and special characters become underscores.

**`W pruning removed all latent factors`**
All W factors correlated with the contrast variable (p < 0.05). The final
model runs without W factors — equivalent to base DESeq2. Consider whether
sample size is sufficient for RUVseq.

**fGSEA returns 0 significant pathways**
Not necessarily an error. Check: (1) the ranking statistic produces a
reasonable spread, (2) gene symbols in the GMT match your gene names,
(3) `fgsea.min_gene_set_size` is not too high.

**Snakemake shows "Nothing to be done"**
All output files already exist. Delete `results/` (or the specific contrast
directory) to force a re-run.

### Reading log files

Steps 02-06 write structured logs to `logs/{contrast_id}/{stratum}.log`:

```
[2026-03-06 14:32:01] [INFO]  contrast=Diabetic_vs_ND | biosample=Beta | step=start | ...
[2026-03-06 14:32:03] [INFO]  contrast=Diabetic_vs_ND | biosample=Beta | step=filter_samples | n_samples=89
[2026-03-06 14:32:15] [WARN]  contrast=Diabetic_vs_ND | biosample=Beta | step=lfc_shrinkage | msg=...
```

To see all warnings and errors after a run:

```bash
grep -rh "WARN\|ERROR\|SKIP" logs/
```

Step 07 also collects all flagged lines into
`results/{contrast}/summary/log_summary.tsv` for easy review.
