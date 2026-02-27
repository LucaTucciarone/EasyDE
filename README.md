# EasyDE

A modular R pipeline for pseudobulk differential expression analysis with
RUVseq correction and gene set enrichment. Designed for scRNA-seq data
(PanKbase and similar resources) but works with any pseudobulk count matrices.

---

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Resource Files](#resource-files)
4. [Project Layout](#project-layout)
5. [Configuration](#configuration)
6. [Running Step by Step (Local)](#running-step-by-step-local)
7. [Running with Snakemake (Coming Soon)](#running-with-snakemake)
8. [Running on SLURM (Coming Soon)](#running-on-slurm)
9. [Output Files](#output-files)
10. [Troubleshooting](#troubleshooting)

---

## Overview

The pipeline runs in seven sequential steps, parallelizable at the stratum
(cell type / tissue) level for steps 02–06:

```
01  validate config + contrasts        (once per run)
02  filter samples, build coldata      (per contrast × stratum)
03  base DESeq2                        (per contrast × stratum)
04  RUVseq normalization + k selection (per contrast × stratum)
05  final DESeq2 with W factors        (per contrast × stratum)
06  fGSEA                              (per contrast × stratum)
07  aggregate results                  (once per contrast)
```

Each step reads from and writes to `results/{contrast_id}/{stratum}/intermediates/`.
No step reaches back two levels — intermediate files are the explicit contract
between scripts.

---

## Installation

### R (≥ 4.2 recommended)

```bash
# Check your R version
R --version
```

### CRAN packages

```r
install.packages(c(
  "yaml",
  "data.table",
  "dplyr",
  "tibble",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "ggcorrplot",
  "viridis",
  "reshape2",
  "optparse",
  "stringr"
))
```

### Bioconductor packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2",      # core differential expression
  "RUVSeq",      # remove unwanted variation
  "EDASeq",      # required by RUVSeq
  "fgsea"        # fast gene set enrichment analysis
))
```

### Verify installation

```r
pkgs <- c("yaml","data.table","dplyr","tibble","tidyr","ggplot2","ggrepel",
          "ggcorrplot","viridis","reshape2","optparse","stringr",
          "DESeq2","RUVSeq","EDASeq","fgsea")

not_installed <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(not_installed) == 0) {
    cat("All packages installed.\n")
} else {
    cat("Missing:", paste(not_installed, collapse = ", "), "\n")
}
```

---

## Resource Files

You need two types of resource files: gene exclusion lists (for fGSEA cleanup)
and a pathway GMT file. Create the resources directory structure first:

```bash
mkdir -p resources/gsea_files
```

### Gene exclusion lists

These are used by step 06 to remove ribosomal and mitochondrial genes before
ranking — these genes are technically variable and inflate pathway results.

Each file should be a CSV with at minimum an **`Approved symbol`** column
(HGNC standard). Download from HGNC BioMart:

**Option A — HGNC BioMart (recommended, always current)**

1. Go to: https://biomart.genenames.org
2. Under **Filters**, set:
   - Gene group: select the group you want (see below)
3. Under **Attributes**, ensure **Approved symbol** is included
4. Click **Results** → download as CSV

Gene groups to download separately:
- **Ribosomal proteins — large subunit**: search `RPL` or filter by gene group "Ribosomal proteins"
- **Ribosomal proteins — small subunit**: search `RPS`
- **Mitochondrial genes**: filter by gene group "Mitochondrially encoded"

Save as:
```
resources/gsea_files/rpl_genes.csv
resources/gsea_files/rps_genes.csv
resources/gsea_files/mtr_genes.csv
```

**Option B — Direct download from HGNC FTP**

```bash
# Ribosomal large subunit (RPL)
wget -O resources/gsea_files/rpl_genes.csv \
  "https://www.genenames.org/cgi-bin/genegroup/download?id=1054&type=branch"

# Ribosomal small subunit (RPS)
wget -O resources/gsea_files/rps_genes.csv \
  "https://www.genenames.org/cgi-bin/genegroup/download?id=1053&type=branch"

# Mitochondrially encoded genes
wget -O resources/gsea_files/mtr_genes.csv \
  "https://www.genenames.org/cgi-bin/genegroup/download?id=639&type=branch"
```

> **Note:** Check that the downloaded CSV contains an `Approved symbol` column.
> If the column header differs, either rename it or update `fgsea.exclude_gene_lists`
> logic in `06_run_fgsea.R` accordingly.

**Minimal manual fallback** — if downloads fail, create placeholder CSVs:

```r
# Run this in R from your pipeline root
write.csv(data.frame("Approved symbol" = c("RPL3","RPL4","RPL5"), check.names=FALSE),
          "resources/gsea_files/rpl_genes.csv", row.names=FALSE)
write.csv(data.frame("Approved symbol" = c("RPS3","RPS4","RPS5"), check.names=FALSE),
          "resources/gsea_files/rps_genes.csv", row.names=FALSE)
write.csv(data.frame("Approved symbol" = c("MT-CO1","MT-ND1","MT-ATP6"), check.names=FALSE),
          "resources/gsea_files/mtr_genes.csv", row.names=FALSE)
```

### Pathway GMT file

**MSigDB (recommended)**

1. Go to: https://www.gsea-msigdb.org/gsea/downloads.jsp
2. Create a free account if needed
3. Under **Gene Set Collections**, download:
   - **C2: curated gene sets → REACTOME** (file: `c2.cp.reactome.v*.symbols.gmt`)
   - **C2: curated gene sets → KEGG** (file: `c2.cp.kegg_legacy.v*.symbols.gmt`)
4. Optionally merge them:

```bash
cat c2.cp.reactome.v*.symbols.gmt \
    c2.cp.kegg_legacy.v*.symbols.gmt \
    > resources/gsea_files/reactome_kegg.gmt.txt
```

Or use them separately by updating `fgsea.gene_sets_file` in `config.yaml`.

**Alternative: MSigDB R package** (no account needed)

```r
# install.packages("msigdbr")
library(msigdbr)

# Get REACTOME + KEGG for human
gsets <- rbind(
    msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME"),
    msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG_LEGACY")
)

# Write as GMT
gmt_lines <- tapply(gsets$gene_symbol, gsets$gs_name, function(genes) {
    paste(c(unique(gsets$gs_name[gsets$gene_symbol == genes[1]]),
            "na", genes), collapse = "\t")
})
# Cleaner approach:
gmt_list  <- split(gsets$gene_symbol, gsets$gs_name)
gmt_lines <- sapply(names(gmt_list), function(nm) {
    paste(c(nm, "na", gmt_list[[nm]]), collapse = "\t")
})
writeLines(gmt_lines, "resources/gsea_files/reactome_kegg.gmt.txt")
cat("GMT written:", length(gmt_lines), "gene sets\n")
```

---

## Project Layout

After setup, your directory should look like this:

```
pipeline/
├── config/
│   └── config.yaml                  ← main configuration
├── contrasts/
│   └── contrasts.csv                ← your contrasts (copy from contrasts_template.csv)
├── data/
│   ├── counts/
│   │   ├── Beta.counts.csv          ← one file per stratum
│   │   ├── Alpha.counts.csv
│   │   └── ...
│   └── sample_metadata.csv          ← one row per biosample
├── resources/
│   └── gsea_files/
│       ├── rpl_genes.csv
│       ├── rps_genes.csv
│       ├── mtr_genes.csv
│       └── reactome_kegg.gmt.txt
├── results/                         ← created automatically on first run
├── logs/                            ← created automatically on first run
└── workflow/
    └── scripts/
        ├── 01_prepare_inputs.R
        ├── 02_prepare_coldata.R
        ├── 03_run_deseq_base.R
        ├── 04_run_ruvseq.R
        ├── 05_run_deseq_final.R
        ├── 06_run_fgsea.R
        ├── 07_aggregate_results.R
        ├── fetch/
        │   ├── pankbase_config.yaml
        │   └── pankbase_helpers.R
        └── utils/
            ├── filter_utils.R
            ├── io_utils.R
            ├── logging_utils.R
            ├── plot_utils.R
            └── stats_utils.R
```

### Count matrix format

Each stratum needs its own CSV at `data/counts/{stratum_safe_name}.counts.csv`.
The safe name is your stratum name with spaces and special characters replaced
by underscores (e.g. `"CD4 T cell"` → `CD4_T_cell.counts.csv`).

Format: genes as rows, samples as columns. First column is gene symbol.

```
gene,SAMN001,SAMN002,SAMN003,...
INS,142,0,88,...
GCG,0,201,0,...
```

### Sample metadata format

`data/sample_metadata.csv` — one row per biosample, columns include at minimum:
- The column named in `biosample_id_col` in your contrasts file (e.g. `sample_accession`)
- The contrast variable (e.g. `disease_status`)
- Any covariates (e.g. `Sex`, `BMI`, `Age`)
- The donor ID column named in `filtering.donor_col` in config (e.g. `donor_id`)

---

## Configuration

1. Copy the template and edit for your project:

```bash
cp contrasts/contrasts_template.csv contrasts/contrasts.csv
```

2. Edit `config/config.yaml`:
   - Update all file paths under `inputs:` and `outputs:`
   - Set `filtering.donor_col` to your donor ID column name
   - Set thresholds as needed (`deseq2.fdr_threshold`, etc.)
   - Set `steps.ruvseq: false` to skip RUVseq on a first exploratory run

3. Edit `contrasts/contrasts.csv`:
   - One row per comparison
   - `strata` field lists cell types / tissues to loop over, pipe-separated

---

## Running Step by Step (Local)

All scripts are called with `Rscript` from the **pipeline root directory**.
For a first test, pick one contrast and one stratum.

### Step 0 — Set your working variables

```bash
# Set these for your test run
PIPELINE_ROOT="/path/to/your/pipeline"
CONFIG="config/config.yaml"
CONTRAST="T1D_vs_ND"
STRATUM="Beta"

cd "$PIPELINE_ROOT"
```

### Step 01 — Validate inputs

Runs once. Checks all files exist, all columns match, all thresholds are sane.
Fix all errors before proceeding.

```bash
Rscript workflow/scripts/01_prepare_inputs.R \
    --config "$CONFIG"
```

Expected output: `All checks passed. Ready to run.`
Creates: `logs/validation.log`

### Step 02 — Prepare coldata

Filters samples, deduplicates donors, builds coldata for this stratum.

```bash
Rscript workflow/scripts/02_prepare_coldata.R \
    --config  "$CONFIG" \
    --contrast "$CONTRAST" \
    --stratum  "$STRATUM"
```

Creates in `results/{contrast}/{stratum}/intermediates/`:
- `coldata.csv`
- `counts_filtered.csv`
- `name_map.csv`

### Step 03 — Base DESeq2

```bash
Rscript workflow/scripts/03_run_deseq_base.R \
    --config  "$CONFIG" \
    --contrast "$CONTRAST" \
    --stratum  "$STRATUM"
```

Creates:
- `results_base.tsv` — all genes, with `log2FoldChange` and `log2FoldChange_shrunk`
- `dds_base.rds`
- `plots/volcano_base.pdf`

### Step 04 — RUVseq

```bash
Rscript workflow/scripts/04_run_ruvseq.R \
    --config  "$CONFIG" \
    --contrast "$CONTRAST" \
    --stratum  "$STRATUM"
```

Creates:
- `ruvseq_best_coldata.csv` — coldata + W factor columns
- `ruvseq_summary.csv` — best k, safe Ws, disease-associated Ws
- `plots/ruvseq_diagnostics.pdf`

### Step 05 — Final DESeq2 (with W factors)

```bash
Rscript workflow/scripts/05_run_deseq_final.R \
    --config  "$CONFIG" \
    --contrast "$CONTRAST" \
    --stratum  "$STRATUM"
```

Creates:
- `results_ruv.tsv`
- `dds_ruv.rds`
- `plots/volcano_ruv.pdf`

### Step 06 — fGSEA

```bash
Rscript workflow/scripts/06_run_fgsea.R \
    --config  "$CONFIG" \
    --contrast "$CONTRAST" \
    --stratum  "$STRATUM"
```

Creates (depending on `fgsea.run_on` config):
- `fgsea_base_all.tsv`, `fgsea_base_sig.tsv`
- `fgsea_ruv_all.tsv`, `fgsea_ruv_sig.tsv`
- `plots/fgsea_base.pdf`, `plots/fgsea_ruv.pdf`

### Step 07 — Aggregate results

Runs once per contrast after all strata complete. Discovers strata from the
filesystem — no `--stratum` argument needed.

```bash
Rscript workflow/scripts/07_aggregate_results.R \
    --config  "$CONFIG" \
    --contrast "$CONTRAST"
```

Creates in `results/{contrast}/summary/`:
- `contrast_summary.csv` — one row per stratum, all DEG and pathway counts
- `log_summary.tsv` — only written if any stratum had warnings/errors

### All steps as a shell script (one contrast, one stratum)

```bash
#!/usr/bin/env bash
set -euo pipefail

CONFIG="config/config.yaml"
CONTRAST="T1D_vs_ND"
STRATUM="Beta"

echo "[$(date)] Step 01: validate"
Rscript workflow/scripts/01_prepare_inputs.R --config "$CONFIG"

echo "[$(date)] Step 02: coldata"
Rscript workflow/scripts/02_prepare_coldata.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

echo "[$(date)] Step 03: base DESeq2"
Rscript workflow/scripts/03_run_deseq_base.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

echo "[$(date)] Step 04: RUVseq"
Rscript workflow/scripts/04_run_ruvseq.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

echo "[$(date)] Step 05: final DESeq2"
Rscript workflow/scripts/05_run_deseq_final.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

echo "[$(date)] Step 06: fGSEA"
Rscript workflow/scripts/06_run_fgsea.R \
    --config "$CONFIG" --contrast "$CONTRAST" --stratum "$STRATUM"

echo "[$(date)] Step 07: aggregate"
Rscript workflow/scripts/07_aggregate_results.R \
    --config "$CONFIG" --contrast "$CONTRAST"

echo "[$(date)] Done. Check results/${CONTRAST}/summary/contrast_summary.csv"
```

Save as `run_one.sh`, make executable with `chmod +x run_one.sh`, then run
as `./run_one.sh`.

---

## Reading the Log Files

Every step 02–06 writes a structured log to `logs/{contrast_id}/{stratum}.log`:

```
[2026-02-27 14:32:01] [INFO]  contrast=T1D_vs_ND | biosample=Beta | step=start | ...
[2026-02-27 14:32:03] [INFO]  contrast=T1D_vs_ND | biosample=Beta | step=filter_samples | n_samples=47→31
[2026-02-27 14:32:15] [WARN]  contrast=T1D_vs_ND | biosample=Beta | step=lfc_shrinkage | msg=lfcShrink failed...
[2026-02-27 14:32:16] [INFO]  contrast=T1D_vs_ND | biosample=Beta | step=complete | n_degs_base=142
```

To see all warnings and errors across all strata after a run:

```bash
grep -h "WARN\|ERROR\|SKIP" logs/{contrast_id}/*.log
```

Step 07 also writes `results/{contrast}/summary/log_summary.tsv` with all
non-INFO lines collected in one place — only written if issues exist.

---

## Running with Snakemake

> Coming in a future version. The Snakemake workflow will handle dependency
> management across all contrasts and strata, with parallel execution of
> steps 02–06 per stratum.

---

## Running on SLURM

> Coming in a future version. The SLURM submission scripts will wrap steps
> 02–06 as array jobs (one job per contrast × stratum), with step 07 as a
> dependent job that waits for all stratum jobs to complete.

---

## Output Files Reference

```
results/
└── {contrast_id}/
    ├── {stratum}/
    │   ├── intermediates/
    │   │   ├── coldata.csv                  step 02
    │   │   ├── counts_filtered.csv          step 02
    │   │   ├── name_map.csv                 step 02
    │   │   ├── results_base.tsv             step 03  ← log2FoldChange + log2FoldChange_shrunk
    │   │   ├── dds_base.rds                 step 03
    │   │   ├── ruvseq_best_coldata.csv      step 04
    │   │   ├── ruvseq_summary.csv           step 04
    │   │   ├── results_ruv.tsv              step 05  ← same two-column LFC format
    │   │   ├── dds_ruv.rds                  step 05
    │   │   ├── fgsea_base_all.tsv           step 06
    │   │   ├── fgsea_base_sig.tsv           step 06
    │   │   ├── fgsea_ruv_all.tsv            step 06
    │   │   └── fgsea_ruv_sig.tsv            step 06
    │   └── plots/
    │       ├── volcano_base.pdf             step 03
    │       ├── volcano_ruv.pdf              step 05
    │       ├── ruvseq_diagnostics.pdf       step 04  ← correlation heatmap + RLE elbow
    │       ├── fgsea_base.pdf               step 06
    │       └── fgsea_ruv.pdf                step 06
    └── summary/
        ├── contrast_summary.csv             step 07  ← always written
        └── log_summary.tsv                  step 07  ← only if issues found
```

**`results_base.tsv` / `results_ruv.tsv` columns:**

| Column | Description |
|---|---|
| `gene` | Gene symbol |
| `baseMean` | Mean normalized count across all samples |
| `log2FoldChange` | DESeq2 MLE estimate |
| `log2FoldChange_shrunk` | apeglm-shrunk estimate (`NA` if shrinkage disabled or failed) |
| `lfcSE` | Standard error of log2FoldChange |
| `stat` | Wald statistic (used as fGSEA ranking statistic by default) |
| `pvalue` | Raw p-value |
| `padj` | BH-adjusted p-value |

---

## Troubleshooting

**`every gene contains at least one zero`**
DESeq2 cannot estimate size factors. Usually caused by too few samples or extreme
sparsity. Reduce `filtering.min_gene_counts` or check that you have enough
samples in both groups.

**`lfcShrink failed`**
apeglm requires at least one sample per group with non-zero counts for the
coefficient. This is a warning, not an error — `log2FoldChange_shrunk` will be
NA in the output but `log2FoldChange` is still valid.

**`no count matrix found for stratum 'X'`**
The count file expected at `data/counts/X.counts.csv` (after safe-name conversion)
does not exist. Check that your `strata` values in contrasts.csv match the actual
filenames in `data/counts/`. Spaces and special characters are converted to
underscores.

**`W pruning removed all latent factors`**
All W factors correlated with the contrast variable (raw p < 0.05). The final
DESeq2 model will run without any W factors — equivalent to the base model.
This is logged as a warning. Consider whether your sample size is sufficient
for RUVseq to find true technical variation.

**fGSEA returns 0 significant pathways**
Not necessarily an error. Check: (1) the ranking statistic produces a reasonable
spread (`fgsea.ranking_stat`), (2) the GMT file contains human gene symbols that
match your gene names, (3) `fgsea.min_gene_set_size` is not set too high.
