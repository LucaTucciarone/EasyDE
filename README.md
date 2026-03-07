# EasyDE

A modular R pipeline for pseudobulk differential expression analysis with
RUVseq correction and gene set enrichment. Designed for scRNA-seq data
([PanKbase](workflow/scripts/fetch/PANKBASE.md) and similar resources)
but works with any pseudobulk count matrices.

Orchestrated by Snakemake for local or SLURM cluster execution.

---

## Pipeline Overview

Seven sequential steps. Steps 02–06 run in parallel across every
**contrast × stratum** (cell type) combination:

```
01  Validate config + contrasts          (once per run — soft gate)
02  Filter samples, build coldata        (per contrast × stratum)
03  Base DESeq2                          (per contrast × stratum)
04  RUVseq normalization + k selection   (per contrast × stratum)  ← skipped in paired mode
05  Final DESeq2 with W factors          (per contrast × stratum)  ← skipped in paired mode
06  fGSEA pathway enrichment             (per contrast × stratum)
07  Aggregate results + status plot      (once per contrast)
```

Step 02 includes **preflight validation** — checks for design matrix rank,
zero-inflation, and data quality issues. Rejected strata are skipped cleanly
through the rest of the pipeline.

---

## Quick Start

```bash
# 1. Install and activate
micromamba env create -f installation/EasyDE_install.yml
micromamba activate EasyDE

# 2. Fetch data from PanKbase (optional — skip if you have your own data)
Rscript workflow/scripts/fetch/pankbase_helpers.R \
    --config config/config_traits.yaml \
    --pankbase-config workflow/scripts/fetch/pankbase_config.yaml \
    --strata "Beta|Alpha|Delta"

# 3. Run the full pipeline
snakemake --profile profiles/local \
    --config pipeline_config=config/config_traits.yaml

# 4. Check results
cat results/Diabetic_vs_ND/summary/contrast_summary.csv | column -t -s,
```

---

## Project Layout

```
EasyDE/
├── Snakefile                             ← Snakemake workflow (do not edit)
├── config/                               ← EDIT THESE
│   ├── config_traits.yaml                ← unpaired analysis config
│   └── config_treatments.yaml            ← paired analysis config
├── contrasts/                            ← EDIT THESE
│   ├── contrasts.csv                     ← your contrast definitions
│   └── contrasts_template.csv            ← documented template to copy from
├── profiles/                             ← Snakemake execution profiles
│   ├── local/config.yaml
│   └── slurm/config.yaml
├── data/                                 ← YOUR DATA (or fetched from PanKbase)
│   ├── counts/                           ← one CSV per cell type
│   └── sample_metadata.csv
├── resources/gsea_files/                 ← shipped with the repo
├── workflow/scripts/                     ← pipeline code (do not edit)
│   ├── 01–07_*.R                         ← pipeline steps
│   ├── fetch/                            ← PanKbase data adapter
│   └── utils/                            ← shared R utilities
├── notebooks/                            ← post-run exploration
│   └── explore_results.ipynb
├── docs/                                 ← detailed documentation
└── installation/
    └── EasyDE_install.yml
```

| What to edit | Where |
|-------------|-------|
| Analysis parameters | `config/config_*.yaml` |
| Contrast definitions | `contrasts/contrasts.csv` |
| Execution settings | `profiles/*/config.yaml` |
| Data files | `data/counts/` + `data/sample_metadata.csv` |

---

## Documentation

| Guide | Contents |
|-------|----------|
| **[Installation](docs/INSTALLATION.md)** | Environment setup, verification, known issues |
| **[Configuration](docs/CONFIGURATION.md)** | Config files, contrast definitions, data formats |
| **[Running](docs/RUNNING.md)** | Snakemake commands, SLURM setup, step-by-step manual run |
| **[Methods](docs/METHODS.md)** | DESeq2, RUVSeq, fGSEA, paired analysis, preflight validation |
| **[Output](docs/OUTPUT.md)** | File tree, column descriptions, status labels |
| **[Troubleshooting](docs/TROUBLESHOOTING.md)** | Common errors, log reading, resource files |

---

## Pipeline Status Values

The `contrast_summary.csv` produced by step 07 categorizes each stratum:

| Status | Meaning |
|--------|---------|
| `success` | All steps completed, DE results produced |
| `skipped_preflight` | Preflight rejected (rank-deficient design, zero-inflated counts) |
| `skipped_filtering` | Sample filtering eliminated all samples (no paired donors, too few samples) |
| `error` | Unexpected failure after preflight passed |

Step 07 also generates a **status tile plot** (`status_tile.pdf`) for quick
visual overview of which strata succeeded or were skipped.

---

## Exploratory Notebook

After running the pipeline, use the Jupyter notebook to visualize results:

```bash
jupyter notebook notebooks/explore_results.ipynb
```

Produces DE boxplots, fGSEA heatmaps, divergent pathway analysis,
and pipeline status overviews across all contrasts and cell types.
