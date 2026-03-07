# Output Files

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
        ├── status_tile.pdf                      status heatmap (cell type × status)
        ├── {contrast}_results_base.tsv          merged base DE (all strata)
        ├── {contrast}_results_ruv.tsv           merged RUV DE (all strata)
        ├── {contrast}_fgsea_base.tsv            merged base fGSEA (all strata)
        ├── {contrast}_fgsea_ruv.tsv             merged RUV fGSEA (all strata)
        ├── errors.log                           errors only (if any)
        ├── warnings.log                         warnings only (if any)
        └── log_summary.tsv                      all flagged log lines
```

> **Paired mode:** In paired analysis, RUVseq outputs (steps 04/05) are
> replaced by skip sentinels. The `results_ruv.tsv`, `fgsea_ruv.tsv`, and
> RUVseq intermediates will not contain results.

---

## DE Result Columns

`results_base.tsv` and `results_ruv.tsv` share this schema:

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

---

## Pipeline Status Values

The `status` column in `contrast_summary.csv` uses these values:

| Status | Meaning | When |
|--------|---------|------|
| `success` | All steps completed, DE results produced | Pipeline ran to completion |
| `skipped_preflight` | Preflight validation rejected this stratum | Design matrix rank-deficient, or data too sparse (100% of genes zero-inflated) |
| `skipped_filtering` | Sample filtering eliminated all samples | No paired donors found, too few samples, or no matching groups |
| `error` | Unexpected failure during analysis | Coldata was produced but DESeq2 or downstream steps crashed |

---

## Preflight Report Columns

The per-stratum `intermediates/preflight.csv` contains:

| Column | Description |
|--------|-------------|
| `status` | `pass`, `warn`, or `fail` |
| `n_genes` | Number of genes after filtering |
| `n_samples` | Total samples |
| `n_control` / `n_experimental` | Group sizes |
| `rank_ok` | `TRUE` if design matrix is full rank |
| `near_singular_vars` | Covariates with >95% one level |
| `zero_inflated` | `TRUE` if every gene has ≥1 zero (DESeq2 cannot normalize) |
| `pct_genes_with_zero` | Percentage of genes with at least one zero |
| `errors` / `warnings` | Pipe-separated diagnostic messages |
