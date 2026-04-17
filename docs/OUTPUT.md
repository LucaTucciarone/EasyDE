# Output Files

Results can be written to any directory by setting `outputs.results_dir` in
your config file. The structure within that directory is:

```
results/
├── pipeline_status_overview.pdf            step 09 — cross-contrast status heatmap
├── pipeline_summary.csv                    step 09 — one row per contrast × stratum
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
    │   ├── plots/
    │   │   ├── volcano_base.pdf                 step 03
    │   │   ├── volcano_ruv.pdf                  step 05
    │   │   ├── ruvseq_diagnostics.pdf           step 04 — correlation + elbow plot
    │   │   ├── fgsea_base.pdf                   step 06
    │   │   └── fgsea_ruv.pdf                    step 06
    │   │
    │   └── benchmarking/                        step 07 — per-stratum signature benchmarking
    │       ├── benchmark_signatures.tsv         FORA + fGSEA results vs curated signatures
    │       ├── benchmark_positive.pdf           positive control heatmap (per-stratum)
    │       ├── benchmark_negative.pdf           negative control heatmap (per-stratum)
    │       ├── benchmark_collapsed.pdf          class-level collapsed heatmap (per-stratum)
    │       └── benchmark.done                   Snakemake touch file
    │
    └── summary/                                 per-contrast aggregation (step 08)
        ├── contrast_summary.csv                 one row per stratum, all counts + status
        ├── preflight_summary.csv                preflight results across strata
        ├── status_overview.pdf                  status heatmap (cell type × status + DEG counts)
        ├── {contrast}_deseq_final.tsv           best-model results (RUV preferred, Base fallback)
        ├── {contrast}_deseq_base.tsv            merged base DE (all strata)
        ├── {contrast}_deseq_ruv.tsv             merged RUV DE (all strata)
        ├── {contrast}_fgsea_base.tsv            merged base fGSEA (all strata)
        ├── {contrast}_fgsea_ruv.tsv             merged RUV fGSEA (all strata)
        ├── {contrast}_fgsea_final.tsv           best-model fGSEA (RUV preferred, Base fallback)
        ├── positive_benchmarking.html           interactive pathway drilldown (fGSEA × mega-sets)
        ├── negative_benchmarking.html           interactive negative controls drilldown (LLM signatures)
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

### Best-model consolidated file

`{contrast}_deseq_final.tsv` contains the same columns as above, plus:

| Column | Description |
|--------|-------------|
| `model` | Which model produced this row: `"RUVseq"` or `"Base"` |
| `formula` | The exact design formula used (e.g. `~ Age + Gender + W2 + disease_status`) |

This file selects RUVseq results where available (`success` strata) and falls
back to Base DESeq2 for `success_base_only` strata. It is the recommended
input for downstream analyses that need a single best-effort result per gene
per cell type.

### fGSEA best-model consolidated file

`{contrast}_fgsea_final.tsv` mirrors the `_deseq_final.tsv` logic for pathway results: RUVseq fGSEA where available, base fGSEA for strata where only base ran. Same schema as `_fgsea_base.tsv` / `_fgsea_ruv.tsv`. Used as the data source for `positive_benchmarking.html`.

---

## Interactive HTML Outputs (step 08)

Two self-contained HTML files are generated per contrast (requires `pathway_drilldown.enabled: true` in config):

### `positive_benchmarking.html`

Three-level interactive drill-down for fGSEA pathway results across all cell types:

- **Level 1 (L1):** Mega-sets — curated groupings of Reactome/KEGG/Hallmark pathways into biological concepts (e.g. `ER_STRESS_UPR`, `APOPTOSIS`, `INSULIN_SECRETION`). One row per mega-set; NES columns per cell type, colored by direction and magnitude.
- **Level 2 (L2):** Member pathways within a mega-set. Expand any L1 row to see individual pathway NES scores.
- **Level 3 (L3):** Leading-edge genes for a specific pathway × cell-type cell. Expand any L2 cell to see gene-level log2FC and padj.

NES cells are colored red (positive) / blue (negative) by default. "Highlight significant" toggle fades non-significant cells and pops significant ones.

Mega-set definitions come from `resources/pathway_drilldown/mega_sets.tsv` (configurable).

### `negative_benchmarking.html`

Same three-level structure, but uses LLM-curated artifact and confounder gene signatures from `resources/benchmarking/master_gene_signatures.tsv` (negative control rows only):

- **L1:** Signature class (e.g. `dissociation_stress`, `sex`, `BMI`, `smoking`). fGSEA run on the union of all class genes.
- **L2:** Individual signatures within the class.
- **L3:** Leading-edge genes with DE stats.

Color scheme: artifact classes in red, subject_confounder classes in purple. Significant highlight is ON by default (unexpected enrichment warrants attention).

---

## Pipeline Status Values

The `status` column in `contrast_summary.csv` uses a fine-grained taxonomy:

| Status | Meaning | When |
|--------|---------|------|
| `success` | RUVseq model ran successfully | Both base and RUV results produced |
| `success_base_only` | Base DESeq2 succeeded, RUV skipped/failed | Paired mode, or RUV crashed (e.g. singular fit), or all W factors excluded |
| `skipped_no_samples` | Zero samples after filtering | Sample filters or group matching eliminated all samples |
| `skipped_min_group` | Too few samples per group | Fewer than 3 control or experimental samples |
| `skipped_no_pairs` | No valid donor pairs | Paired contrasts: no donors present in both conditions |
| `skipped_preflight` | Design matrix rejected | Rank-deficient design, or 100% of genes zero-inflated |
| `skipped` | Other skip reason | Unclassified — check `error_message` for details |
| `error` | Analysis crashed | DESeq2 or downstream step threw an error — `error_message` populated |
| `not_run` | No data found | No log entries or output files for this stratum |

The `error_message` column is always populated for non-success statuses,
extracted from structured log entries by `extract_log_message()`.

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
