# Methods

## DESeq2 (steps 03 and 05)

Standard DESeq2 pseudobulk analysis. Step 03 fits the base model with
user-specified covariates. If RUVseq is enabled (unpaired mode only), step 05
re-fits the model with additional W latent factors to control for unwanted
technical variation. Both steps produce log2 fold changes (MLE and
apeglm-shrunk) and BH-adjusted p-values.

For **continuous** contrasts (e.g. Age, BMI, HbA1c), the pipeline detects
that `control_grp` and `experimental_grp` are blank and treats the contrast
variable as numeric. Log2 fold changes represent the effect per unit increase.

---

## Paired Analysis

When `paired_analysis: true` is set in the config, the pipeline:

1. **Matches samples by donor** across treatment and control conditions
   using the `donor_col` specified in the config.
2. **Includes donor as a blocking factor** in the DESeq2 formula
   (e.g. `~ donor_accession + treatment`), which accounts for individual
   variability.
3. **Skips RUVseq (steps 04 and 05)** because the donor blocking factor
   already captures inter-individual variation that RUVseq's latent
   factors would model. Running RUVseq on top of donor blocking is
   redundant and can absorb genuine biological signal, reducing the
   number of detected DEGs.

The paired-mode pipeline flow is: **01 → 02 → 03 → 06 → 07**.

---

## RUVseq (step 04)

> **Note:** RUVseq is automatically skipped in paired analysis mode.

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

---

## fGSEA (step 06)

Fast Gene Set Enrichment Analysis using the Wald test statistic as ranking
metric. Ribosomal (RPL/RPS) and mitochondrial (MT-) genes are excluded
before ranking to prevent these highly variable housekeeping genes from
dominating pathway results.

Runs on base results, RUV-corrected results, or both (set `fgsea.run_on`
in config). In paired mode, fGSEA runs on base results only (since RUVseq
is skipped). Default gene set: merged Reactome + KEGG pathways from MSigDB.

---

## Preflight Validation (steps 01–02)

### Step 01 — Structural validation

Checks that all files exist, required columns are present, contrasts are
well-formed, and groups have matching samples. Per-contrast issues are
**flagged but not fatal** — only structural problems (missing CSV columns,
duplicate contrast IDs) halt the pipeline.

### Step 02 — Data-level validation

For each stratum, step 02 performs:

1. **Count matrix checks** — no negative values, warns on all-zero gene rows
2. **Zero-inflation detection** — if 100% of genes have ≥1 zero across
   samples, the stratum is rejected (no DESeq2 normalization method can
   handle this extreme sparsity)
3. **Design matrix rank check** — detects rank-deficient models
4. **Near-singularity detection** — warns on covariates with >95% one level

Results are written to `intermediates/preflight.csv` with columns including
`zero_inflated`, `pct_genes_with_zero`, `rank_ok`, and `near_singular_vars`.
Steps 03–06 read the preflight status and skip rejected strata gracefully.
