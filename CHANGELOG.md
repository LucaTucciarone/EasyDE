# Changelog

All notable changes to EasyDE are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versions follow [Semantic Versioning](https://semver.org/).

---

## [v1.1.0] — 2026-04-20

### Added

- **NA-covariate handling — single authority pattern.**
  Step 02 (`02_prepare_coldata.R`) is now the sole decision-maker for which
  covariates enter the model. Any covariate that has one or more missing
  values (`NA`) in the filtered sample set is dropped from the formula
  automatically, with a warning logged and the reason recorded in
  `preflight.csv`. The decision is written to a new intermediate file,
  `active_covariates.csv`, which steps 03, 04, and 05 read directly instead
  of re-parsing `contrasts.csv`. This replaces the previous approach that let
  DESeq2 crash with "variables in design formula cannot contain NA".

- **`active_covariates.csv` intermediate file.**
  Written by step 02 to `intermediates/active_covariates.csv`; contains two
  columns: `original_name` (as entered in `contrasts.csv`) and
  `sanitized_name` (R-safe version used in the formula). Steps 03–05 read
  this file to build their design formulas, eliminating redundant covariate
  parsing across pipeline steps.

- **`drop_reasons` tracking in preflight.**
  The `dropped_covariates` field in `preflight.csv` now records structured
  reason strings (e.g., `Age_dropped_NA_3_of_47`, `Sex_dropped_single_value`)
  rather than bare covariate names, making it easier to diagnose why a
  covariate was excluded.

- **NA-covariate edge case documented** in `docs/METHODS.md` under
  *Edge Cases and Workarounds*.

### Fixed

- **`anova()` size-mismatch crash in `correlate_latent_vars()`** (RUVseq step 04).
  When a covariate used for W-factor correlation testing had missing values,
  `lm()` with `na.action = na.omit` could produce models fitted on different
  subsets, causing `anova()` to error with "models were not all fitted to the
  same size of dataset". Fixed by pre-filtering `meta` to non-NA rows for
  each variable before fitting both the null and full models.

### Changed

- **Snakefile dependency graph** updated: `active_covariates.csv` is declared
  as an output of `prepare_coldata` and an explicit input of `deseq_base`,
  `ruvseq`, and `deseq_final`, enforcing correct Snakemake ordering.

- **Documentation — Sample and Gene Filtering** (`docs/METHODS.md`): new
  section describing QC filters (min cells per sample, donor deduplication,
  sample-size gates) and analytical restrictions (contrast group restriction,
  metadata filters, paired-design matching) applied in step 02.

---

## [v1.0.0] — 2026-04-03

Initial stable release. Features:

- Nine-step Snakemake pipeline: validate → coldata → DESeq2 base → RUVseq →
  DESeq2 final → fGSEA → benchmark → aggregate → summary.
- Unpaired and paired (donor-matched) contrast modes.
- Per-stratum preflight validation: rank checks, zero-inflation detection,
  sample-size gates, and design matrix diagnostics.
- RUVseq W-factor selection with automated safe-W identification (excludes Ws
  correlated with the contrast variable at p < 0.05).
- LFC shrinkage via apeglm with fallback.
- fGSEA on base and/or RUV results; configurable gene set files.
- Step 07 optional benchmarking against curated gene signatures.
- Step 08 per-contrast aggregation, HTML pathway drilldowns (positive and
  negative benchmarking), and status-overview heatmap.
- Step 09 pipeline-level status heatmap across all contrasts × strata.
- PanKbase data-fetching adapter (`workflow/scripts/fetch/`).
- SLURM and local Snakemake execution profiles.
