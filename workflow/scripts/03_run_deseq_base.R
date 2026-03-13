# =============================================================================
# 03_run_deseq_base.R
# EasyDE - Base DESeq2 Analysis
#
# PURPOSE:
#   Runs the initial DESeq2 model on the filtered coldata and count matrix
#   produced by 02_prepare_coldata.R. This is the "baseline" model - no
#   RUVseq correction yet. Results are used directly if ruvseq is disabled,
#   and as the empirical gene set input for RUVseq if it is enabled.
#
# INPUTS:
#   - results/{contrast_id}/{stratum}/intermediates/coldata.csv
#   - results/{contrast_id}/{stratum}/intermediates/counts_filtered.csv
#   - results/{contrast_id}/{stratum}/intermediates/name_map.csv
#
# OUTPUTS:
#   - intermediates/dds_base.rds             DESeqDataSet object
#   - intermediates/results_base.tsv         Full DESeq2 results table
#   - intermediates/results_base_shrunk.tsv  LFC-shrunk results (if enabled)
#   - plots/volcano_base.pdf                 Volcano plot
#
# USAGE:
#   Rscript workflow/scripts/03_run_deseq_base.R \
#       --config config/config.yaml \
#       --contrast T1D_vs_ND \
#       --stratum Beta
# =============================================================================

suppressMessages({
    library(yaml)
    library(data.table)
    library(dplyr)
    library(DESeq2)
    library(optparse)
})

# Resolve utils directory relative to this script
.script_dir <- tryCatch({
    args     <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        dirname(normalizePath(sub("--file=", "", file_arg[1])))
    } else {
        dirname(normalizePath(sys.frame(1)$ofile))
    }
}, error = function(e) getwd())


source(file.path(.script_dir, "utils", "io_utils.R"))
source(file.path(.script_dir, "utils", "logging_utils.R"))
source(file.path(.script_dir, "utils", "stats_utils.R"))
source(file.path(.script_dir, "utils", "plot_utils.R"))


# =============================================================================
# Script-specific helpers
# =============================================================================

#' Build the DESeq2 model formula string
#'
#' Constructs the formula from the active covariates and contrast variable.
#' Covariates that were dropped (single-value) are excluded.
#' Format: ~ covariate_1 + covariate_2 + contrast_var
#'
#' @param covariates   Character vector of active sanitized covariate names
#' @param contrast_var Sanitized contrast variable name
#' @return Formula string (e.g. "~ Sex + BMI + disease_status")
build_deseq_formula <- function(covariates, contrast_var) {
    terms <- c(covariates, contrast_var)
    paste0("~ ", paste(terms, collapse = " + "))
}


#' Extract DESeq2 results for a contrast or continuous variable
#'
#' Handles both categorical (contrast vector) and numeric (name-based) contrasts.
#' If lfc_shrinkage is enabled, runs lfcShrink (apeglm) and merges the shrunk
#' estimates as an additional column `log2FoldChange_shrunk` in the same output
#' data.frame. This produces ONE file with both LFC estimates side by side,
#' rather than two separate files with ambiguous column names.
#' pvalue, padj, and stat are always from the regular (unshrunk) results.
#'
#' @param dds              DESeqDataSet after running DESeq()
#' @param contrast_var     Sanitized contrast variable name
#' @param control_grp      Control group value (or "" for numeric)
#' @param experimental_grp Experimental group value (or "" for numeric)
#' @param lfc_shrinkage    Logical: attempt lfcShrink and add shrunk column?
#' @param fdr_threshold    Numeric: FDR cutoff for counting DEGs
#' @param l2fc_threshold   Numeric: minimum absolute LFC for counting DEGs
#' @param logger           Logger instance (optional, for shrinkage warnings)
#' @return Named list: $results (data.frame with log2FoldChange and
#'         log2FoldChange_shrunk), $n_degs (integer), $formula (string),
#'         $shrinkage_ran (logical)
extract_deseq_results <- function(dds,
                                  contrast_var,
                                  control_grp,
                                  experimental_grp,
                                  lfc_shrinkage  = TRUE,
                                  fdr_threshold  = 0.05,
                                  l2fc_threshold = 0,
                                  logger         = NULL) {

    is_categorical <- !is.na(control_grp) && control_grp != ""

    if (is_categorical) {
        result <- DESeq2::results(dds,
                                  contrast = c(contrast_var, experimental_grp, control_grp))
    } else {
        result <- DESeq2::results(dds, name = contrast_var)
    }

    result_df        <- as.data.frame(result)
    result_df$gene   <- rownames(result_df)
    result_df        <- result_df[order(result_df$pvalue, na.last = TRUE), ]

    # Count DEGs using both thresholds
    n_degs <- sum(
        !is.na(result_df$padj) &
        result_df$padj < fdr_threshold &
        abs(result_df$log2FoldChange) > l2fc_threshold,
        na.rm = TRUE
    )

    # LFC shrinkage: merged as a second column, not a separate file.
    # log2FoldChange    = regular DESeq2 MLE estimate (always present)
    # log2FoldChange_shrunk = apeglm-shrunk estimate (NA if shrinkage disabled/failed)
    shrinkage_ran <- FALSE
    result_df$log2FoldChange_shrunk <- NA_real_

    if (lfc_shrinkage && is_categorical) {
        # DESeq2 uses make.names() internally to sanitize variable/level names.
        # Replicate that exactly to match resultsNames(dds).
        coef_name <- paste0(
            make.names(contrast_var), "_",
            make.names(experimental_grp), "_vs_",
            make.names(control_grp)
        )
        available_coefs <- DESeq2::resultsNames(dds)
        if (!coef_name %in% available_coefs) {
            # Fallback: grep for the experimental group name
            match <- available_coefs[grepl(make.names(experimental_grp), available_coefs, fixed = TRUE)]
            if (length(match) == 1) coef_name <- match
        }
        result_shrunk <- tryCatch(
            DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE),
            error = function(e) {
                if (!is.null(logger)) {
                    logger$warn("lfc_shrinkage",
                        sprintf("lfcShrink failed (coef=%s): %s", coef_name, e$message))
                }
                NULL
            }
        )
        if (!is.null(result_shrunk)) {
            shrunk_df <- as.data.frame(result_shrunk)
            idx <- match(result_df$gene, rownames(shrunk_df))
            result_df$log2FoldChange_shrunk <- shrunk_df$log2FoldChange[idx]
            shrinkage_ran <- TRUE
        }
    }

    list(
        results       = result_df,
        n_degs        = n_degs,
        formula       = deparse(DESeq2::design(dds)),
        shrinkage_ran = shrinkage_ran
    )
}


# =============================================================================
# Main
# =============================================================================

main <- function(config_path, contrast_id, stratum) {

    loaded       <- load_config(config_path)
    cfg          <- loaded$cfg
    contrasts_df <- loaded$contrasts_df
    contrast_row <- contrasts_df[contrasts_df$contrast_id == contrast_id, ]

    stratum_safe <- sanitize_names(stratum)
    log_file     <- file.path(cfg$outputs$logs_dir, contrast_id,
                               paste0(stratum_safe, ".log"))
    logger       <- make_logger(contrast_id = contrast_id,
                                biosample   = stratum,
                                log_file    = log_file,
                                level       = cfg$logging$level)

    logger$info("start", "script=03_run_deseq_base")

    # Define inter_dir and finals_dir early and register skip sentinels
    inter_dir  <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "intermediates")
    finals_dir <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "finals")
    register_skip_sentinel(finals_dir, c("results_base.tsv"))
    register_skip_sentinel(inter_dir,  c("model_info_base.csv"))

    # --- Sentinel fires on any skip ---
    on.exit({
        for (.f in c("results_base.tsv")) {
            .p <- file.path(finals_dir, .f)
            if (!file.exists(.p)) {
                dir.create(dirname(.p), recursive = TRUE, showWarnings = FALSE)
                writeLines("skipped=TRUE", .p)
            }
        }
        for (.f in c("model_info_base.csv")) {
            .p <- file.path(inter_dir, .f)
            if (!file.exists(.p)) {
                dir.create(dirname(.p), recursive = TRUE, showWarnings = FALSE)
                writeLines("skipped=TRUE", .p)
            }
        }
    }, add = TRUE)


    # --- Load intermediates from 02 ---

    # Bail out gracefully if upstream step was skipped
    if (is_skip_sentinel(file.path(inter_dir, "coldata.csv"))) {
        logger$skip("upstream_skip", "upstream step was skipped — propagating skip")
        return(invisible(NULL))
    }

    # Check preflight status from 02
    preflight_status <- read_preflight_status(inter_dir)
    if (preflight_status == "fail") {
        logger$skip("preflight_fail", "preflight validation failed in step 02 — propagating skip")
        return(invisible(NULL))
    }
    if (preflight_status == "warn") {
        logger$warn("preflight_status", "preflight reported warnings — proceeding with caution")
    }

    coldata  <- try_logged(
        read.csv(file.path(inter_dir, "coldata.csv"), check.names = FALSE, row.names = NULL),
        logger = logger, step = "load_coldata", success = "coldata loaded"
    )
    if (is.null(coldata)) return(invisible(NULL))
    rownames(coldata) <- coldata[[1]]

    counts_df <- try_logged(
        fread(file.path(inter_dir, "counts_filtered.csv")) %>% as.data.frame(),
        logger = logger, step = "load_counts", success = "counts loaded"
    )
    if (is.null(counts_df)) return(invisible(NULL))
    gene_names         <- counts_df[[1]]
    counts_mat         <- as.matrix(counts_df[, -1])
    storage.mode(counts_mat) <- "integer"
    rownames(counts_mat)     <- gene_names

    name_map_df <- fread(file.path(inter_dir, "name_map.csv")) %>% as.data.frame()
    name_map    <- setNames(name_map_df$sanitized, name_map_df$original)

    # --- Resolve sanitized names ---
    contrast_var     <- trimws(contrast_row$contrast_var)
    control_grp      <- trimws(contrast_row$control_grp)
    experimental_grp <- trimws(contrast_row$experimental_grp)
    covariates_raw   <- parse_pipe_field(contrast_row$covariates)

    safe_contrast <- name_map[[contrast_var]]
    safe_covs     <- sanitize_names(covariates_raw)

    # Only keep covariates that survived into coldata (not dropped for single value)
    safe_covs <- safe_covs[safe_covs %in% colnames(coldata)]

    # --- Build formula and run DESeq2 ---
    formula_str <- build_deseq_formula(safe_covs, safe_contrast)
    logger$info("deseq_base", sprintf("formula: %s", formula_str))

    # --- Check preflight for zero-inflation flag ---
    preflight <- read_preflight_report(inter_dir)
    use_poscounts <- !is.null(preflight) &&
        "zero_inflated" %in% colnames(preflight) &&
        isTRUE(preflight$zero_inflated[1])

    # Align column order
    counts_mat <- counts_mat[, rownames(coldata), drop = FALSE]

    dds <- try_logged({
        dds <- suppressWarnings(suppressMessages(DESeq2::DESeqDataSetFromMatrix(
            countData = counts_mat,
            colData   = coldata,
            design    = as.formula(formula_str)
        )))
        if (use_poscounts) {
            logger$info("deseq_base", "zero-inflated counts — using sfType='poscounts'")
            suppressMessages(DESeq2::DESeq(dds, sfType = "poscounts"))
        } else {
            suppressMessages(DESeq2::DESeq(dds))
        }
    },
    logger  = logger,
    step    = "deseq_base",
    success = "DESeq2 completed"
    )
    if (is.null(dds)) return(invisible(NULL))

    logger$info("deseq_base", "extracting results")

    # --- Extract results ---
    deseq_out <- try_logged(
        extract_deseq_results(
            dds              = dds,
            contrast_var     = safe_contrast,
            control_grp      = control_grp,
            experimental_grp = experimental_grp,
            lfc_shrinkage    = isTRUE(cfg$deseq2$lfc_shrinkage),
            fdr_threshold    = cfg$deseq2$fdr_threshold,
            l2fc_threshold   = cfg$deseq2$l2fc_threshold,
            logger           = logger
        ),
        logger  = logger,
        step    = "extract_results",
        success = "results extracted"
    )
    if (is.null(deseq_out)) return(invisible(NULL))

    logger$info("deseq_base",
        sprintf("DEGs (FDR<%.2f, |LFC|>%.1f): %d",
                cfg$deseq2$fdr_threshold, cfg$deseq2$l2fc_threshold,
                deseq_out$n_degs))
    if (deseq_out$shrinkage_ran) {
        logger$info("deseq_base", "log2FoldChange_shrunk column added to results")
    }

    # --- Write outputs ---
    out_dir   <- inter_dir
    plot_dir  <- file.path(cfg$outputs$results_dir, contrast_id, stratum_safe, "plots")
    dir.create(finals_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

    if (isTRUE(cfg$write_outputs$deseq_base_results)) {
        write_result(deseq_out$results,
                     file.path(finals_dir, "results_base.tsv"))
        logger$info("write_outputs",
            sprintf("results_base.tsv written (shrunk_col=%s)",
                    ifelse(deseq_out$shrinkage_ran, "yes", "no")))
    }

    # Save dds object for downstream RUVseq step
    saveRDS(dds, file.path(out_dir, "dds_base.rds"))
    logger$info("write_outputs", "dds_base.rds saved")

    # Write lightweight model info for downstream aggregation (avoids loading .rds)
    model_info_base <- data.frame(
        formula = formula_str,
        n_degs  = deseq_out$n_degs,
        stringsAsFactors = FALSE
    )
    write_result(model_info_base, file.path(out_dir, "model_info_base.csv"), sep = ",")
    logger$info("write_outputs", "model_info_base.csv written")

    if (isTRUE(cfg$write_outputs$volcano_plots)) {
        # For volcano: prefer shrunk LFC for x-axis if available, else regular LFC
        plot_df <- deseq_out$results
        if (deseq_out$shrinkage_ran) {
            plot_df$log2FoldChange <- plot_df$log2FoldChange_shrunk
        }

        p <- plot_volcano(
            deseq_df = plot_df,
            fdr      = cfg$deseq2$fdr_threshold,
            title    = formula_str,
            subtitle = sprintf("DEGs: %d | %s - %s vs %s",
                               deseq_out$n_degs, stratum, experimental_grp, control_grp)
        )
        pdf(file.path(plot_dir, "volcano_base.pdf"), width = 8, height = 6)
        print(p)
        dev.off()
        logger$info("write_outputs", "volcano_base.pdf written")
    }

    logger$info("complete",
        sprintf("n_degs_base=%d | formula=%s", deseq_out$n_degs, formula_str))

    invisible(list(dds = dds, results = deseq_out))
}


# =============================================================================
# CLI entry point
# =============================================================================

if (!interactive() && identical(sys.nframe(), 0L)) {

    opt_list <- list(
        make_option("--config",   type = "character", default = "config/config.yaml"),
        make_option("--contrast", type = "character"),
        make_option("--stratum",  type = "character")
    )
    opts <- parse_args(OptionParser(option_list = opt_list))
    if (is.null(opts$contrast)) stop("--contrast is required")
    if (is.null(opts$stratum))  stop("--stratum is required")

    main(opts$config, opts$contrast, opts$stratum)
}