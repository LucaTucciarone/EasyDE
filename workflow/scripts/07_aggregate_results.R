# =============================================================================
# 07_aggregate_results.R
# EasyDE -- Results Aggregation and Run Summary
#
# PURPOSE:
#   Runs ONCE per contrast (not per stratum). Walks all stratum subdirectories
#   under results/{contrast_id}/, collects per-stratum outputs, and assembles:
#
#     1. contrast_summary.csv       -- one row per stratum, audit trail
#     2. preflight_summary.csv      -- concatenated preflight reports
#     3. {contrast}_deseq_base.tsv  -- all strata DESeq2 base results
#     4. {contrast}_deseq_ruv.tsv   -- all strata DESeq2 RUV results
#     5. {contrast}_fgsea_base.tsv  -- all strata fGSEA base results
#     6. {contrast}_fgsea_ruv.tsv   -- all strata fGSEA RUV results
#     7. errors.log / warnings.log  -- split consolidated logs
#     8. log_summary.tsv            -- aggregated flagged log lines
#
# INPUTS:
#   - results/{contrast_id}/{stratum}/finals/*.tsv
#   - results/{contrast_id}/{stratum}/intermediates/*.csv
#   - logs/{contrast_id}/*.log
#
# USAGE:
#   Rscript workflow/scripts/07_aggregate_results.R \
#       --config config/config.yaml \
#       --contrast T1D_vs_ND
# =============================================================================

suppressMessages({
    library(yaml)
    library(data.table)
    library(dplyr)
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


# =============================================================================
# Script-level helpers
# =============================================================================

#' Safe file reader -- returns NULL if file missing or is a skip sentinel
safe_read <- function(path) {
    if (!file.exists(path)) return(NULL)
    if (is_skip_sentinel(path)) return(NULL)
    tryCatch(fread(path) %>% as.data.frame(), error = function(e) NULL)
}


#' Collect summary statistics for one stratum
#'
#' Reads whatever output files exist -- missing files just become NA values.
#' This means a partially-completed stratum still contributes a row to the
#' summary, showing exactly how far it got before failing.
#'
#' @param inter_dir     Path to {contrast_id}/{stratum}/intermediates/
#' @param finals_dir    Path to {contrast_id}/{stratum}/finals/
#' @param contrast_id   Contrast identifier string
#' @param stratum       Stratum name
#' @param fdr_threshold FDR threshold for counting DEGs and pathways
#' @param l2fc_threshold Minimum |LFC| for counting DEGs
#' @param fgsea_fdr     FDR threshold for counting significant pathways
#' @return Single-row data.frame (a summary row)
collect_stratum_summary <- function(inter_dir,
                                     finals_dir,
                                     contrast_id,
                                     stratum,
                                     fdr_threshold,
                                     l2fc_threshold,
                                     fgsea_fdr) {

    # --- coldata: sample counts (from intermediates) ---
    coldata      <- safe_read(file.path(inter_dir, "coldata.csv"))
    n_total      <- if (!is.null(coldata)) nrow(coldata) else NA

    # --- base DESeq2 results (from finals) ---
    base_results <- safe_read(file.path(finals_dir, "results_base.tsv"))

    n_degs_base <- if (!is.null(base_results)) {
        sum(!is.na(base_results$padj) &
            base_results$padj < fdr_threshold &
            abs(base_results$log2FoldChange) > l2fc_threshold, na.rm = TRUE)
    } else NA

    # --- model info: formulas (from intermediates) ---
    model_base   <- safe_read(file.path(inter_dir, "model_info_base.csv"))
    base_formula <- if (!is.null(model_base) && "formula" %in% colnames(model_base)) {
        model_base$formula[1]
    } else NA

    model_ruv    <- safe_read(file.path(inter_dir, "model_info_ruv.csv"))
    ruv_formula  <- if (!is.null(model_ruv) && "formula" %in% colnames(model_ruv)) {
        model_ruv$formula[1]
    } else NA

    # --- RUVseq summary (from intermediates) ---
    ruv_summary  <- safe_read(file.path(inter_dir, "ruvseq_summary.csv"))
    best_k       <- if (!is.null(ruv_summary)) ruv_summary$best_k[1]    else NA
    safe_ws      <- if (!is.null(ruv_summary)) ruv_summary$safe_Ws[1]   else NA
    disease_ws   <- if (!is.null(ruv_summary)) ruv_summary$disease_associated_Ws[1] else NA

    # --- RUV DESeq2 results (from finals) ---
    ruv_results  <- safe_read(file.path(finals_dir, "results_ruv.tsv"))

    n_degs_ruv <- if (!is.null(ruv_results)) {
        sum(!is.na(ruv_results$padj) &
            ruv_results$padj < fdr_threshold &
            abs(ruv_results$log2FoldChange) > l2fc_threshold, na.rm = TRUE)
    } else NA

    # --- fGSEA results (from finals -- count significant from full table) ---
    fgsea_base <- safe_read(file.path(finals_dir, "fgsea_base.tsv"))
    fgsea_ruv  <- safe_read(file.path(finals_dir, "fgsea_ruv.tsv"))

    n_pathways_base <- if (!is.null(fgsea_base)) {
        sum(!is.na(fgsea_base$padj) & fgsea_base$padj < fgsea_fdr, na.rm = TRUE)
    } else NA
    n_pathways_ruv <- if (!is.null(fgsea_ruv)) {
        sum(!is.na(fgsea_ruv$padj) & fgsea_ruv$padj < fgsea_fdr, na.rm = TRUE)
    } else NA

    # --- Read preflight report (from intermediates) ---
    preflight      <- read_preflight_report(inter_dir)
    pf_status      <- if (!is.null(preflight) && "status" %in% colnames(preflight)) preflight$status[1] else NA
    pf_dropped     <- if (!is.null(preflight) && "dropped_covariates" %in% colnames(preflight)) preflight$dropped_covariates[1] else NA
    pf_rank_ok     <- if (!is.null(preflight) && "rank_ok" %in% colnames(preflight)) preflight$rank_ok[1] else NA
    pf_near_sing   <- if (!is.null(preflight) && "near_singular_vars" %in% colnames(preflight)) preflight$near_singular_vars[1] else NA
    pf_warnings    <- if (!is.null(preflight) && "warnings" %in% colnames(preflight)) preflight$warnings[1] else NA
    pf_n_genes     <- if (!is.null(preflight) && "n_genes" %in% colnames(preflight)) preflight$n_genes[1] else NA
    pf_n_control   <- if (!is.null(preflight) && "n_control" %in% colnames(preflight)) preflight$n_control[1] else NA
    pf_n_experimental <- if (!is.null(preflight) && "n_experimental" %in% colnames(preflight)) preflight$n_experimental[1] else NA

    # --- Determine status ---
    status <- if (!is.na(pf_status) && pf_status == "fail") {
        "skipped_preflight"    # Preflight explicitly rejected (e.g. rank-deficient design)
    } else if (!is.null(base_results)) {
        "success"              # Pipeline completed
    } else if (!is.null(coldata)) {
        "error"                # Coldata produced but DESeq2 failed unexpectedly
    } else {
        "skipped_filtering"    # Filtering eliminated all samples (e.g. no paired donors)
    }

    # Build summary row
    row <- make_summary_row(
        biosample            = stratum,
        contrast_id          = contrast_id,
        status               = status,
        n_control_samps      = pf_n_control,
        n_experimental_samps = pf_n_experimental,
        n_samples_total      = n_total,
        base_formula         = base_formula,
        n_degs_base          = n_degs_base,
        ruv_formula          = ruv_formula,
        n_degs_ruv           = n_degs_ruv,
        n_pathways_base      = n_pathways_base,
        n_pathways_ruv       = n_pathways_ruv
    )

    # Append preflight columns
    row$preflight_status     <- pf_status
    row$n_genes_preflight    <- pf_n_genes
    row$dropped_covariates   <- pf_dropped
    row$rank_ok              <- pf_rank_ok
    row$near_singular_vars   <- pf_near_sing
    row$preflight_warnings   <- pf_warnings

    row
}


#' Collect error messages from a stratum's error log file (if it exists)
get_error_log_path <- function(log_dir, stratum) {
    stratum_safe <- sanitize_names(stratum)
    path <- file.path(log_dir, paste0(stratum_safe, ".errors.log"))
    if (file.exists(path) && file.size(path) > 0) path else NULL
}


#' Collect and concatenate a result file across all strata, adding metadata columns
collect_results <- function(strata, results_dir, finals_file, inter_file) {
    rows <- lapply(strata, function(s) {
        finals_dir <- file.path(results_dir, s, "finals")
        inter_dir  <- file.path(results_dir, s, "intermediates")
        df <- safe_read(file.path(finals_dir, finals_file))
        if (is.null(df)) return(NULL)
        mi <- safe_read(file.path(inter_dir, inter_file))
        df$stratum <- s
        df$formula <- if (!is.null(mi) && "formula" %in% colnames(mi)) mi$formula[1] else NA
        df
    })
    do.call(rbind, Filter(Negate(is.null), rows))
}


#' Collect and concatenate fGSEA results across all strata
collect_fgsea <- function(strata, results_dir, file_name) {
    rows <- lapply(strata, function(s) {
        df <- safe_read(file.path(results_dir, s, "finals", file_name))
        if (is.null(df)) return(NULL)
        df$stratum <- s
        df
    })
    do.call(rbind, Filter(Negate(is.null), rows))
}


# =============================================================================
# Main
# =============================================================================

main <- function(config_path, contrast_id) {

    loaded       <- load_config(config_path)
    cfg          <- loaded$cfg
    contrasts_df <- loaded$contrasts_df

    cat(sprintf(
        "\n[%s] [INFO ] contrast=%s | step=07_aggregate_results | starting\n",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"), contrast_id
    ))

    results_dir  <- file.path(cfg$outputs$results_dir, contrast_id)
    log_dir      <- file.path(cfg$outputs$logs_dir, contrast_id)
    summary_dir  <- file.path(results_dir, "summary")
    dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

    # --- Discover all strata that were run ---
    all_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = FALSE)
    strata   <- all_dirs[all_dirs != "summary"]

    # --- Discover expected strata from contrasts.csv ---
    contrast_row <- contrasts_df[contrasts_df$contrast_id == contrast_id, ]
    if (nrow(contrast_row) > 0) {
        expected_strata_safe <- sanitize_names(parse_pipe_field(contrast_row$strata[1]))
    } else {
        expected_strata_safe <- character(0)
    }
    missing_strata <- setdiff(expected_strata_safe, strata)

    if (length(strata) == 0 && length(missing_strata) == 0) {
        cat(sprintf("[WARN] No stratum directories found under: %s\n", results_dir))
        return(invisible(NULL))
    }

    cat(sprintf("[INFO] Found %d strata with results: %s\n",
                length(strata), paste(strata, collapse = ", ")))
    if (length(missing_strata) > 0) {
        cat(sprintf("[INFO] %d expected strata not run: %s\n",
                    length(missing_strata), paste(missing_strata, collapse = ", ")))
    }

    # --- Collect per-stratum summary rows ---
    fgsea_fdr <- cfg$fgsea$fdr_threshold %||% 0.10

    summary_rows <- lapply(strata, function(stratum) {
        inter_dir  <- file.path(results_dir, stratum, "intermediates")
        finals_dir <- file.path(results_dir, stratum, "finals")

        collect_stratum_summary(
            inter_dir      = inter_dir,
            finals_dir     = finals_dir,
            contrast_id    = contrast_id,
            stratum        = stratum,
            fdr_threshold  = cfg$deseq2$fdr_threshold,
            l2fc_threshold = cfg$deseq2$l2fc_threshold,
            fgsea_fdr      = fgsea_fdr
        )
    })

    if (length(summary_rows) > 0) {
        contrast_summary <- do.call(rbind, summary_rows)
    } else {
        contrast_summary <- NULL
    }

    # --- Add rows for strata that were expected but never ran ---
    if (length(missing_strata) > 0) {
        missing_rows <- lapply(missing_strata, function(s) {
            row <- make_summary_row(
                biosample   = s,
                contrast_id = contrast_id,
                status      = "not_run"
            )
            row$preflight_status   <- NA
            row$n_genes_preflight  <- NA
            row$dropped_covariates <- NA
            row$rank_ok            <- NA
            row$near_singular_vars <- NA
            row$preflight_warnings <- NA
            row
        })
        missing_df <- do.call(rbind, missing_rows)
        contrast_summary <- if (!is.null(contrast_summary)) {
            rbind(contrast_summary, missing_df)
        } else {
            missing_df
        }
    }

    # =====================================================================
    # Write summary outputs
    # =====================================================================

    # --- 1. contrast_summary.csv ---
    summary_path <- file.path(summary_dir, "contrast_summary.csv")
    write_result(contrast_summary, summary_path, sep = ",")

    # --- 2. preflight_summary.csv ---
    preflight_rows <- lapply(strata, function(s) {
        read_preflight_report(file.path(results_dir, s, "intermediates"))
    })
    preflight_all <- do.call(rbind, Filter(Negate(is.null), preflight_rows))
    if (!is.null(preflight_all) && nrow(preflight_all) > 0) {
        write_result(preflight_all, file.path(summary_dir, "preflight_summary.csv"), sep = ",")
    }

    # --- 3. Consolidated DESeq2 tables ---
    deseq_base <- collect_results(strata, results_dir, "results_base.tsv", "model_info_base.csv")
    if (!is.null(deseq_base) && nrow(deseq_base) > 0) {
        write_result(deseq_base, file.path(summary_dir, paste0(contrast_id, "_deseq_base.tsv")))
    }

    deseq_ruv <- collect_results(strata, results_dir, "results_ruv.tsv", "model_info_ruv.csv")
    if (!is.null(deseq_ruv) && nrow(deseq_ruv) > 0) {
        write_result(deseq_ruv, file.path(summary_dir, paste0(contrast_id, "_deseq_ruv.tsv")))
    }

    # --- 4. Consolidated fGSEA tables ---
    fgsea_base <- collect_fgsea(strata, results_dir, "fgsea_base.tsv")
    if (!is.null(fgsea_base) && nrow(fgsea_base) > 0) {
        write_result(fgsea_base, file.path(summary_dir, paste0(contrast_id, "_fgsea_base.tsv")))
    }

    fgsea_ruv <- collect_fgsea(strata, results_dir, "fgsea_ruv.tsv")
    if (!is.null(fgsea_ruv) && nrow(fgsea_ruv) > 0) {
        write_result(fgsea_ruv, file.path(summary_dir, paste0(contrast_id, "_fgsea_ruv.tsv")))
    }

    # =====================================================================
    # Console summary
    # =====================================================================

    cat("\n-- Run Summary -------------------------------------------------------\n")
    cat(sprintf("  Contrast: %s\n", contrast_id))
    cat(sprintf("  Strata:   %d total\n\n", nrow(contrast_summary)))

    status_counts <- table(contrast_summary$status)
    for (s in names(status_counts)) {
        cat(sprintf("    %-15s %d\n", s, status_counts[[s]]))
    }

    if (length(missing_strata) > 0) {
        cat(sprintf("\n  Not run:  %d strata (%s)\n",
                    length(missing_strata), paste(missing_strata, collapse = ", ")))
    }

    # --- Preflight summary ---
    if ("preflight_status" %in% colnames(contrast_summary)) {
        cat("\n  Preflight:\n")
        pf_vals <- contrast_summary$preflight_status
        pf_vals[is.na(pf_vals)] <- "not available"
        pf_counts <- table(pf_vals)
        for (s in names(pf_counts)) {
            cat(sprintf("    %-15s %d\n", s, pf_counts[[s]]))
        }
    }

    cat("\n  Per-stratum DEG counts:\n")
    for (i in seq_len(nrow(contrast_summary))) {
        row <- contrast_summary[i, ]
        cat(sprintf("    %-20s  status=%-15s  DEGs_base=%-6s  DEGs_ruv=%-6s  pathways_base=%-6s  pathways_ruv=%s\n",
                    row$biosample,
                    row$status,
                    ifelse(is.na(row$n_degs_base), "--", row$n_degs_base),
                    ifelse(is.na(row$n_degs_ruv),  "--", row$n_degs_ruv),
                    ifelse(is.na(row$n_pathways_base), "--", row$n_pathways_base),
                    ifelse(is.na(row$n_pathways_ruv),  "--", row$n_pathways_ruv)))
    }

    # =====================================================================
    # Log aggregation
    # =====================================================================

    cat("\n-- Log Summary -------------------------------------------------------\n")
    log_summary_path <- file.path(summary_dir, "log_summary.tsv")

    strata_with_errors <- Filter(Negate(is.null),
        lapply(strata, get_error_log_path, log_dir = log_dir))

    if (length(strata_with_errors) > 0) {
        cat(sprintf("  %d stratum/strata had warnings or errors.\n", length(strata_with_errors)))

        summarize_logs(log_dir, log_summary_path)
        cat(sprintf("  Log summary: %s\n", log_summary_path))

        # Split into errors-only and warnings-only
        error_log_files <- list.files(log_dir, pattern = "\\.errors\\.log$",
                                      full.names = TRUE, recursive = FALSE)
        if (length(error_log_files) > 0) {
            all_lines   <- unlist(lapply(error_log_files, readLines, warn = FALSE))
            error_lines <- all_lines[grepl("\\[(ERROR|SKIP )\\]", all_lines)]
            warn_lines  <- all_lines[grepl("\\[WARN \\]", all_lines)]

            if (length(error_lines) > 0) {
                writeLines(error_lines, file.path(summary_dir, "errors.log"))
                cat(sprintf("  Errors:   %s\n", file.path(summary_dir, "errors.log")))
            }
            if (length(warn_lines) > 0) {
                writeLines(warn_lines, file.path(summary_dir, "warnings.log"))
                cat(sprintf("  Warnings: %s\n", file.path(summary_dir, "warnings.log")))
            }
        }
    } else {
        cat("  No warnings or errors found across all strata.\n")
    }

    cat("\n-- Complete ----------------------------------------------------------\n")
    cat(sprintf("  Results: %s\n\n", summary_dir))

    invisible(contrast_summary)
}


# =============================================================================
# CLI entry point
# =============================================================================

if (!interactive() && identical(sys.nframe(), 0L)) {

    opt_list <- list(
        make_option("--config",   type = "character", default = "config/config.yaml"),
        make_option("--contrast", type = "character",
                    help = "contrast_id to aggregate (e.g. T1D_vs_ND)")
    )
    opts <- parse_args(OptionParser(option_list = opt_list))
    if (is.null(opts$contrast)) stop("--contrast is required")

    main(opts$config, opts$contrast)
}
