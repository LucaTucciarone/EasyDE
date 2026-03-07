# =============================================================================
# stats_utils.R
# Small statistical helper functions
#
# These are pure utility functions with no side effects.
# =============================================================================


#' Convert a p-value to an asterisk significance label
#'
#' Used for annotating correlation plots with significance stars.
#'
#' @param x Numeric p-value (or vector of p-values)
#' @return Character string: "", "*", "**", or "***"
pvalue_to_stars <- function(x) {
    dplyr::case_when(
        x >= 0.05             ~ "",
        x < 0.05 & x >= 0.01 ~ "*",
        x < 0.01 & x >= 0.001 ~ "**",
        x < 0.001             ~ "***"
    )
}


#' Correlate RUVseq latent variables against known metadata variables
#'
#' Runs linear regression (latent_var ~ covariate) for each combination
#' of latent variable and covariate. Returns t-statistics and raw p-values
#' as matrices (covariates x latent vars).
#'
#' @param meta         data.frame of sample metadata (must contain latent var columns W_1, W_2, ...)
#' @param latent_names Character vector of latent variable column names (e.g. c("W_1", "W_2"))
#' @param correlate_vars Character vector of metadata column names to correlate against
#' @return Named list with two matrices: $t_mat and $p_mat (covariates x latent vars)
correlate_latent_vars <- function(meta, latent_names, correlate_vars) {

    n_lat  <- length(latent_names)
    n_vars <- length(correlate_vars)

    t_mat <- matrix(nrow = n_vars, ncol = n_lat,
                    dimnames = list(correlate_vars, latent_names))
    p_mat <- matrix(nrow = n_vars, ncol = n_lat,
                    dimnames = list(correlate_vars, latent_names))

    for (v in correlate_vars) {
        for (x in latent_names) {
            fit     <- lm(as.formula(paste0(x, " ~ ", v)), data = meta, na.action = na.omit)
            coefs   <- summary(fit)$coefficients
            # Guard: if the variable was dropped (e.g. only 1 level), record NA
            if (nrow(coefs) < 2) {
                t_mat[v, x] <- NA
                p_mat[v, x] <- NA
            } else {
                t_mat[v, x] <- coefs[2, 3]
                p_mat[v, x] <- coefs[2, 4]
            }
        }
    }

    return(list(t_mat = t_mat, p_mat = p_mat))
}


#' Find the best number of RUV latent factors (k) from ANOVA results
#'
#' Uses the "elbow" method: maximizes the second derivative of the
#' F-statistic curve across k values. This finds the point of
#' diminishing returns in variance reduction.
#'
#' @param anova_df data.frame with columns k_num (integer) and f_statistic (numeric)
#' @return Integer: the best k value
find_best_k <- function(anova_df) {

    df     <- anova_df[order(anova_df$k_num), ]
    f_stat <- df$f_statistic

    d1 <- c(NA, diff(f_stat))
    d2 <- c(diff(d1), NA)

    best_k <- df$k_num[which.max(d2)]

    return(best_k)
}


#' Count differentially expressed genes using FDR and LFC thresholds
#'
#' Centralized DEG counting used by 03, 05, and 07. Avoids duplicating
#' the same threshold logic in multiple scripts.
#'
#' @param result_df      data.frame with columns `padj` and `log2FoldChange`
#' @param fdr_threshold  Numeric: maximum adjusted p-value
#' @param l2fc_threshold Numeric: minimum absolute log2 fold change
#' @return Integer count of significant genes
count_degs <- function(result_df, fdr_threshold, l2fc_threshold) {
    sum(!is.na(result_df$padj) &
        result_df$padj < fdr_threshold &
        abs(result_df$log2FoldChange) > l2fc_threshold,
        na.rm = TRUE)
}


#' Sanitize column/variable names for use in R formulas
#'
#' Replaces non-alphanumeric characters with underscores.
#' We use underscores, NOT dots — dots have implicit method-dispatch meaning
#' in R (e.g. print.data.frame) and can cause subtle bugs in model formulas.
#'
#' @param x Character vector of names to sanitize
#' @return Character vector with non-alphanumeric characters replaced by "_"
sanitize_names <- function(x) {
    gsub("[^[:alnum:]]", "_", x)
}
