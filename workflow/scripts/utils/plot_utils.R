# =============================================================================
# plot_utils.R
# All plotting functions for EasyDE
#
# Every function returns a ggplot object - callers decide where to send it
# (pdf, png, print to screen). No dev.open/dev.off here.
# =============================================================================

suppressMessages({
    library(ggplot2)
    library(ggrepel)
    library(ggcorrplot)
    library(viridis)
    library(reshape2)
    library(dplyr)
})

# Resolve utils directory relative to this file's location, robustly
# stats_utils.R is a sibling - sourced via .script_dir set by the calling script
source(file.path(.script_dir, "utils", "stats_utils.R"))


#' Volcano plot from DESeq2 results
#'
#' @param deseq_df   data.frame with columns: log2FoldChange, pvalue, padj, gene
#' @param fdr        FDR threshold for coloring significant genes (default 0.05)
#' @param title      Plot title (e.g. the formula used)
#' @param subtitle   Plot subtitle (e.g. "number DEGs: 42")
#' @param n_label    Number of top significant genes to label (default 20)
#' @return ggplot object
plot_volcano <- function(deseq_df, fdr = 0.05, title = "", subtitle = "", n_label = 20) {

    # Separate significant genes for labeling
    sig_df <- deseq_df[!is.na(deseq_df$padj) & deseq_df$padj < fdr, ]
    top_df <- sig_df[order(sig_df$pvalue), ]
    top_df <- head(top_df, n_label)

    plot_df         <- deseq_df[!is.na(deseq_df$padj), ]
    plot_df$signif  <- ifelse(plot_df$padj < fdr, "Signif.", "N.S.")

    p <- ggplot(plot_df, aes(log2FoldChange, -log10(pvalue))) +
        geom_point(aes(col = signif), size = 0.5) +
        scale_color_manual(values = c("N.S." = "gray", "Signif." = "firebrick")) +
        labs(col = "", title = title, subtitle = subtitle) +
        theme_bw() +
        theme(plot.title = element_text(size = 8))

    if (nrow(top_df) > 0) {
        p <- p + ggrepel::geom_text_repel(
            data = top_df,
            aes(x = log2FoldChange, y = -log10(pvalue), label = gene),
            size = 3
        )
    }

    return(p)
}


#' Correlation heatmap: RUV latent vars vs known covariates
#'
#' Displays t-statistics as cell values with significance stars based on the
#' two-tier W exclusion strategy:
#'   - Contrast variable row: raw p < 0.05 (most aggressive, protects biology)
#'   - Other covariates: BH-adjusted per-W p < 0.05 (controls FDR)
#' Excluded Ws are marked with a red "X" header to make the exclusion visible.
#'
#' @param t_mat          Numeric matrix of t-statistics (covariates x latent vars)
#' @param p_mat_display  Numeric matrix of p-values used for exclusion decisions:
#'                       raw p for the contrast_var row, BH per-W for other rows
#' @param disease_ws     Character vector of excluded W names (default NULL)
#' @param contrast_var   Sanitized contrast variable name (for caption, default NULL)
#' @param title          Plot title string
#' @param t_limits       Symmetric limits for color scale (default c(-15, 15))
#' @return ggplot object
plot_latent_correlation <- function(t_mat, p_mat_display,
                                     disease_ws   = NULL,
                                     contrast_var = NULL,
                                     title = "", t_limits = c(-15, 15)) {

    # Get significance asterisks from the display p-value matrix
    p_labs        <- as.data.frame(apply(p_mat_display, 2, pvalue_to_stars))
    p_labs$Var1   <- as.factor(rownames(p_labs))
    p_labs_long   <- reshape2::melt(p_labs, id.vars = "Var1",
                                     variable.name = "Var2", value.name = "lab")

    cor_plot <- ggcorrplot(t_mat,
                           hc.order = FALSE,
                           lab      = TRUE,
                           title    = title,
                           lab_size = 2.25) +
        scale_fill_gradient2(
            low      = "#3A86FF",
            mid      = "white",
            high     = "#FD2244",
            limits   = t_limits,
            midpoint = 0,
            oob      = scales::squish
        )

    # Build caption explaining the two-tier significance thresholds
    contrast_label <- if (!is.null(contrast_var)) contrast_var else "contrast var"
    caption_text <- paste0(
        "* p<0.05  ** p<0.01  *** p<0.001\n",
        contrast_label, " row: raw p-value | other covariates: BH-adjusted per-W"
    )
    if (!is.null(disease_ws) && length(disease_ws) > 0) {
        caption_text <- paste0(caption_text,
            "\n\u2717 = excluded from final model (associated with known variable)")
    }
    cor_plot <- cor_plot + labs(fill = "t-statistic", caption = caption_text)

    # Align asterisks to the ggcorrplot data coordinates
    p_labs_long$in_df <- ifelse(
        is.na(match(
            paste0(p_labs_long$Var1, p_labs_long$Var2),
            paste0(cor_plot[["data"]]$Var1, cor_plot[["data"]]$Var2)
        )), "No", "Yes"
    )
    p_labs_long <- dplyr::filter(p_labs_long, in_df == "Yes")

    cor_plot <- cor_plot +
        geom_text(
            data    = p_labs_long,
            aes(x = Var1, y = Var2, label = lab),
            nudge_y = 0.25,
            size    = 5,
            inherit.aes = FALSE
        ) +
        guides(fill = guide_legend(title = "t-statistic"))

    # Mark excluded Ws with a red "X" marker above the column
    if (!is.null(disease_ws) && length(disease_ws) > 0) {
        # ggcorrplot maps columns to y-axis positions
        all_ws    <- colnames(t_mat)
        excl_pos  <- which(all_ws %in% disease_ws)
        if (length(excl_pos) > 0) {
            excl_df <- data.frame(
                x     = 0.4,                   # left edge of the plot
                y     = excl_pos,              # y positions of excluded Ws
                label = "\u2717"               # ✗ symbol
            )
            cor_plot <- cor_plot +
                geom_text(
                    data = excl_df,
                    aes(x = x, y = y, label = label),
                    color = "red", size = 5, fontface = "bold",
                    inherit.aes = FALSE
                )
        }
    }

    return(cor_plot)
}


#' RLE elbow plot: F-statistic vs number of RUV factors (k)
#'
#' Visualizes the ANOVA-based k selection. The chosen best_k is annotated.
#'
#' @param anova_df  data.frame with columns: k (factor), k_num, f_statistic, residual_variance
#' @param best_k    Integer: the selected best k value
#' @param title     Plot title string
#' @return ggplot object
plot_rle_elbow <- function(anova_df, best_k, title = "") {

    p <- ggplot(anova_df, aes(x = k, y = f_statistic, group = 1)) +
        geom_line() +
        geom_point(
            aes(size = residual_variance, fill = residual_variance),
            shape = 21, color = "black"
        ) +
        scale_fill_viridis(name = "Residual variance", option = "D") +
        scale_size_continuous(name = "Residual variance") +
        theme_classic() +
        labs(
            x     = "Number of RUV factors (k)",
            y     = "F-statistic\n(ANOVA: RLE ~ sample)",
            title = paste0(title, "\nBest-k -> ", best_k)
        ) +
        theme(
            axis.title  = element_text(size = 14),
            axis.text   = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
        )

    return(p)
}


#' fGSEA barplot: top N significant pathways by NES
#'
#' @param fgsea_sig_df  data.frame of significant fGSEA results (already filtered by FDR)
#'                      Must have columns: pathway, NES, padj
#' @param top_n         Max pathways to show (default 20, split top/bottom by NES)
#' @param title         Plot title
#' @return ggplot object, or NULL if no significant pathways
plot_fgsea_barplot <- function(fgsea_sig_df, top_n = 20, title = "Top significant pathways") {

    if (is.null(fgsea_sig_df) || nrow(fgsea_sig_df) == 0) {
        return(NULL)
    }

    n_show   <- min(top_n, nrow(fgsea_sig_df))
    plot_df  <- fgsea_sig_df[order(fgsea_sig_df$NES, decreasing = TRUE), ]
    plot_df  <- head(plot_df, n_show)

    p <- ggplot(plot_df, aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill = padj)) +
        coord_flip() +
        scale_fill_viridis(option = "C", direction = -1, name = "FDR") +
        labs(x = "Pathway", y = "Normalized Enrichment Score", title = title) +
        theme_minimal() +
        theme(text = element_text(size = 6))

    return(p)
}


#' Pipeline status overview heatmap
#'
#' Tile plot showing the status of each stratum (cell type) across the contrast.
#' Success statuses are shown in green shades, skips in orange shades,
#' errors in red, and not_run in grey. Cell labels show DEG counts where
#' available (RUV if present, else base).
#'
#' @param summary_df     data.frame: contrast_summary with columns biosample,
#'                       status, n_degs_base, n_degs_ruv
#' @param contrast_id    Character: contrast name for the title
#' @return ggplot object
plot_status_heatmap <- function(summary_df, contrast_id = "") {

    if (is.null(summary_df) || nrow(summary_df) == 0) return(NULL)

    # Choose the best DEG count to display: RUV > base
    summary_df$n_degs_display <- ifelse(
        !is.na(summary_df$n_degs_ruv),
        summary_df$n_degs_ruv,
        ifelse(!is.na(summary_df$n_degs_base), summary_df$n_degs_base, NA)
    )
    summary_df$deg_label <- ifelse(
        is.na(summary_df$n_degs_display), "",
        as.character(summary_df$n_degs_display)
    )

    # Order: status categories for consistent legend ordering
    status_levels <- c("success", "success_base_only",
                       "skipped_no_samples", "skipped_min_group",
                       "skipped_no_pairs", "skipped_preflight", "skipped",
                       "error", "not_run")
    summary_df$status <- factor(summary_df$status,
                                 levels = intersect(status_levels, unique(summary_df$status)))

    # Sort biosamples alphabetically
    summary_df$biosample <- factor(summary_df$biosample,
                                    levels = sort(unique(summary_df$biosample)))

    # Color palette
    status_colors <- c(
        "success"            = "#2ca02c",   # green
        "success_base_only"  = "#98df8a",   # light green
        "skipped_no_samples" = "#ff7f0e",   # orange
        "skipped_min_group"  = "#ffbb78",   # light orange
        "skipped_no_pairs"   = "#f0b27a",   # peach
        "skipped_preflight"  = "#e59866",   # tan
        "skipped"            = "#fad7a0",   # pale orange
        "error"              = "#d62728",   # red
        "not_run"            = "grey80"     # grey
    )
    # Only include colors for statuses present in the data
    used_colors <- status_colors[names(status_colors) %in% levels(summary_df$status)]

    p <- ggplot(summary_df, aes(x = biosample, y = 1, fill = status)) +
        geom_tile(color = "white", linewidth = 0.8) +
        geom_text(aes(label = deg_label), size = 2.5, color = "black") +
        scale_fill_manual(values = used_colors, name = "Status", drop = FALSE) +
        labs(
            title    = paste0("Pipeline Status: ", contrast_id),
            subtitle = sprintf("%d strata | %d success | %d skipped | %d errors",
                               nrow(summary_df),
                               sum(grepl("^success", summary_df$status)),
                               sum(grepl("^skipped", summary_df$status)),
                               sum(summary_df$status == "error")),
            x = "Cell type",
            y = ""
        ) +
        theme_minimal() +
        theme(
            axis.text.x     = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
            axis.text.y      = element_blank(),
            axis.ticks.y     = element_blank(),
            panel.grid       = element_blank(),
            legend.position  = "bottom",
            legend.text      = element_text(size = 7),
            plot.title       = element_text(size = 12, face = "bold"),
            plot.subtitle    = element_text(size = 9, color = "grey40")
        ) +
        guides(fill = guide_legend(nrow = 2))

    return(p)
}