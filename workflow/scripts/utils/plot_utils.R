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
#' @param t_mat          Numeric matrix of t-statistics (covariates x latent vars)
#' @param p_mat_adj      Numeric matrix of Bonferroni-adjusted p-values (same dims)
#' @param title          Plot title string
#' @param t_limits       Symmetric limits for color scale (default c(-15, 15))
#' @return ggplot object
plot_latent_correlation <- function(t_mat, p_mat_adj, title = "", t_limits = c(-15, 15)) {

    # Get significance asterisks
    p_labs        <- as.data.frame(apply(p_mat_adj, 2, pvalue_to_stars))
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
        ) +
        labs(
            fill    = "t-statistic",
            caption = "Bonferroni corrected p-values: * p<0.05, ** p<0.01, *** p<0.001"
        )

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