---
title: "EasyDE Walkthrough — 02 DESeq2 and fGSEA"
output: html_document
---

# Step 2: Running DESeq2 and Pathway Analysis (fGSEA)

This notebook takes the cleaned count matrix and metadata from Notebook 01,
runs a standard DESeq2 differential expression analysis, and then performs
Gene Set Enrichment Analysis (fGSEA) to find activated/suppressed pathways.

```{r setup, message=FALSE, warning=FALSE}
pipeline_root <- normalizePath("../../")
knitr::opts_knit$set(root.dir = pipeline_root)
setwd(pipeline_root)

suppressMessages({
    library(DESeq2)
    library(fgsea)
    library(ggplot2)
    library(dplyr)
})

source("workflow/scripts/utils/stats_utils.R")
source("workflow/scripts/utils/plot_utils.R")
```

## 1. Prepare DESeq2 Inputs

Let's assume our metadata and counts were prepared in Step 01. We'll use the
pipeline's intermediate files for `T1D_vs_ND` in `Beta` cells to save time.

```{r load-data}
# Load the cleaned coldata and count matrix
coldata <- read.csv("results/T1D_vs_ND/Beta/intermediates/coldata.csv", row.names = 1)
counts_df <- read.csv("results/T1D_vs_ND/Beta/intermediates/counts_filtered.csv", row.names = 1)
counts_mat <- as.matrix(counts_df)

# Show the metadata we're feeding to DESeq2
print(head(coldata[, c("disease", "age", "sex")]))
```

## 2. Run DESeq2

We build a `DESeqDataSet` object. We will control for `age` and `sex`,
and test for differences between `disease` groups.

```{r run-deseq, message=FALSE}
# Create the dataset
# Note: The variable of interest (disease) goes LAST in the design formula.
dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = coldata,
    design    = ~ age + sex + disease
)

# Run the DESeq2 pipeline (normalization, dispersion, model fitting)
dds <- DESeq(dds)

cat("\nDESeq2 complete!\n")
```

## 3. View Results & Volcano Plot

We extract the results, comparing T1D vs ND. We'll use a False Discovery Rate (FDR)
cutoff of 0.05.

```{r view-results}
# Extract results for T1D vs ND
res <- results(dds, contrast = c("disease", "T1D", "ND"), alpha = 0.05)
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Count significant Differential Expressed Genes (DEGs)
n_degs <- sum(!is.na(res_df$padj) & res_df$padj < 0.05)
cat(sprintf("Found %d DEGs (FDR < 0.05)\n\n", n_degs))

# Look at the top 5 most significant genes
top_genes <- res_df %>% arrange(padj) %>% head(5)
print(top_genes[, c("gene", "log2FoldChange", "pvalue", "padj")])
```

```{r plot-volcano, fig.width=8, fig.height=6}
# Draw a volcano plot
plot_volcano(
    df = res_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_col = "gene",
    title = "Volcano Plot: Beta Cells (T1D vs ND)",
    fdr_line = 0.05,
    l2fc_line = 0
)
```

## 4. Run fGSEA Pathway Analysis

Now we rank all genes by their statistical significance (Wald statistic) and see
which biological pathways are enriched at the top (upregulated) or bottom (downregulated).

```{r run-fgsea}
# Load the Reactome/KEGG pathways GMT file
pathways <- fgsea::gmtPathways("resources/c2.cp.v7.5.1.symbols.gmt")

# Rank genes by DESeq2 'stat' (Wald statistic)
# Removing NA stats first
res_clean <- res_df[!is.na(res_df$stat), ]
ranks <- setNames(res_clean$stat, res_clean$gene)
ranks <- sort(ranks, decreasing = TRUE)

# Run fGSEA
# Note: In real pipeline runs, we also exclude ribosomal/mitochondrial genes first!
set.seed(42) # For reproducibility
fgsea_res <- fgseaMultilevel(pathways = pathways, stats = ranks,
                             minSize = 10, maxSize = 500)

# Filter to significant pathways (FDR < 0.1)
sig_pathways <- fgsea_res[fgsea_res$padj < 0.10, ]
sig_pathways <- sig_pathways[order(sig_pathways$pval), ]

cat(sprintf("\nFound %d significant pathways (FDR < 0.10)\n\n", nrow(sig_pathways)))

# Look at the top 5 pathways
print(head(sig_pathways[, c("pathway", "pval", "padj", "NES")]))
```

```{r plot-fgsea, fig.width=10, fig.height=8}
# Plot the top pathways
plot_fgsea_barplot(
    fgsea_sig_df = as.data.frame(sig_pathways),
    top_n = 20,
    title = "Top Enriched Pathways (T1D vs ND)"
)
```
