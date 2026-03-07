---
title: "EasyDE Walkthrough — 03 Batch Correction with RUVSeq"
output: html_document
---

# Step 3: Removing Unwanted Variation (RUVseq)

Technical differences between samples (library prep batches, sequencing lanes) can
mask out the real biological signal we care about. In this notebook, we use the
results from Notebook 02 to discover these hidden technical factors, then we
add them to our DESeq2 formula to "correct" the data.

```{r setup, message=FALSE, warning=FALSE}
pipeline_root <- normalizePath("../../")
knitr::opts_knit$set(root.dir = pipeline_root)
setwd(pipeline_root)

suppressMessages({
    library(DESeq2)
    library(RUVSeq)
    library(EDASeq)
    library(ggplot2)
    library(dplyr)
})

source("workflow/scripts/utils/stats_utils.R")
```

## 1. Identify "Negative Controls"

RUVseq assumes that genes which are *not* biologically changing between groups
can act as anchors. If these "negative control" genes vary, it must be due to
technical noise.

We find them by taking the uncorrected DESeq2 results (from Step 02) and selecting
genes with very high p-values (p > 0.5), meaning they are definitely not
differentially expressed.

```{r find-controls}
# Load the uncorrected DESeq2 results from Step 02
base_res <- read.delim("results/T1D_vs_ND/Beta/finals/results_base.tsv")

# Find empirical control genes: high p-value (p > 0.5)
empirical_controls <- base_res$gene[!is.na(base_res$pvalue) & base_res$pvalue > 0.5]

cat(sprintf("Found %d empirical negative control genes\n", length(empirical_controls)))
```

## 2. Run RUVg

We run the `RUVg` function on our raw counts, using these control genes.
This function estimates hidden technical factors (we call them `W`).
We have to tell it how many factors to find (`k`). Let's find 2 factors (`k=2`).

```{r run-ruvg}
# Load the filtered counts & coldata again
counts_df <- read.csv("results/T1D_vs_ND/Beta/intermediates/counts_filtered.csv", row.names = 1)
counts_mat <- as.matrix(counts_df)
coldata <- read.csv("results/T1D_vs_ND/Beta/intermediates/coldata.csv", row.names = 1)

# RUVg requires raw integer counts
# The output contains the W factors
ruvg_res <- RUVg(counts_mat, empirical_controls, k = 2)

cat("\nComputed W factors (k=2):\n")
w_factors <- ruvg_res$W
print(head(w_factors))

# Add these W factors to our sample metadata
coldata$W_1 <- w_factors[, 1]
coldata$W_2 <- w_factors[, 2]
```

## 3. Protect Biological Signal

Sometimes, a technical factor (`W`) accidentally captures real biological disease
signal. If we correct for it, we wipe out the biology! So we must check if
either `W_1` or `W_2` strongly correlates with our `disease` variable.

```{r check-w}
# In the pipeline, we use `identify_disease_ws()` for this correlation check.
# Let's see if either W factor correlates with disease status.
# (For a two-group factor like T1D/ND, this is like running a t-test).
cat("Checking W factors against disease:\n")

t.test.w1 <- t.test(coldata$W_1 ~ coldata$disease)
cat(sprintf("  W_1 vs Disease: p = %.3f\n", t.test.w1$p.value))

t.test.w2 <- t.test(coldata$W_2 ~ coldata$disease)
cat(sprintf("  W_2 vs Disease: p = %.3f\n", t.test.w2$p.value))

if (t.test.w1$p.value < 0.05) cat("  WARNING: W_1 correlates with disease. It should be dropped.\n")
if (t.test.w2$p.value < 0.05) cat("  WARNING: W_2 correlates with disease. It should be dropped.\n")

# Normally, if a W factor correlates with disease, we EXCLUDE it from the formula below.
# Here we'll assume both are safe to use.
```

## 4. Re-run DESeq2 with W factors

Now we add the safe W factors to our DESeq2 design formula. By including them,
DESeq2 effectively "subtracts" the technical noise out before testing the
disease effect.

```{r run-deseq-corrected, message=FALSE}
# New formula: ~ W_1 + W_2 + age + sex + disease
dds_ruv <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = coldata,
    design    = ~ W_1 + W_2 + age + sex + disease
)

dds_ruv <- DESeq(dds_ruv)

# Extract new results
res_ruv <- results(dds_ruv, contrast = c("disease", "T1D", "ND"), alpha = 0.05)
res_ruv_df <- as.data.frame(res_ruv)
res_ruv_df$gene <- rownames(res_ruv_df)

n_degs_ruv <- sum(!is.na(res_ruv_df$padj) & res_ruv_df$padj < 0.05)
cat(sprintf("\nFound %d DEGs after RUVseq correction (FDR < 0.05)\n", n_degs_ruv))
```

By removing technical noise, the statistical power increases — usually resulting
in many more biologically meaningful DEGs being discovered!

We could now take `res_ruv_df` and run volcano plots and `fGSEA` exactly like we
did in Notebook 02.
