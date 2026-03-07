---
title: "EasyDE Walkthrough — 02_alt Paired Samples"
output: html_document
---

# Step 2 (Alternative): Paired Sample DESeq2

This notebook replaces Notebook 02 & 03 for **paired experimental designs**
— where you have multiple samples from the exact same donor (e.g., PBMC cells
before and after treatment from Patients 1, 2, and 3).

In a paired design, the differences *between* individual donors are usually much larger
than the differences caused by the treatment. We must "block" on the donor to let DESeq2
subtract out each person's baseline. **Because we explicitly model the donor baseline,
we skip RUVSeq batch correction entirely.**

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

In this example, imagine we have `CD4_T` cells treated with `Stimulated` vs `Unstimulated`,
with pairs drawn from 5 distinct `donor_id` accessions.

```{r load-paired-data}
# (Simulating loading the 'CD4_T' stratum from the 'Treatment_Effect' contrast)
# coldata <- read.csv("results/Treatment_Effect/CD4_T/intermediates/coldata.csv", row.names = 1)
# counts_mat <- as.matrix(read.csv("results/Treatment_Effect/CD4_T/intermediates/counts_filtered.csv", row.names = 1))

# Let's peek at a hypothetical coldata layout for this:
hypothetical_coldata <- data.frame(
    sample_id = paste0("Sample_", 1:6),
    donor_id  = c("Donor_A", "Donor_B", "Donor_C", "Donor_A", "Donor_B", "Donor_C"),
    treatment = c("Unstimulated", "Unstimulated", "Unstimulated", "Stimulated", "Stimulated", "Stimulated")
)
cat("Notice how each donor appears exactly TWICE — once in each treatment arm:\n")
print(hypothetical_coldata)
```

## 2. Run Paired DESeq2

To perform paired analysis, we just add the donor column to the beginning of the
design formula!

`design = ~ donor_id + treatment`

```{r run-paired-deseq, eval=FALSE}
# Create the dataset (we place the blocking factor FIRST, and the condition LAST)
dds_paired <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = coldata,
    design    = ~ donor_id + treatment
)

# Run the DESeq2 pipeline
dds_paired <- DESeq(dds_paired)

cat("\nPaired DESeq2 complete!\n")
```

## 3. Why Not RUVSeq?

In Notebook 03, we used RUVSeq to find hidden factors of unwanted variation (`W_1`, `W_2`).
However, in a paired design:

1. The largest source of technical/unwanted variation is almost always **inter-patient genetics and baseline state.**
2. By including `donor_id` in the formula, we explicitly provide DESeq2 with the 
   exact boundaries of this variation.
3. If we also ran RUVSeq, the `W` factors would mostly just rediscover the `donor_id` groupings!

So, for paired pipelines:
- **Skip RUVSeq**.
- Proceed directly to fGSEA pathway analysis using the paired DESeq2 results.

```{r paired-results, eval=FALSE}
# Extract paired results
res_paired <- results(dds_paired, contrast = c("treatment", "Stimulated", "Unstimulated"), alpha = 0.05)
res_paired_df <- as.data.frame(res_paired)
res_paired_df$gene <- rownames(res_paired_df)

n_degs <- sum(!is.na(res_paired_df$padj) & res_paired_df$padj < 0.05)
cat(sprintf("\nFound %d DEGs controlling for donor baseline (FDR < 0.05)\n\n", n_degs))

# Run Volcano plot
plot_volcano(
    df = res_paired_df,
    lfc_col = "log2FoldChange",
    padj_col = "padj",
    gene_col = "gene",
    title = "Volcano Plot: CD4_T (Stimulated vs Unstimulated, Paired)",
    fdr_line = 0.05,
    l2fc_line = 0
)

# ... and run fGSEA just like in Notebook 02!
```
