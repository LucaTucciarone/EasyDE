---
title: "EasyDE Walkthrough — 01 Preparation and Preflight"
output: html_document
---

# Step 1: Data Preparation and Preflight QC

This notebook walks through the very first steps of the pipeline for a single comparison:
merging the count matrix with the metadata, filtering to the right samples, and checking
for zero-inflation.

We strip away the pipeline loops to focus purely on the data!

```{r setup, message=FALSE, warning=FALSE}
# Replace with the path to your pipeline root
pipeline_root <- normalizePath("../../")
knitr::opts_knit$set(root.dir = pipeline_root)
setwd(pipeline_root)

suppressMessages({
    library(data.table)
    library(dplyr)
})

# Load the utility functions we need
source("workflow/scripts/utils/io_utils.R")
source("workflow/scripts/utils/filter_utils.R")
source("workflow/scripts/utils/validation_utils.R")
```

## 1. Load the data

Let's pretend we are running the `T1D_vs_ND` contrast on the `Beta` cell stratum.

```{r load-data}
# 1. Load the pseudobulk counts for Beta cells
counts_df <- fread("resources/counts/Beta.counts.csv") %>% as.data.frame()
gene_names <- counts_df[[1]]
counts_mat <- as.matrix(counts_df[, -1])
storage.mode(counts_mat) <- "integer"
rownames(counts_mat) <- gene_names

cat("Count matrix loaded:\n")
cat("  Genes:", nrow(counts_mat), "\n")
cat("  Samples:", ncol(counts_mat), "\n\n")

print(head(counts_mat[, 1:5]))  # Show the first 5 samples and 6 genes
```

```{r load-metadata}
# 2. Load the sample metadata
metadata <- fread("resources/metadata.csv") %>% as.data.frame()

cat("\nMetadata loaded:\n")
cat("  Total rows:", nrow(metadata), "\n\n")

# Show the relevant columns for a few rows
cols_to_show <- c("biosample_id", "disease", "author_cell_type", "age", "sex")
print(head(metadata[, cols_to_show]))
```

## 2. Filter Samples

The count matrix contains *all* samples that had Beta cells. The metadata contains *all*
samples, even non-Beta ones. Let's merge them and filter down to just what we want to compare:
Type 1 Diabetes (T1D) vs Non-Diabetic (ND).

```{r filter-samples}
# Keep only metadata rows that match the columns in our count matrix
meta_filtered <- metadata[metadata$biosample_id %in% colnames(counts_mat), ]

# Filter to the contrast groups
meta_filtered <- meta_filtered[meta_filtered$disease %in% c("ND", "T1D"), ]

# Make sure the counts matrix only has these samples, and they are in the exact same order
counts_mat_filtered <- counts_mat[, meta_filtered$biosample_id, drop = FALSE]

cat("After filtering:\n")
cat("  Samples:", nrow(meta_filtered), "\n")
cat("  ND samples:", sum(meta_filtered$disease == "ND"), "\n")
cat("  T1D samples:", sum(meta_filtered$disease == "T1D"), "\n\n")

# Look at our clean dataset:
cat("Counts (filtered):\n")
print(head(counts_mat_filtered[, 1:4]))

cat("\nMetadata (filtered):\n")
print(head(meta_filtered[, cols_to_show]))
```

## 3. Filter Genes and Preflight Checks

Now that we have matched counts and metadata, we drop genes that are barely expressed,
and run the **Zero-Inflation Check** to ensure the data is suitable for DESeq2.

```{r filter-genes}
# Drop genes with less than 10 total counts across all remaining samples
min_counts <- 10
genes_to_keep <- rownames(counts_mat_filtered)[rowSums(counts_mat_filtered) >= min_counts]
counts_mat_filtered <- counts_mat_filtered[genes_to_keep, , drop = FALSE]

cat("After gene filtering:\n")
cat("  Genes retained:", nrow(counts_mat_filtered), "\n\n")
```

```{r preflight}
# The pipeline's check_zero_inflation utility
zi_result <- check_zero_inflation(counts_mat_filtered)

cat("Zero-Inflation Check:\n")
cat("  Threshold for rejection:", zi_result$threshold, "%\n")
cat("  Actual % genes with >=1 zero:", round(zi_result$pct_genes_with_zero, 1), "%\n")
cat("  Status:", ifelse(zi_result$zero_inflated, "REJECTED", "PASSED"), "\n")
```

If the status is `PASSED`, this dataset is ready for DESeq2 (Notebook 02)!
