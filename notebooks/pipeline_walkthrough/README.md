# EasyDE Walkthrough — Educational Notebooks

This folder contains a stripped-down, educational interactive version of the
EasyDE pipeline. 

Instead of dealing with Snakemake, complex logging, or bulk processing loops,
these notebooks walk through **a single contrast and stratum** (e.g., examining
Type 1 Diabetes vs Non-Diabetic in Beta cells). 

They are designed for students to **see the data at every step**, understand
the core concepts of RNA-Seq analysis, and experiment with the parameters.

## Notebook Sequence

| Notebook | Purpose |
|----------|---------|
| [01_Prep_and_Preflight.md](01_Prep_and_Preflight.md) | **Start here.** Shows how pseudo-bulk counts and metadata are joined. Demonstrates filtering and the zero-inflation preflight check. |
| [02_DESeq_and_fGSEA.md](02_DESeq_and_fGSEA.md) | Runs a standard DESeq2 differential expression test. Extracts results, draws a volcano plot, and runs fGSEA pathway analysis. |
| [03_DESeq_RUVSeq_fGSEA.md](03_DESeq_RUVSeq_fGSEA.md) | Explains **Batch Correction**. Shows how to find empirical control genes, use RUVSeq to identify technical noise factors, and add them to the DESeq2 formula to rescue biological signal. |
| [02_alt_Paired_Samples.md](02_alt_Paired_Samples.md) | **Alternative path.** Shows how to handle *paired* experimental designs using donor-blocking in the formula, explaining why RUVSeq is skipped for these datasets. |

## How to Use

1. Open these notebooks in RStudio or Jupyter.
2. Edit the `setwd()` path in the first code chunk to point to your EasyDE project root.
3. Run the chunks step-by-step. Don't just click "Run All" — look at the table outputs
   (`head(...)`) and read the explanations to understand *why* the pipeline makes
   certain decisions.
