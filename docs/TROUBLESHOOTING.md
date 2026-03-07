# Troubleshooting

## Common Errors

**`skipped_preflight` with "100% of genes have >= 1 zero"**
The count matrix is extremely sparse — every gene has at least one sample
with zero counts. DESeq2 cannot estimate size factors with any normalization
method (ratio, poscounts, iterate). This typically occurs in small or rare
cell types in pseudobulk data. The stratum is cleanly rejected at preflight.

**`lfcShrink failed`**
apeglm requires at least one sample per group with non-zero counts. This is
a warning — `log2FoldChange_shrunk` will be NA but `log2FoldChange` is valid.

**`no count matrix found for stratum 'X'`**
The count file at `data/counts/X.counts.csv` (after safe-name conversion)
does not exist. Check that `strata` values in contrasts.csv match filenames
in `data/counts/`. Spaces and special characters become underscores.

**`W pruning removed all latent factors`**
All W factors correlated with the contrast variable (p < 0.05). The final
model runs without W factors — equivalent to base DESeq2. Consider whether
sample size is sufficient for RUVseq.

**fGSEA returns 0 significant pathways**
Not necessarily an error. Check: (1) the ranking statistic produces a
reasonable spread, (2) gene symbols in the GMT match your gene names,
(3) `fgsea.min_gene_set_size` is not too high.

**Snakemake shows "Nothing to be done"**
All output files already exist. Delete `results/` (or the specific contrast
directory) to force a re-run.

---

## Reading Log Files

Steps 02–06 write structured logs to `logs/{contrast_id}/{stratum}.log`:

```
[2026-03-06 14:32:01] [INFO]  contrast=Diabetic_vs_ND | biosample=Beta | step=start | ...
[2026-03-06 14:32:03] [INFO]  contrast=Diabetic_vs_ND | biosample=Beta | step=filter_samples | n_samples=89
[2026-03-06 14:32:15] [WARN]  contrast=Diabetic_vs_ND | biosample=Beta | step=lfc_shrinkage | msg=...
```

To see all warnings and errors after a run:

```bash
grep -rh "WARN\|ERROR\|SKIP" logs/
```

Step 07 also collects all flagged lines into
`results/{contrast}/summary/log_summary.tsv` for easy review.

---

## Resource Files

All resources are included in the repository — no additional downloads needed.

### Gene exclusion lists

Used by step 06 to remove ribosomal and mitochondrial genes before fGSEA ranking:

```
resources/gsea_files/rpl_genes.csv     ribosomal proteins, large subunit
resources/gsea_files/rps_genes.csv     ribosomal proteins, small subunit
resources/gsea_files/mtr_genes.csv     mitochondrially encoded genes
```

### Pathway GMT files

Three MSigDB gene set collections, pre-merged into a single default file:

```
resources/gsea_files/c2.cp.reactome.v*.symbols.gmt.txt    Reactome
resources/gsea_files/c2.cp.kegg_medicus.v*.symbols.gmt.txt KEGG Medicus
resources/gsea_files/c2.cp.kegg_legacy.v*.symbols.gmt.txt  KEGG Legacy
resources/gsea_files/reactome_kegg.gmt.txt                 merged (default)
```

To use a different collection, set `fgsea.gene_sets_file` in your config.

**Updating pathways:** Download new GMT files from
[MSigDB](https://www.gsea-msigdb.org/gsea/downloads.jsp) and re-merge:

```bash
cat c2.cp.reactome.*.gmt c2.cp.kegg_medicus.*.gmt c2.cp.kegg_legacy.*.gmt \
    > resources/gsea_files/reactome_kegg.gmt.txt
```
