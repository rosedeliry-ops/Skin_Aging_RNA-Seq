# Skin RNA-seq Analysis – Young vs Old

## Overview
Differential gene expression analysis of **Young** vs **Old** skin samples using RNA-seq. Includes DESeq2 analysis and GSEA for pathway enrichment.

## Data
- `GSE249225_gene_counts.tsv` (25 samples: 10 Young, 15 Old)
- Gene IDs: ENSG

## Workflow
- DESeq2: Volcano plot & Heatmap (top 20 genes)
- GSEA: GO Biological Processes enrichment
- Plots: Dotplot, Barplot (top 10 pathways), Cnetplot

## Output
- `GSEA_results.csv`
- Volcano plot, Heatmap, Dotplot, Barplot, Cnetplot (PDF/PNG)

## Dependencies
- R packages: `DESeq2`, `clusterProfiler`, `org.Hs.eg.db`, `ggplot2`, `pheatmap`, `dplyr`
