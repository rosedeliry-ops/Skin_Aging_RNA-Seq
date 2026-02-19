# Skin Aging RNA-Seq Analysis

This repository contains the analysis of RNA-Seq data from epidermal tissues of visually younger and older women in the Youngster_odlie (Y_O) study. The analysis focuses on gene expression profiling and visualization of top genes.

## Files

- `GeneCounts.R` : R script with the complete analysis pipeline.
- `Figures/` : Folder containing all the plots generated from the analysis.

## Figures

### Bar Plots
![Top 10 Genes by Total Counts](Figures/Barplot_Top10Genes.png)
![Top Genes Counts](Figures/Barplot_TopGenesCount.png)
### Heatmaps
![X1_002 Heatmap](Figures/heatmap_X1_002.png)
![X1_002 Heatmap with Counts Scale](Figures/heatmap_X1_002_scale.png)
![Gene ID Heatmap](Figures/heatmap_gene_id.png)
![Top 20 Genes Heatmap](Figures/Heatmap_Top20Genes.png)
![Top Genes Heatmap Log Counts](Figures/top_genes_heatmap_log.png)

### Volcano Plot
![Volcano Plot](Figures/Volcanoplot.png)

## Description

This analysis includes:

- Identification of top genes based on total counts and log2 fold change.
- Generation of bar plots for top genes.
- Heatmaps showing gene expression across samples.
- Volcano plot highlighting significantly differentially expressed genes.