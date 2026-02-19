# Skin Aging RNA-Seq Analysis

This repository contains the analysis of RNA-Seq data from epidermal tissues of visually younger and older women in the Youngster_odlie (Y_O) study. The analysis focuses on gene expression profiling and visualization of top genes.

## Files

- `GeneCounts.R` : R script with the complete analysis pipeline.
- `Figures/` : Folder containing all the plots generated from the analysis.

## Figures

### Bar Plots
![Top 10 Genes by Total Counts]
 figures/Barplot_Top10Genes.png
![Top Genes Counts] 
figures/Barplot_TopGenesCount.png
### Heatmaps
![X1_002 Heatmap] 
figures/Heatmap_X1-002.png
![X1_002 Heatmap with Counts Scale] 
figures/Heatmap_X1_002_scale.png
![Gene ID Heatmap]
figures/heatmap_gene_id.png
![Top 20 Genes Heatmap]  
figures/Heatmap_Top20Genes.png
![Top Genes Heatmap Log Counts] 
figures/top_genes_heatmap_log.png
### Volcano Plot
![Volcano Plot]
figures/Volcanoplot.png
## Description

This analysis includes:

- Identification of top genes based on total counts and log2 fold change.
- Generation of bar plots for top genes.
- Heatmaps showing gene expression across samples.
- Volcano plot highlighting significantly differentially expressed genes.
