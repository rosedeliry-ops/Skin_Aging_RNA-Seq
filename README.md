Project: Skin RNA-seq Analysis – Young vs Old
Overview
This project analyzes differential gene expression between Young and Old skin samples using RNA-seq data. The workflow includes DESeq2 analysis, visualization of significant genes, and pathway enrichment analysis using GSEA.
Data
Input file: GSE249225_gene_counts.tsv (gene counts matrix)
Samples: 25 total (10 Young, 15 Old)
Gene IDs: ENSG identifiers
Workflow
DESeq2 Analysis
Differential expression between Young and Old groups
Volcano plot to visualize up- and down-regulated genes
Heatmap of top 20 genes by adjusted p-value
GSEA (Gene Set Enrichment Analysis)
Ranked gene list based on log2 fold change
Enrichment of GO Biological Processes
Dotplot visualization for top pathways
Barplot of top 10 enriched pathways
Cnetplot for top pathways showing gene-to-pathway relationships
Output
GSEA_results.csv: Full table of enriched pathways with NES, p-value, and q-value
Volcano plot (PDF/PNG)
Heatmap of top genes (PDF/PNG)
GSEA dotplot (PDF/PNG)
Barplot of top pathways (PDF/PNG)
Cnetplot (PDF/PNG)
Dependencies
R packages: DESeq2, clusterProfiler, org.Hs.eg.db, ggplot2, pheatmap, dplyr
Notes
Thresholds for differential expression: padj < 0.05, |log2FoldChange| > 1 (adjustable)
GSEA is recommended over classical enrichment for small numbers of significant genes
All plots are generated in R using ggplot2 and clusterProfiler
اگر بخوای، می‌تونم همینو برات نسخه Markdown آماده برای GitHub یا PDF هم بسازم تا مستقیم آپلودش کنی.
میخوای بسازم؟
