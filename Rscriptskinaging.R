library(DESeq2)
library(ggplot2)
library(pheatmap)

gene_counts <- read.delim("C:/Users/roset/OneDrive/Desktop/skin/GSE249225_gene_counts.tsv", header=TRUE)

rownames(gene_counts) <- gene_counts$gene_id
gene_counts <- gene_counts[ , -1]

sample_info <- data.frame(
  row.names = colnames(gene_counts),
  group = c(rep("Young",10), rep("Old",15))
)

dds <- DESeqDataSetFromMatrix(
  countData = gene_counts,
  colData = sample_info,
  design = ~ group
)


dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj), ]

res_clean <- res[!is.na(res$padj), ]

res_clean$threshold <- as.factor(
  res_clean$padj < 0.05 & abs(res_clean$log2FoldChange) > 1
)

ggplot(res_clean, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
  geom_point(alpha=0.6, size=1.5) +
  scale_color_manual(values=c("grey","red")) +
  theme_minimal() +
  labs(title="Volcano Plot",
       x="log2 Fold Change",
       y="-log10 Adjusted P-value")


vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup="group")

mat <- assay(vsd)

top_genes <- head(rownames(res_clean), 20)
mat_top <- mat[top_genes, ]

pheatmap(mat_top,
         scale="row",
         annotation_col = sample_info,
         main="Top 20 Differentially Expressed Genes")

plotMA(res, ylim=c(-5,5))





up_genes <- rownames(res_clean[
  res_clean$padj < 0.1 & abs(res_clean$log2FoldChange) > 0.5, ])
length(up_genes)




library(clusterProfiler)
library(org.Hs.eg.db)


gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!is.na(gene_list)]

gsea_res <- gseGO(
  geneList     = gene_list,
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENSEMBL",
  ont          = "BP",
  minGSSize    = 10,
  pvalueCutoff = 0.05
)
head(as.data.frame(gsea_res))
write.csv(as.data.frame(gsea_res), "GSEA_results.csv")

library(ggplot2)
library(clusterProfiler)

p <- dotplot(gsea_res, showCategory = 20) 
p <- p + ggtitle("Top Enriched GO Biological Processes")
p <- p + theme_minimal(base_size = 12)
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

print(p)


top10 <- head(as.data.frame(gsea_res)[order(as.data.frame(gsea_res)$p.adjust), ], 10)
top10_IDs <- top10$ID
barplot_data <- top10
barplot_data$Description <- factor(barplot_data$Description, levels = rev(barplot_data$Description))

ggplot(barplot_data, aes(x = NES, y = Description, fill = NES)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal(base_size = 12) +
  labs(title = "Top 10 Enriched Pathways (GSEA)", x = "Normalized Enrichment Score (NES)", y = "Pathway") +
  theme(axis.text.y = element_text(size = 10))

cnetplot(gsea_res, showCategory = 10, foldChange = gene_list)