"C:/Users/roset/Desktop/Skin/GeneCounts.tsv.gz"

file.exists("C:/Users/roset/OneDrive/Desktop/skin/GSE249225_gene_counts.tsv")

gene_counts <- read.delim("C:/Users/roset/OneDrive/Desktop/skin/GSE249225_gene_counts.tsv", header=TRUE)

head(gene_counts)

dim(gene_counts)

colnames(gene_counts)

rownames(gene_counts) <- gene_counts$gene_id

gene_counts <- gene_counts[ , -1]

length(colnames((gene_counts)))

sample_info <- data.frame(

    sample = colnames(gene_counts),
  group  = c(rep("Young", 10), rep("Old", 15))
)

sample_info

head(sample_info)

str(sample_info)

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = gene_counts,
                              colData = sample_info,
                              design = ~ group)

dds <- DESeq(dds)

res <- results(dds)

res <- res[order(res$pvalue), ]

head(res)

sig_genes <- res[which(res$padj < 0.05), ]

head(sig_genes)

library(ggplot2)

res_clean <- res[!is.na(res$padj), ]

res_clean$threshold <- as.factor(res_clean$padj < 0.05 & abs(res_clean$log2FoldChange) > 1)

ggplot(res_clean, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
  geom_point(alpha=0.5, size=1.5) +
  scale_color_manual(values=c("grey","red")) +
  theme_minimal() +
  xlab("log2 Fold Change") 

top_genes <- head(res[order(-abs(res$log2FoldChange)), ], 10)

top_genes

rownames(gene_counts) <- gene_counts$gene_id
head(gene_counts$GeneID, 20)
top_genes_ids <- gene_counts$GeneID[1:10]
heat_data <- gene_counts[top_genes_ids, -1, drop=FALSE]
heat_data <- apply(heat_data, 2, function(x) as.numeric(as.character(x)))
colnames(gene_counts)
top_genes <- gene_counts[order(rowSums(gene_counts[,-1]), decreasing=TRUE), ][1:10, ]
gene_names <- top_genes$GeneID
counts <- rowSums(top_genes[,-1])
barplot(counts,
        names.arg = gene_names,
        las=2,
        col = "steelblue",
        main = "Top 10 Genes by Total Counts",
        ylab = "Sum of Counts")
library(ggplot2)
sample_cols <- colnames(gene_counts)[-1]

if(nrow(gene_counts) == 0 || length(sample_cols) == 0){
  stop("gene_counts empty or no sample columns")
}

if(nrow(top_genes) == 0){
  stop("no valid genes to plot")
}

top_genes <- gene_counts[order(rowSums(gene_counts[, sample_cols]), decreasing=TRUE), ][1:min(10, nrow(gene_counts)), ]

gene_names <- top_genes$GeneID
counts <- rowSums(top_genes[, sample_cols])
counts <- rowSums(top_genes[, sample_cols, drop=FALSE])
gene_names <- top_genes$GeneID

if(length(gene_names) != length(counts) || length(gene_names) == 0){
  gene_names <- "dummy_gene"
  counts <- 1e-6
}

df_plot <- data.frame(Gene = gene_names, Counts = counts)
library(ggplot2)


ggplot(df_plot, aes(x = Gene, y = Counts, fill = Counts)) +
  geom_bar(stat="identity") +
  scale_fill_gradient(low="white", high="firebrick3") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Top Genes Counts", y="Sum of Counts")

if(!require(tidyr)) install.packages("tidyr")
library(tidyr)
colnames(top_genes)
valid_cols <- intersect(c("GeneID", sample_cols), colnames(top_genes))


df_melt <- pivot_longer(top_genes[, valid_cols, drop=FALSE],
                        cols = all_of(sample_cols[sample_cols %in% colnames(top_genes)]),
                        names_to = "Sample",
                        values_to = "Counts")

df_melt$Counts <- as.numeric(as.character(df_melt$Counts))
df_melt$Counts[is.na(df_melt$Counts)] <- 0
gene_col <- colnames(top_genes)[1]


sample_cols <- colnames(top_genes)[-1]
df_melted <- pivot_longer(top_genes, cols = all_of(sample_cols), names_to = "sample", values_to = "counts")
head(df_melted)

colnames(df_melted)



ggplot(df_melted, aes(x = sample, y = "X1-002", fill = counts)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_gradient(low = "white", high = "red")

library(dplyr)

gene_ordered <- df_melted %>%
  group_by(X1_002) %>%
  summarize(total_counts = sum(counts)) %>%
  arrange(desc(total_counts))


df_long <- df_melted %>%
  pivot_longer(cols = counts, names_to = "variable", values_to = "counts")

ggplot(df_melted, aes(x = sample, y = X1_002, fill = counts)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_gradient(low = "white", high = "red")


library(dplyr)
library(ggplot2)
top_genes <- df_melted %>%
  group_by(X1_002) %>%
  summarize(total_counts = sum(counts)) %>%
  arrange(desc(total_counts)) %>%
  slice_head(n = 20) %>%
  pull(X1_002)

df_top <- df_melted %>%
  filter(X1_002 %in% top_genes)

ggplot(df_top, aes(x = sample, y = X1_002, fill = counts)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Sample", y = "Gene ID", fill = "Counts")



df_top$X1_002 <- factor(df_top$X1_002, levels = rev(unique(df_top$X1_002)))


ggplot(df_top, aes(x = sample, y = X1_002, fill = counts)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_gradient(low = "white", high = "red") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8)
  ) +
  labs(title = "Top 20 Genes Heatmap",
       x = "Sample",
       y = "Gene",
       fill = "Counts")



ggplot(df_top, aes(x = sample, y = X1_002, fill = log1p(counts))) +
  geom_tile() +
  theme_minimal() +
  scale_fill_gradient(low = "white", high = "red") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y = element_text(size = 8)
  ) +
  labs(title = "Top Genes Heatmap",
       x = "Sample",
       y = "Gene",
       fill = "log(Counts)")




