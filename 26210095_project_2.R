library( ggplot2)
library( dplyr)
library( DESeq2)
library( pheatmap)
setwd("~/be623_projet2")
counts <- read.csv("airway_counts.csv", row.names = 1)
samples <- read.csv("airway_samples.csv")
dim(counts)
str(samples)
samples$lib_size <- colSums(counts)
samples %>%
  +     group_by(group) %>%
  +   summarize(mean_lib_size = mean(lib_size), n_samples = n())
rownames(samples) <- samples$run
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                                +                               colData = samples, design = ~ group)
dds$group <- relevel(dds$group, ref = "untrt")
levels(dds$group)
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
nrow(dds)
vsd <- vst(dds, blind = TRUE)
p2 <- plotPCA(vsd, intgroup = "group")
ggsave("plots/PCA.png", plot = p2)
vsd <- vst(dds, blind = TRUE)
dir.create("plots", showWarnings = FALSE)
p_pca <- plotPCA(vsd, intgroup = "group")
p_pca
ggsave("plots/PCA.png", p_pca, width = 7, height = 5, dpi = 300)
numerator   <- "untrt"
denominator <- "trt" 
paste0(numerator, "_vs_", denominator)
dds <- DESeq(dds)
numerator   <- "untrt"
denominator <- "trt"
comparison  <- paste0(numerator, "_vs_", denominator)
res <- results(dds, contrast = c("group", numerator, denominator))
head(res)
summary(res)
res <- res[order(res$padj), ]
es_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df$status <- "Not significant"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange >=  1] <- "Up"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- "Down"
table(res_df$status)
degs <- subset(res_df, status c("Up", "Down"))
nrow(degs)
head(degs)
write.csv(res_df, "results/all_results.csv", row.names = FALSE)
write.csv(deg_df, "results/DEGs.csv", row.names = FALSE)
deg_check <- read.csv("results/DEGs.csv")
head(deg_check, 5)
p3 <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = status)) +
  geom_point(alpha = 0.6) +
  labs(title = "trt_vs_untrt", x = "log2 Fold Change", y = "-log10 Adjusted P-value")
ggsave("plots/volcano.png", plot = p3)


