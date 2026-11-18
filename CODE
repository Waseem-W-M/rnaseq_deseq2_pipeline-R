# ==============================================================================
# End-to-End RNA-Seq Differential Gene Expression Analysis using DESeq2
# Author: Waseem Moshtaha
# Description: Simulates a realistic RNA-Seq count dataset (Control vs Treated),
#              performs QC, runs DESeq2 statistical modeling, and generates
#              publication-ready PCA, Volcano, and Heatmap visualizations.
# ==============================================================================

# 1. Install and Load Required Packages ----------------------------------------
required_cran <- c("ggplot2", "pheatmap", "ggrepel")
for (pkg in required_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2", update = FALSE, ask = FALSE)
}

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)

# Create output directory for plots and results
if (!dir.exists("results")) {
  dir.create("results")
}

# 2. Generate Realistic Simulated RNA-Seq Count Matrix -------------------------
set.seed(42)
n_genes <- 1000
n_samples <- 6 # 3 Control vs 3 Treated

gene_ids <- sprintf("GENE_%04d", 1:n_genes)
sample_ids <- c("Ctrl_1", "Ctrl_2", "Ctrl_3", "Treat_1", "Treat_2", "Treat_3")

# Base mean expression drawn from negative binomial distribution
base_means <- round(rlnorm(n_genes, meanlog = 4.5, sdlog = 1.2))
count_matrix <- matrix(0, nrow = n_genes, ncol = n_samples)
rownames(count_matrix) <- gene_ids
colnames(count_matrix) <- sample_ids

for (i in 1:n_genes) {
  mu <- max(base_means[i], 10)
  # Introduce true differential expression in first 100 genes
  if (i <= 50) {
    # Upregulated in Treated (4-fold increase)
    ctrl_counts <- rnbinom(3, mu = mu, size = 10)
    treat_counts <- rnbinom(3, mu = mu * 4, size = 10)
  } else if (i <= 100) {
    # Downregulated in Treated (4-fold decrease)
    ctrl_counts <- rnbinom(3, mu = mu * 4, size = 10)
    treat_counts <- rnbinom(3, mu = mu, size = 10)
  } else {
    # Non-differentially expressed genes
    ctrl_counts <- rnbinom(3, mu = mu, size = 10)
    treat_counts <- rnbinom(3, mu = mu, size = 10)
  }
  count_matrix[i, ] <- c(ctrl_counts, treat_counts)
}

# Construct Sample Metadata (colData)
col_data <- data.frame(
  row.names = sample_ids,
  condition = factor(c(rep("Control", 3), rep("Treated", 3)), levels = c("Control", "Treated"))
)

message("[1/5] Count matrix and metadata initialized successfully.")

# 3. Run DESeq2 Differential Expression Pipeline -------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ condition
)

# Filter out low-count genes (at least 10 reads across all samples)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2 normalization, dispersion estimation, and Wald test
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Treated", "Control"), alpha = 0.05)
res_ordered <- res[order(res$padj), ]

# Export full results table to CSV
res_df <- as.data.frame(res_ordered)
res_df$gene_id <- rownames(res_df)
write.csv(res_df, file = "results/differential_expression_results.csv", row.names = FALSE)

message("[2/5] DESeq2 statistical analysis completed and saved to CSV.")

# 4. Quality Control: Variance Stabilizing Transformation & PCA Plot -----------
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5, show.legend = FALSE) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("Principal Component Analysis (VST Normalized Counts)") +
  theme_bw(base_size = 12) +
  scale_color_manual(values = c("Control" = "#2b5c8f", "Treated" = "#d95f02"))

ggsave("results/pca_plot.png", plot = pca_plot, width = 7, height = 5, dpi = 300)
message("[3/5] PCA plot saved to results/pca_plot.png")

# 5. Visualization: Volcano Plot -----------------------------------------------
res_df$significance <- "Not Significant"
res_df$significance[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >= 1] <- "Upregulated"
res_df$significance[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- "Downregulated"

top_labelled <- head(res_df[!is.na(res_df$padj), ], 10)

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(data = top_labelled, aes(label = gene_id), size = 3, max.overlaps = 15, color = "black") +
  scale_color_manual(values = c("Upregulated" = "#d73027", "Downregulated" = "#4575b4", "Not Significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(
    title = "Volcano Plot: Treated vs Control",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value (FDR)",
    color = "Regulation Status"
  ) +
  theme_minimal(base_size = 12)

ggsave("results/volcano_plot.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)
message("[4/5] Volcano plot saved to results/volcano_plot.png")

# 6. Visualization: Heatmap of Top 30 Differentially Expressed Genes -----------
top30_genes <- head(rownames(res_ordered), 30)
mat_top30 <- assay(vsd)[top30_genes, ]
# Z-score scaling across rows (genes)
mat_scaled <- t(scale(t(mat_top30)))

annotation_col <- data.frame(Condition = col_data$condition)
rownames(annotation_col) <- rownames(col_data)

png("results/top30_degs_heatmap.png", width = 2000, height = 2400, res = 300)
pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  color = colorRampPalette(c("#313695", "white", "#a50026"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  fontsize_row = 8,
  main = "Top 30 Differentially Expressed Genes (Z-Score Normalized)"
)
dev.off()

message("[5/5] Heatmap saved to results/top30_degs_heatmap.png")
message("Pipeline execution finished successfully!")
