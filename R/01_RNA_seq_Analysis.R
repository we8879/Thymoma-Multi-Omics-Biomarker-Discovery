# ===========================================================================
# Script: 01_RNA_seq_Analysis.R
# Project: Thymoma Multi-Omics Biomarker Discovery
# Description: Transcriptomic DEG analysis of thymoma using GSE177522 + GSE79978
#              - Merge raw counts and correct batch effect with ComBat-seq
#              - limma-voom differential expression (Tumor vs Normal)
#              - PCA, volcano plot, DEG heatmap and sample-correlation heatmap
# Date: 2026-09-05
# ===========================================================================
# Transcriptome Differential Expression Analysis (GSE177522 + GSE79978)
library(limma)
library(edgeR)
library(sva)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(pheatmap)
library(patchwork)

# ------------------------------------------------------------------------------
# 1. Data Loading and Preprocessing
# ------------------------------------------------------------------------------

setwd("D:/Thymoma_MultiOmics/data/rna_seq/")

data22 <- read.csv("22counts_self.csv", row.names = 1)
data78 <- read.csv("78counts.csv", row.names = 1)

# Retain only common genes
common_genes <- intersect(rownames(data22), rownames(data78))
combined_tpm <- cbind(data22[common_genes, ], data78[common_genes, ])

# Define batch and group variables
batch <- factor(c(rep("GSE177522", ncol(data22)), rep("GSE79978", ncol(data78))))

normal_samples <- c("GSM2109588", "GSM2109589", "GSM2109590",
                    "SRR14794296", "SRR14794297", "SRR14794298",
                    "SRR14794299", "SRR14794301")
group <- factor(ifelse(colnames(combined_tpm) %in% normal_samples, "Normal", "Tumor"))

sample_info <- data.frame(
  sample = colnames(combined_tpm),
  group = group,
  batch = batch,
  row.names = colnames(combined_tpm)
)

# ------------------------------------------------------------------------------
# 2. Batch Correction and Differential Expression (limma-voom)
# ------------------------------------------------------------------------------

# Create DGEList and filter low-expression genes
dge <- DGEList(counts = combined_tpm, group = group)
keep <- rowSums(cpm(dge) > 1) >= 3
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# ComBat-seq batch correction
combat_corrected <- ComBat_seq(
  counts = dge$counts,
  batch = batch,
  group = group
)

# Voom transformation and linear model
dge_combat <- dge
dge_combat$counts <- combat_corrected

design <- model.matrix(~ group)
colnames(design) <- make.names(colnames(design))

v <- voom(dge_combat, design, plot = TRUE, normalize = "quantile")
fit <- lmFit(v, design)
contrast_matrix <- makeContrasts(Tumor_vs_Normal = groupTumor, levels = design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit)

# Extract results
all_results <- topTable(fit, coef = "Tumor_vs_Normal", number = Inf, adjust.method = "BH")
all_results$Gene <- rownames(all_results)

# Define DE thresholds
padj_threshold <- 0.05
logFC_threshold <- 1

sig_genes <- all_results %>%
  filter(abs(logFC) > logFC_threshold & adj.P.Val < padj_threshold)
up_genes <- sig_genes %>% filter(logFC > 0)
down_genes <- sig_genes %>% filter(logFC < 0)

cat("Up:", nrow(up_genes), "Down:", nrow(down_genes), "Total:", nrow(sig_genes), "\n")

# ------------------------------------------------------------------------------
# 3. PCA: Before vs After Batch Correction
# ------------------------------------------------------------------------------

logcpm_raw <- cpm(dge, log = TRUE, prior.count = 1)
pca_raw <- prcomp(t(logcpm_raw))
pca_df_raw <- data.frame(PC1 = pca_raw$x[, 1], PC2 = pca_raw$x[, 2],
                         Group = group, Batch = batch)

pca_corrected <- prcomp(t(v$E))
pca_df_corrected <- data.frame(PC1 = pca_corrected$x[, 1], PC2 = pca_corrected$x[, 2],
                               Group = group, Batch = batch)

p1 <- ggplot(pca_df_raw, aes(PC1, PC2, color = Group, shape = Batch)) +
  geom_point(size = 3) + ggtitle("Before ComBat-seq") + theme_minimal()

p2 <- ggplot(pca_df_corrected, aes(PC1, PC2, color = Group, shape = Batch)) +
  geom_point(size = 3) + ggtitle("After ComBat-seq") + theme_minimal()

p1 + p2

# ------------------------------------------------------------------------------
# 4. Volcano Plot
# ------------------------------------------------------------------------------

volcano_data <- all_results %>%
  mutate(
    diffexpressed = case_when(
      logFC > logFC_threshold & adj.P.Val < padj_threshold ~ "UP",
      logFC < -logFC_threshold & adj.P.Val < padj_threshold ~ "DOWN",
      TRUE ~ "NOT SIG"
    )
  )

top_up <- volcano_data %>% filter(diffexpressed == "UP") %>% arrange(adj.P.Val) %>% head(10)
top_down <- volcano_data %>% filter(diffexpressed == "DOWN") %>% arrange(adj.P.Val) %>% head(10)
genes_to_label <- bind_rows(top_up, top_down)

ggplot(volcano_data, aes(logFC, -log10(P.Value), color = diffexpressed)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("UP" = "#E64B35", "DOWN" = "#3182bd", "NOT SIG" = "grey60")) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", linewidth = 0.5) +
  geom_text_repel(data = genes_to_label, aes(label = Gene), size = 3, max.overlaps = 50,
                  segment.color = "grey50", segment.size = 0.2) +
  labs(x = expression(log[2] * "Fold Change"), y = expression(-log[10] * "P-value"),
       title = "Differential Expression Analysis (Tumor vs Normal)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# ------------------------------------------------------------------------------
# 5. Top 50 Genes Heatmap (ComplexHeatmap)
# ------------------------------------------------------------------------------

library(ComplexHeatmap)
library(circlize)

top50_genes <- all_results %>% arrange(desc(abs(logFC))) %>% head(50) %>% pull(Gene)
expr_matrix <- v$E[top50_genes, ]
expr_scaled <- t(scale(t(expr_matrix)))

# Order samples by group
group_order <- order(group)
expr_scaled <- expr_scaled[, group_order]
sample_info_ordered <- sample_info[group_order, ]

ha_col <- HeatmapAnnotation(
  Group = sample_info_ordered$group,
  Batch = sample_info_ordered$batch,
  col = list(
    Group = c("Normal" = "#53AA56", "Tumor" = "#F3CA39"),
    Batch = c("GSE177522" = "#97C3D7", "GSE79978" = "#BBA7CD")
  )
)

Heatmap(
  expr_scaled,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#7EB5EA", "white", "#F57A7A")),
  top_annotation = ha_col,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_title = "Top 50 DEGs by |logFC| (Group-ordered)"
)

# ------------------------------------------------------------------------------
# 6. Sample Correlation Heatmap
# ------------------------------------------------------------------------------

cor_matrix <- cor(v$E, method = "pearson")
cor_matrix_ordered <- cor_matrix[group_order, group_order]

annotation_col_ordered <- data.frame(
  Group = factor(group[group_order], levels = c("Normal", "Tumor")),
  Batch = factor(batch[group_order], levels = c("GSE177522", "GSE79978")),
  row.names = colnames(cor_matrix_ordered)
)

ann_colors <- list(
  Group = c("Normal" = "#27AE60", "Tumor" = "#F1C40F"),
  Batch = c("GSE177522" = "#5DADE2", "GSE79978" = "#BB8FCE")
)

pheatmap(
  cor_matrix_ordered,
  annotation_col = annotation_col_ordered,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Pearson Correlation Heatmap (Group-ordered)",
  fontsize_row = 8,
  fontsize_col = 8
)
