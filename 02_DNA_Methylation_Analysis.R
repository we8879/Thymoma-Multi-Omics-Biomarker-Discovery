# ===========================================================================
# Script: 02_DNA_Methylation_Analysis.R
# Project: Thymoma Multi-Omics Biomarker Discovery
# Description: DNA methylation (Illumina EPIC) analysis
#              - ChAMP normalization (run once, section 1)
#              - PCA-based detection and removal of outlier samples
#              - Differentially methylated position (DMP) analysis with ChAMP
# Date: 2026-09-05
# ===========================================================================
# Load packages
suppressPackageStartupMessages({
  library(minfi)
  library(ChAMP)
  library(limma)
  library(ChAMPdata)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
})

# Set working directory
setwd("D:/Thymoma_MultiOmics/data/methylation/")

#=============================================================================
# 1. Data normalization (run once)
#=============================================================================
myLoad$beta <- myLoad$beta[!rowSums(is.na(myLoad$beta)) > 0, ]

champ.obj <- champ.norm(beta = myLoad$beta, arraytype = "EPIC", cores = 5)
save(champ.obj, file = "./CHAMP_Normalization/normalized_data.RData")

#=============================================================================
# 2. Load normalized data
#=============================================================================
load("./CHAMP_Normalization/normalized_data.RData")
beta_matrix <- champ.obj

#=============================================================================
# 3. Define sample groups
#=============================================================================
group_vector <- c(rep("Tumor", 26), rep("Normal", 9), rep("Tumor", 87))
phenoData <- data.frame(
  row.names = colnames(beta_matrix),
  group = factor(group_vector, levels = c("Tumor", "Normal"))
)

#=============================================================================
# 4. PCA analysis - detect outliers
#=============================================================================
pca_input <- t(beta_matrix)

# Remove zero-variance probes
zero_var_probes <- which(apply(pca_input, 2, var) == 0)
if (length(zero_var_probes) > 0) {
  pca_input <- pca_input[, -zero_var_probes]
}

pca_result <- prcomp(pca_input, scale. = TRUE, center = TRUE)

pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = phenoData[rownames(pca_result$x), "group"],
  Sample = rownames(pca_result$x)
)

var_explained <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

# Identify tumor samples within Normal 95% confidence ellipse
normal_samples <- pca_data[pca_data$Group == "Normal", ]
tumor_samples <- pca_data[pca_data$Group == "Tumor", ]

normal_mean <- c(mean(normal_samples$PC1), mean(normal_samples$PC2))
normal_cov <- cov(cbind(normal_samples$PC1, normal_samples$PC2))
chi_sq_critical <- qchisq(0.95, df = 2)

tumor_mahalanobis <- apply(tumor_samples[, c("PC1", "PC2")], 1, function(x) {
  diff <- x - normal_mean
  sqrt(t(diff) %*% solve(normal_cov) %*% diff)
})

tumor_in_normal_ellipse <- tumor_samples[tumor_mahalanobis <= sqrt(chi_sq_critical), ]

#=============================================================================
# 5. PCA plot with outlier labeling
#=============================================================================
pca_data$SpecialLabel <- ifelse(
  pca_data$Sample %in% tumor_in_normal_ellipse$Sample,
  "Tumor_in_Normal_Ellipse",
  as.character(pca_data$Group)
)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.8) +
  geom_point(
    data = subset(pca_data, SpecialLabel == "Tumor_in_Normal_Ellipse"),
    aes(shape = SpecialLabel),
    size = 5, color = "purple", fill = "purple"
  ) +
  geom_text_repel(
    data = subset(pca_data, SpecialLabel == "Tumor_in_Normal_Ellipse"),
    aes(label = Sample),
    color = "purple", size = 3,
    box.padding = 0.8, max.overlaps = 50
  ) +
  scale_color_manual(
    values = c("Tumor" = "#F55A5A", "Normal" = "#65B2FA"),
    labels = c(paste0("Normal (n=", sum(phenoData$group == "Normal"), ")"), 
               paste0("Tumor (n=", sum(phenoData$group == "Tumor"), ")"))
  ) +
  scale_shape_manual(
    name = "Special samples",
    values = c("Tumor_in_Normal_Ellipse" = 23),
    labels = paste0("Tumor in Normal ellipse (n=", nrow(tumor_in_normal_ellipse), ")")
  ) +
  labs(
    title = "PCA of DNA Methylation (EPIC Array)",
    subtitle = paste0("Purple diamonds: tumor samples within Normal 95% CI ellipse (n=", 
                      nrow(tumor_in_normal_ellipse), ")"),
    x = paste0("PC1 (", var_explained[1], "% Variance)"),
    y = paste0("PC2 (", var_explained[2], "% Variance)")
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "purple"),
    panel.grid = element_blank()
  )

print(pca_plot)

#=============================================================================
# 6. Remove outlier samples
#=============================================================================
samples_to_remove_1 <- tumor_in_normal_ellipse$Sample

additional_samples_to_remove <- c(
  "GSM6751739_205555380022_R03C01", "GSM6751755_205549600090_R02C01",
  "GSM6751748_205549600078_R04C01", "GSM6751745_205549600070_R03C01",
  "GSM6751746_205549600070_R03C01", "GSM6751707_205549600070_R06C01",
  "GSM6751718_205549610016_R07C01", "GSM6751724_205549610017_R05C01",
  "GSM6751762_205549600094_R01C01", "GSM6751769_205549600122_R01C01",
  "GSM6751709_205549600078_R01C01", "GSM6751819_205549600149_R07C01",
  "GSM6751711_205549600078_R03C01", "GSM6751710_205549600078_R02C01",
  "GSM6751708_205549600070_R08C01", "GSM6751746_205549600070_R04C01",
  "GSM6751749_205549600078_R05C01", "GSM6751750_205549610176_R05C01",
  "GSM6751742_205549600070_R05C01", "GSM6751744_205549600070_R02C01",
  "GSM6751757_205549600090_R04C01", "GSM6751772_205549610176_R03C01",
  "GSM6751811_205549600149_R02C01", "GSM6751715_205549610016_R04C01",
  "GSM6751730_205549610175_R05C01", "GSM6751828_205549610015_R08C01",
  "GSM6751824_205549610015_R04C01", "GSM6751818_205549600149_R06C01"
)

all_samples_to_remove <- c(samples_to_remove_1, additional_samples_to_remove)

final_keep <- !rownames(pca_input) %in% all_samples_to_remove
pca_input_cleaned <- pca_input[final_keep, ]

zero_var_probes_cleaned <- which(apply(pca_input_cleaned, 2, var, na.rm = TRUE) == 0)
if (length(zero_var_probes_cleaned) > 0) {
  pca_input_cleaned <- pca_input_cleaned[, -zero_var_probes_cleaned]
}

pca_result_cleaned <- prcomp(pca_input_cleaned, scale. = TRUE, center = TRUE)

pca_data_cleaned <- data.frame(
  PC1 = pca_result_cleaned$x[, 1],
  PC2 = pca_result_cleaned$x[, 2],
  Group = phenoData[rownames(pca_result_cleaned$x), "group"],
  Sample = rownames(pca_result_cleaned$x)
)

var_explained_cleaned <- round(100 * summary(pca_result_cleaned)$importance[2, 1:2], 1)

#=============================================================================
# 7. Final PCA plot with sample labels
#=============================================================================
pca_plot_final <- ggplot(pca_data_cleaned, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = "dashed", linewidth = 0.8) +
  geom_text_repel(
    aes(label = Sample),
    size = 3, box.padding = 0.5,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("Tumor" = "#F55A5A", "Normal" = "#65B2FA"),
    labels = c(paste0("Normal (n=", sum(pca_data_cleaned$Group == "Normal"), ")"), 
               paste0("Tumor (n=", sum(pca_data_cleaned$Group == "Tumor"), ")"))
  ) +
  labs(
    title = "PCA of DNA Methylation - Final Cleaned Data",
    subtitle = paste0("Total: ", nrow(pca_data_cleaned), 
                      " (Tumor: ", sum(pca_data_cleaned$Group == "Tumor"), 
                      ", Normal: ", sum(pca_data_cleaned$Group == "Normal"), ")"),
    x = paste0("PC1 (", var_explained_cleaned[1], "% Variance)"),
    y = paste0("PC2 (", var_explained_cleaned[2], "% Variance)")
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

print(pca_plot_final)

#=============================================================================
# 8. Prepare filtered data for DMP analysis
#=============================================================================
filtered_samples <- rownames(pca_input_cleaned)
beta_matrix_filtered <- beta_matrix[, filtered_samples]
phenoData_filtered <- phenoData[filtered_samples, , drop = FALSE]

#=============================================================================
# 9. DMP analysis
#=============================================================================
DMP_filtered <- champ.DMP(
  beta = beta_matrix_filtered,
  pheno = phenoData_filtered$group,
  compare.group = c("Tumor", "Normal"),
  adjust.method = "BH",
  adjPVal = 1,
  arraytype = "EPIC"
)

DMP_results <- DMP_filtered$Tumor_to_Normal

# Correct deltaBeta sign (champ.DMP computes Normal - Tumor by default)
DMP_results$deltaBeta <- -DMP_results$deltaBeta

# Filter significant DMPs
significant_DMPs <- subset(DMP_results, 
                           abs(deltaBeta) > 0.2 & adj.P.Val < 0.05)

# Deduplicate by gene (keep first probe per gene)
significant_DMPs_unique <- significant_DMPs[!duplicated(significant_DMPs$gene), ]
significant_DMPs_unique <- significant_DMPs_unique[significant_DMPs_unique$gene != "", ]

final_DMP_results <- significant_DMPs_unique

#=============================================================================
# 10. Export results
#=============================================================================
write.csv(final_DMP_results, file = "DNA_methylation_DMPs_gene.csv", row.names = TRUE)

# Statistics
upgene <- subset(final_DMP_results, deltaBeta > 0.2 & adj.P.Val < 0.05)
downgene <- subset(final_DMP_results, deltaBeta < -0.2 & adj.P.Val < 0.05)

cat("\n=== DMP Statistics ===\n")
cat("Total significant DMPs:", nrow(significant_DMPs), "\n")
cat("Unique genes after dedup:", nrow(significant_DMPs_unique), "\n")
cat("Hypermethylated (deltaBeta > 0.2):", nrow(upgene), "\n")
cat("Hypomethylated (deltaBeta < -0.2):", nrow(downgene), "\n")
cat("Final output genes:", nrow(final_DMP_results), "\n")

#=============================================================================
# 11. Heatmap - Top 100 probes
#=============================================================================
sig_genes_sorted <- final_DMP_results[order(-abs(final_DMP_results$deltaBeta)), ]
top100_probes <- head(sig_genes_sorted, 100)

heatmap_data <- as.matrix(beta_matrix_filtered[rownames(top100_probes), ])
rownames(heatmap_data) <- top100_probes$gene

annotation_col <- data.frame(Group = phenoData_filtered$group)
rownames(annotation_col) <- colnames(heatmap_data)

group_order <- order(phenoData_filtered$group)
heatmap_data_ordered <- heatmap_data[, group_order]
annotation_col_ordered <- annotation_col[group_order, , drop = FALSE]

group_colors <- list(Group = c(Normal = "#1F77B4", Tumor = "#FF7F0E"))

pheatmap(
  heatmap_data_ordered,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = annotation_col_ordered,
  annotation_colors = group_colors,
  color = colorRampPalette(c("#65B2FA", "white", "#F55A5A"))(100),
  main = "Top 100 DM Probes After Sample Removal",
  border_color = NA,
  treeheight_row = 17
)

#=============================================================================
# 12. Correlation heatmap
#=============================================================================
cor_matrix <- cor(as.matrix(beta_matrix_filtered), method = "pearson")

group_counts <- table(phenoData_filtered$group)
gaps_col <- cumsum(group_counts)[-length(group_counts)]

pheatmap(
  cor_matrix,
  annotation_col = annotation_col,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  gaps_col = gaps_col,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Pearson Correlation Heatmap - All Probes"
)

#=============================================================================
# 13. Volcano plot
#=============================================================================
DMP_results$Significance <- "Not significant"
DMP_results$Significance[DMP_results$deltaBeta < -0.2 & DMP_results$adj.P.Val < 0.05] <- "Hypomethylated"
DMP_results$Significance[DMP_results$deltaBeta > 0.2 & DMP_results$adj.P.Val < 0.05] <- "Hypermethylated"

# Top genes for labeling
top_hyper <- DMP_results %>%
  filter(Significance == "Hypermethylated") %>%
  arrange(P.Value) %>%
  head(10)

top_hypo <- DMP_results %>%
  filter(Significance == "Hypomethylated") %>%
  arrange(P.Value) %>%
  head(10)

top_genes <- rbind(top_hyper, top_hypo)

my_colors <- c("Hypomethylated" = "#65B2FA",
               "Not significant" = "#D3D3D3",
               "Hypermethylated" = "#F55A5A")

volcano_plot <- ggplot(DMP_results, aes(x = deltaBeta, y = -log10(P.Value))) +
  geom_point(aes(color = Significance, fill = Significance), 
             alpha = 0.8, size = 2.5, shape = 21, stroke = 0.3) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
             color = "#666666", alpha = 0.8, linewidth = 0.6) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed", 
             color = "#666666", alpha = 0.8, linewidth = 0.6) +
  geom_text_repel(data = top_genes, 
                  aes(label = gene, color = Significance),
                  size = 3.2, fontface = "bold",
                  max.overlaps = 25,
                  box.padding = 0.6,
                  point.padding = 0.3,
                  show.legend = FALSE) +
  labs(
    x = "Methylation Difference (deltaBeta)", 
    y = expression(-log[10](P-value)),
    title = "Differential DNA Methylation Analysis - After Sample Removal"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )

print(volcano_plot)

cat("\n=== Analysis Complete ===", Sys.time(), "\n")