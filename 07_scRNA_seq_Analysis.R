# ===========================================================================
# Script: 07_scRNA_seq_Analysis.R
# Project: Thymoma Multi-Omics Biomarker Discovery
# Description: Single-cell RNA-seq validation
#              - QC, normalization, Harmony integration, clustering, UMAP/tSNE
#              - Cell type annotation (Pre B cell, Cytotoxic T, DP cell,
#                Macrophage, B cell, Tumor cell)
#              - Marker discovery and target-gene feature plots
#                (FOXN1, PKP1, IRF6, TRIP6, CDH1, PBX1, KRT1)
#              - scTenifoldKnk in-silico PKP1 knockout on the tumor cluster
# Date: 2026-09-05
# ===========================================================================
# 1. Load data and create Seurat object
rm(list = ls())
getwd()
setwd("/public1/home/l0s001024/Normal/")

# Load expression matrix
cc <- read.csv("expression.csv")
save(cc, file = "scRNA_seq.Rds")

library(Seurat)
rownames(cc) <- cc$X
cc <- cc[, -1]

# Filter cells by sample names
table(scobj$seurat_clusters)
new_cc <- cc[, colnames(cc) %in% name]

# Create initial Seurat object
final_seurat <- CreateSeuratObject(
  counts = new_cc,
  project = "scRNA_PROJECT",
  min.cells = 3,
  min.features = 200
)

table(final_seurat@meta.data$orig.ident)
table(final_seurat$orig.ident)

# Alternative object creation
new_seurat <- CreateSeuratObject(
  counts = dd,
  project = "scRNA_PROJECT",
  min.cells = 3,
  min.features = 200
)

colnames(seurat@meta.data)
table(seurat$orig.ident)

#=============================================================================
# 2. Load required packages
#=============================================================================
library(harmony)
library(ggsci)
library(dplyr) 
library(future)
library(Seurat)
library(clustree)
library(cowplot)
library(data.table)
library(ggplot2)
library(patchwork)
library(stringr)
library(Matrix)

#=============================================================================
# 3. Load and prepare Seurat object
#=============================================================================
getwd()
setwd("/public1/home/l0s001024/Normal/")
load("final_seurat.RData")

nrow(final_seurat)

# Remove duplicate gene names (suffix pattern)
new_names <- sub("\\.[0-9]+$", "", rownames(final_seurat))
sum(duplicated(new_names))

name <- new_names[!duplicated(new_names)]
table(duplicated(name))
final <- final_seurat[name, ]

scobj <- final

str(scobj)
gc()
table(scobj$orig.ident)
nrow(scobj)
ncol(scobj)

#=============================================================================
# 4. Quality Control
#=============================================================================
# Calculate mitochondrial and ribosomal percentages
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
scobj[["percent.rb"]] <- PercentageFeatureSet(scobj, pattern = "^RP[SL]")

# Visualize QC metrics
VlnPlot(scobj, features = c("nFeature_RNA", "percent.mt", "percent.rb", "nCount_RNA"), 
        ncol = 4, pt.size = 0, group.by = "orig.ident")

# Remove ribosomal and mitochondrial genes
ribosomal_genes <- grep("^RP[SL]", rownames(scobj), value = TRUE)
mitochondrial_genes <- grep("^MT-", rownames(scobj), value = TRUE)
genes_to_remove <- c(ribosomal_genes, mitochondrial_genes)
scobj <- scobj[!rownames(scobj) %in% genes_to_remove, ]

print(dim(scobj))

# QC thresholds
quantile(scobj$nFeature_RNA, seq(0.05, 0.95, 0.05))
quantile(scobj$nCount_RNA, seq(0.05, 0.95, 0.05))
quantile(scobj$percent.mt, seq(0.05, 0.95, 0.05))
quantile(scobj$percent.rb, seq(0.05, 0.95, 0.05))

# Filter cells based on QC thresholds
scobj <- subset(scobj, 
                subset = nFeature_RNA > 574.0 & 
                  nFeature_RNA < 4728.8 & 
                  nCount_RNA < 3659.517 & 
                  percent.mt < 2.0819036 & 
                  percent.rb < 13.971344)

#=============================================================================
# 5. Normalization and dimensionality reduction
#=============================================================================
scobj <- NormalizeData(scobj, normalization.method = "LogNormalize", scale.factor = 10000)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)
scobj <- ScaleData(scobj, features = rownames(scobj))

# PCA
scobj <- RunPCA(scobj, features = VariableFeatures(scobj), reduction.name = "pca")
ElbowPlot(scobj)

# Harmony batch correction
scobj <- RunHarmony(scobj, reduction = "pca", 
                    group.by.vars = "orig.ident", 
                    reduction.save = "harmony")

# tSNE and UMAP
scobj <- RunTSNE(scobj, dims = 1:20, reduction = "pca", reduction.name = "tsne")

scobj <- FindNeighbors(scobj, reduction = "pca", dims = 1:20)
scobj <- FindClusters(scobj, resolution = 0.03)
scobj <- RunUMAP(scobj, reduction = "pca", dims = 1:20, reduction.name = "umap")

Idents(scobj) <- "seurat_clusters"
scobj <- JoinLayers(scobj, name = "RNA")

#=============================================================================
# 6. Visualization of clusters
#=============================================================================
DimPlot(scobj, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(scobj, reduction = "tsne", label = TRUE, repel = TRUE)
table(scobj$seurat_clusters)

DimPlot(scobj, reduction = "umap", group.by = "orig.ident", label = TRUE)

#=============================================================================
# 7. Find marker genes
#=============================================================================
markers <- FindAllMarkers(scobj, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.5)
table(markers$cluster)
table(scobj$seurat_clusters)
table(markers$cluster)

# save(scobj, file = "scobj.RDS")

# Top 10 markers per cluster
library(dplyr)
top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

top10 <- top10$gene
table(scobj$seurat_clusters)
DoHeatmap(scobj, features = top10, size = 4, label = TRUE)

#=============================================================================
# 8. Cell type annotation
#=============================================================================
# Dot plot for marker genes
d <- DotPlot(scobj, features = c("VPREB1", "IGLL1",
                                 "NKG7", "IFNG", "GZMA", "GZMK",
                                 "RAG1", "RAG2",
                                 "CD68", "CD14",
                                 "CD79A", "CD79B", "CD19", "MS4A1",
                                 "ADH1B", "PDGFRA")) +
  RotatedAxis() +
  scale_size(range = c(0, 100)) +
  scale_size_continuous(limits = c(0, 100)) +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(axis.text.x = element_text(size = 8),
        axis.test.y = element_text(size = 5))

# Rename clusters
new.cluster.ids <- c("Pre B cell", "Cytotoxic T", "DP cell", 
                     "Macrophage", "B cell", "Tumor cell")
names(new.cluster.ids) <- levels(scobj)
scobj <- RenameIdents(scobj, new.cluster.ids)
table(scobj$seurat_clusters)
saveRDS(scobj, file = "scobj.rds")

#=============================================================================
# 9. Feature plots for target genes
#=============================================================================
FeaturePlot(scobj, features = "FOXN1", pt.size = 0.4, cols = c("lightgrey", "blue"))
FeaturePlot(scobj, features = "PKP1", pt.size = 0.4, cols = c("lightgrey", "blue"))
FeaturePlot(scobj, features = "IRF6", pt.size = 0.4, cols = c("lightgrey", "blue"))
FeaturePlot(scobj, features = "TRIP6", pt.size = 0.4, cols = c("lightgrey", "blue"))
FeaturePlot(scobj, features = "CDH1", pt.size = 0.4, cols = c("lightgrey", "blue"))
FeaturePlot(scobj, features = "PBX1", pt.size = 0.4, cols = c("lightgrey", "blue"))
FeaturePlot(scobj, features = "KRT1", pt.size = 0.4, cols = c("lightgrey", "blue"))

# Dot plot for final genes
pdf("finalgene_dotplot.pdf", width = 8, height = 8)
DotPlot(scobj, features = c("FOXN1", "IRF6", "TRIP6", "PKP1", "CDH1", "PBX1", "KRT1")) +
  RotatedAxis() +
  scale_size(range = c(0, 100)) +
  scale_size_continuous(limits = c(0, 100)) +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  theme(axis.text.x = element_text(size = 8),
        axis.test.y = element_text(size = 5))
dev.off()

#=============================================================================
# 10. Extract cluster 5 for perturbation analysis
#=============================================================================
cluster5 <- subset(scobj, subset = seurat_clusters == 5)
class(cluster5)
saveRDS(cluster5, file = "cluster5.rds")

#=============================================================================
# 11. PKP1 perturbation analysis using scTenifoldKnk
#=============================================================================
library(scTenifoldKnk)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(dplyr)

set.seed(1234)
getwd()
setwd("/public1/home/l0s001024/Normal")
pbmc <- readRDS("cluster5.rds")

# Extract count matrix
scRNAseq = as.matrix(GetAssayData(object = pbmc@assays$RNA, layer = "counts"))
scRNAseq = as.data.frame(scRNAseq)
dim(scRNAseq)

# Filter by target genes
target_gene <- read.csv("6cluster_marker.csv")
target_gene <- subset(target_gene, cluster == 5)
getwd()

length(target_gene$gene)
scRNAseq <- scRNAseq[target_gene$gene, ]

# Filter genes with expression in at least 6% of cells
expr_ratio <- rowSums(scRNAseq > 0) / ncol(scRNAseq)
scRNAseq <- scRNAseq[expr_ratio >= 0.06, ]
dim(scRNAseq)

# Run perturbation analysis
PKP1_expression <- scRNAseq["PKP1", ]
PKP1_expression <- as.data.frame(t(PKP1_expression))

results = scTenifoldKnk(countMatrix = scRNAseq, gKO = 'PKP1')

# Extract differential regulation results
generesults = results[["diffRegulation"]]
generesults$log2FC = log2(generesults$FC)

diffgeneresults = generesults[generesults$p.adj < 0.05, ]
write.csv(diffgeneresults, "final_result.csv")
diffgeneresults <- diffgeneresults[diffgeneresults$gene != "PKP1", ]

#=============================================================================
# 12. Visualization - Z-score bar plot
#=============================================================================
pdf('type3_scTenifoldKnk.pdf', height = 8.5, width = 8)
ggplot(diffgeneresults, aes(x = reorder(gene, Z), y = Z)) +
  geom_bar(stat = 'identity', aes(fill = Z > 0), width = 0.7) +
  scale_fill_manual(values = c("#4198AC", "#ED8D5A")) +
  coord_flip() +
  labs(x = "Gene", y = "Z-score") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.8)
  )
dev.off()

write.csv(generesults, file = "scTenifoldKnk_results.csv")

#=============================================================================
# 13. Visualization - FC bar plot
#=============================================================================
library(dplyr)
library(ggplot2)
library(forcats)
library(ggrepel)
library(glue)

list.files()
results <- read.csv("final_result.csv", row.names = 1)
top_genes <- head(results[order(-results$FC), ], 86)
top_genes <- top_genes[-1, ]

p_bar <- top_genes %>%
  mutate(gene = fct_reorder(gene, FC)) %>%
  ggplot(aes(x = gene, y = FC, fill = FC)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  coord_flip() +
  geom_text(aes(label = round(FC, 1)),
            hjust = -0.15, size = 3.6, color = "grey20") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  scale_fill_gradient(low = "#9ecae1", high = "#08519c") +
  labs(
    title = "Differentially Regulated Genes",
    subtitle = paste0("Ranked by fold change (n = ", nrow(top_genes), ")"),
    x = NULL, y = "Fold Change (FC)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "grey40"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.8)
  )

print(p_bar)

#=============================================================================
# 14. Visualization - Volcano plot (Z-score vs p-value)
#=============================================================================
df <- results
df <- df[-1, ]

pretty_z_p_plot <- function(df,
                            z_col = "Z",
                            p_col = "p.adj",
                            gene_col = "gene",
                            z_cut = 2,
                            p_cut = 0.01,
                            p_line = 0.05,
                            max_label = 60) {
  
  df <- df %>%
    mutate(
      log_pval = -log10(pmax(.data[[p_col]], 1e-300)),
      group = case_when(
        .data[[p_col]] < p_cut & abs(.data[[z_col]]) >= z_cut ~ "Hit (|Z|>=2 & adj.p<0.01)",
        abs(.data[[z_col]]) >= z_cut ~ "|Z|>=2 only",
        .data[[p_col]] < p_cut ~ "adj.p<0.01 only",
        TRUE ~ "Not significant"
      )
    )
  
  label_df <- df %>%
    filter(.data[[p_col]] < p_cut, abs(.data[[z_col]]) >= z_cut)
  
  pal <- c("Not significant" = "#B8B8B8",
           "|Z|>=2 only"     = "#7AA6C2",
           "adj.p<0.01 only"    = "#C49A6C",
           "Hit (|Z|>=2 & adj.p<0.01)" = "#E15759")
  
  y_thr <- -log10(p_cut)
  x_thr <- z_cut
  
  p <- ggplot(df, aes(x = .data[[z_col]], y = log_pval, color = group)) +
    annotate("rect",
             xmin =  x_thr, xmax =  Inf,
             ymin =  y_thr, ymax =  Inf,
             alpha = 0.06, fill = "#E15759") +
    annotate("rect",
             xmin = -Inf, xmax = -x_thr,
             ymin =  y_thr, ymax =  Inf,
             alpha = 0.06, fill = "#E15759") +
    geom_point(size = 1.9, alpha = 0.75, stroke = 0) +
    geom_hline(yintercept = -log10(p_line), linetype = "dashed",
               color = "#E15759", linewidth = 0.6) +
    geom_hline(yintercept =  y_thr, linetype = "dotted",
               color = "#E15759", linewidth = 0.6) +
    geom_vline(xintercept =  c(-x_thr, x_thr), linetype = "dashed",
               color = "#2B6CB0", linewidth = 0.6) +
    ggrepel::geom_text_repel(
      data = head(label_df %>% arrange(desc(log_pval)), max_label),
      aes(label = .data[[gene_col]]),
      size = 3.3, seed = 123,
      box.padding = 0.35, point.padding = 0.2,
      max.overlaps = Inf, min.segment.length = 0,
      segment.color = alpha("grey20", 0.6), segment.size = 0.3
    ) +
    labs(
      title = "Z-score vs -log10(adjusted p-value)",
      subtitle = glue("Thresholds: |Z| >= {z_cut}, adjusted p-value < {p_cut}  .  Hits = {nrow(label_df)}"),
      x = "Z-score",
      y = expression(-log[10](adjusted~p-value))
    ) +
    scale_color_manual(values = pal, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title.position = "plot",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      axis.title = element_text(face = "bold"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.1, "cm"),
      axis.ticks.x = element_line(color = "black", linewidth = 0.8)
    )
  
  return(p)
}

p_scatter <- pretty_z_p_plot(df,
                             z_col = "Z",
                             p_col = "p.adj",
                             gene_col = "gene",
                             z_cut = 2,
                             p_cut = 0.01,
                             p_line = 0.05,
                             max_label = 10)

print(p_scatter)

#=============================================================================
# 15. Save plots
#=============================================================================
pdf("volcanic.pdf", width = 7, height = 7)
print(p_scatter)
dev.off()

pdf("allgene_FC_barplot.pdf", width = 10, height = 11)
print(p_bar)
dev.off()

#=============================================================================
# 16. Top 20 genes Z-score bar plot
#=============================================================================
diffgeneresults <- top_genes
diffgeneresults <- diffgeneresults[c(1:20), ]

p <- ggplot(diffgeneresults, aes(x = reorder(gene, Z), y = Z)) +
  geom_bar(stat = 'identity', aes(fill = Z > 0), width = 0.7) +
  scale_fill_manual(values = c("#65B2FA", "#ED8D5A")) +
  coord_flip(ylim = c(min(diffgeneresults$Z) * 1, max(diffgeneresults$Z) * 1)) +
  labs(
    title = "Top 20 Regulated Genes",
    x = "Gene",
    y = "Z score"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 10),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.1, "cm"),
    axis.ticks.x = element_line(color = "black", linewidth = 0.8)
  )

pdf("top20_regulated_genes.pdf", width = 10, height = 9.5)
print(p)
dev.off()

cat("\n=== Analysis Complete ===", Sys.time(), "\n")