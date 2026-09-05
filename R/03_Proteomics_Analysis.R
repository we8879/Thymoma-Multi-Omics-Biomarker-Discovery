# ===========================================================================
# Script: 03_Proteomics_Analysis.R
# Project: Thymoma Multi-Omics Biomarker Discovery
# Description: Proteomic differential expression and functional interpretation
#              - PCA (before DE)
#              - limma moderated t-test (Tumor vs Normal)
#              - Volcano plot, DEG heatmap, correlation heatmap, PCA
#              - GO + KEGG enrichment for up- and down-regulated proteins
# Date: 2026-09-05
# ===========================================================================

# 1. Load and prepare protein data
rm(list = ls())
library(readxl)
library(stringr)

getwd()
setwd("D:/Thymoma_MultiOmics/data/proteomics/")

# Load protein expression data
protein <- read_excel("proteomics.xlsx")
range(protein$`squamous 7_THYMOMA_61`)
protein <- as.data.frame(protein)
rownames(protein) <- protein$geneSymbol
protein <- protein[,-c(1,2)]
colnames(protein)

# Remove unwanted columns
protein <- protein[,-22]
protein <- protein[,-c(16,29)]
protein <- protein[,-21]

colnames(protein)
ncol(protein)

# ===========================================================================
# 2. PCA analysis 
# ===========================================================================
library(ggplot2)
library(ggrepel)

# Prepare PCA data
pca_data <- t(protein)
pca_data <- pca_data[, colSums(is.na(pca_data)) == 0]
col_var <- apply(pca_data, 2, var, na.rm = TRUE)
pca_data_filtered <- pca_data[, col_var > 0 & !is.na(col_var)]
pca_data <- pca_data_filtered

# Define sample groups
cc <- str_detect(colnames(protein), "para")
pro_group <- ifelse(cc, "Normal", "Tumor")
print(pro_group)
table(pro_group)

# Perform PCA
pca_result <- prcomp(pca_data, scale. = TRUE)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$Group <- pro_group
pca_scores$Sample <- rownames(pca_scores)

variance_explained <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 2)

# PCA with confidence ellipses
pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  stat_ellipse(type = "t", linetype = "dashed", size = 0.8, alpha = 0.8, fill = NA) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Normal" = "#65B2FA", "Tumor" = "#F55A5A")) +
  scale_fill_manual(values = c("Normal" = "#65B2FA", "Tumor" = "#F55A5A")) +
  labs(x = paste0("PC1 (", variance_explained[1], "%)"),
       y = paste0("PC2 (", variance_explained[2], "%)"),
       title = "PCA Plot: Tumor vs Normal Samples",
       color = "Group") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )

print(pca_plot)

# ===========================================================================
# 3. Differential expression analysis 
# ===========================================================================
protein$gene <- rownames(protein)
protein <- protein[,-27]

library(limma)

# Design matrix
design <- model.matrix(~0 + pro_group)
colnames(design) <- c("Normal", "Tumor")

# Fit linear model
fit <- lmFit(protein, design)

# Contrast matrix for Tumor vs Normal
contrast.matrix <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract results
limma_results <- topTable(fit2, coef = "Tumor_vs_Normal", number = Inf, adjust.method = "BH")
head(limma_results, 10)

# Filter significant proteins
upgene <- subset(limma_results, adj.P.Val < 0.05 & logFC > 1)
downgene <- subset(limma_results, adj.P.Val < 0.05 & logFC < -1)

nrow(upgene)
nrow(downgene)

protein_gene <- subset(limma_results, adj.P.Val < 0.05 & abs(logFC) > 1)
protein_gene$gene <- rownames(protein_gene)
nrow(protein_gene)

# Export results
# write.csv(protein_gene, file = "proteomics_DEGs.csv")

# ===========================================================================
# 4. Volcano plot
# ===========================================================================
library(ggplot2)
library(ggrepel)
library(dplyr)

res <- as.data.frame(limma_results)
res <- res %>% rename(p.adjust = adj.P.Val)
colnames(res)[1] <- "log2FoldChange"

log2FC_threshold <- 1
p_adj_threshold <- 0.05

res$gene <- rownames(res)

# Classify genes by expression change
res <- res %>%
  mutate(
    expression = case_when(
      is.na(p.adjust) | is.na(log2FoldChange) ~ "Not significant",
      p.adjust < p_adj_threshold & log2FoldChange > log2FC_threshold ~ "Up",
      p.adjust < p_adj_threshold & log2FoldChange < -log2FC_threshold ~ "Down",
      TRUE ~ "Not significant"
    ),
    expression = factor(expression, levels = c("Up", "Down", "Not significant"))
  )

# Volcano plot
volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = -log10(p.adjust))) +
  geom_point(aes(color = expression, alpha = expression), size = 2) +
  scale_color_manual(values = c("Up" = "#F55A5A",
                                "Down" = "#65B2FA",
                                "Not significant" = "#B3B3B3")) +
  scale_alpha_manual(values = c("Up" = 0.8, "Down" = 0.8, "Not significant" = 0.4),
                     guide = "none") +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold),
             linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(p_adj_threshold),
             linetype = "dashed", color = "black", alpha = 0.5)

sig_genes <- res %>% filter(expression != "Not significant")
if (nrow(sig_genes) > 0) {
  volcano_plot <- volcano_plot +
    geom_text_repel(
      data = sig_genes %>% arrange(p.adjust) %>% head(min(10, nrow(sig_genes))),
      aes(label = gene),
      size = 3, max.overlaps = 20, box.padding = 0.5,
      segment.color = "grey50", segment.size = 0.2
    )
}

volcano_plot <- volcano_plot +
  labs(
    title = "Volcano Plot - Tumor vs Normal",
    x = expression(log[2] * "(Fold Change)"),
    y = expression(-log[10] * "(p.adjust)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(size = 12),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(size = 10, color = "black"),
    legend.title = element_text(face = "bold")
  ) +
  scale_x_continuous(
    breaks = seq(-30, 30, by = 2),
    limits = c(-max(abs(res$log2FoldChange), na.rm = TRUE) - 1,
               max(abs(res$log2FoldChange), na.rm = TRUE) + 1)
  )

print(volcano_plot)

# ===========================================================================
# 5. Heatmap - Top 100 differentially expressed proteins
# ===========================================================================
library(pheatmap)

top_100_genes <- limma_results %>%
  arrange(desc(abs(logFC))) %>%
  head(100)

heatmap_data <- protein[rownames(top_100_genes), ]
heatmap_data_scaled <- t(scale(t(heatmap_data)))

annotation_col <- data.frame(Group = pro_group)
rownames(annotation_col) <- colnames(heatmap_data_scaled)

annotation_colors <- list(Group = c(Normal = "#1F77B4", Tumor = "#FF7F0E"))

pheatmap(heatmap_data_scaled,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("#65B2FA", "white", "#F55A5A"))(100),
         fontsize_row = 6,
         fontsize_col = 8,
         main = "Top 100 Differentially Expressed Proteins",
         border_color = NA,
         treeheight_row = 15,
         treeheight_col = 11)

# ===========================================================================
# 6. PKP1 expression bar plot
# ===========================================================================
library(ggplot2)
library(tidyr)

pkp1_data <- protein["PKP1", ]
pkp1_long <- data.frame(
  Sample = colnames(pkp1_data),
  Expression = as.numeric(pkp1_data)
)

ggplot(pkp1_long, aes(x = Sample, y = Expression, fill = Sample)) +
  geom_col(color = "black", alpha = 0.8) +
  geom_text(aes(label = round(Expression, 2)), vjust = -0.5, size = 4) +
  labs(title = "PKP1 Gene Expression Across Samples",
       x = "Samples",
       y = "Expression Level") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# ===========================================================================
# 7. Sample correlation heatmap
# ===========================================================================
sample_correlation <- cor(protein, use = "complete.obs", method = "pearson")

annotation_col <- data.frame(Group = pro_group)
rownames(annotation_col) <- colnames(protein)

annotation_colors <- list(Group = c(Normal = "#1F77B4", Tumor = "#FF7F0E"))

normal_samples <- which(pro_group == "Normal")
tumor_samples <- which(pro_group == "Tumor")
sample_order <- c(normal_samples, tumor_samples)

sample_correlation_ordered <- sample_correlation[sample_order, sample_order]
annotation_col_ordered <- annotation_col[sample_order, , drop = FALSE]

pheatmap(sample_correlation_ordered,
         annotation_col = annotation_col_ordered,
         annotation_colors = annotation_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("#65B2FA", "white", "#F55A5A"))(100),
         show_rownames = TRUE,
         show_colnames = TRUE,
         display_numbers = FALSE,
         number_color = "black",
         fontsize_number = 6,
         main = "Sample Correlation Heatmap (Grouped)",
         fontsize_row = 8,
         fontsize_col = 8,
         border_color = NA)

# ===========================================================================
# 8. Functional enrichment analysis - Up-regulated genes
# ===========================================================================
rm(list = ls())
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)

R.utils::setOption("clusterProfiler.download.method", "auto")
setwd("D:/Thymoma_MultiOmics/data/proteomics/")

pvalueFilter <- 0.05
qvalueFilter <- 0.05
colorSel <- "qvalue"
if (qvalueFilter > 0.05) colorSel <- "pvalue"
ontology.col <- c("#00AFBB", "#E7B800", "#90EE90")

# Read protein data
rt <- read.csv("protein_gene.csv")
rownames(rt) <- rt$gene

# Up-regulated genes only
rt <- rt[rt$logFC > 0, ]

# Convert gene symbols to Entrez IDs
genes <- unique(as.vector(rt[,1]))
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezIDs)
gene <- entrezIDs[entrezIDs != "NA"]
gene <- gsub("c\\(\"(\\d+)\".*", "\\1", gene)

if (length(gene) == 0) stop("Gene list is empty, please check input")

# GO enrichment
kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, 
               qvalueCutoff = 0.05, ont = "all", readable = TRUE)

GO <- as.data.frame(kk)
GO <- GO[(GO$pvalue < pvalueFilter & GO$qvalue < qvalueFilter), ]
write.csv(GO, file = "up_protein_GO.csv")

# GO bar plot
showNum <- min(10, nrow(GO))
pdf(file = "GObarplot.pdf", width = 10, height = 7)
bar <- barplot(kk, drop = TRUE, showCategory = showNum, label_format = 130, 
               split = "ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bar)
dev.off()

# GO bubble plot
pdf(file = "GObubble.pdf", width = 10, height = 7)
bub <- dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", 
               label_format = 130, split = "ONTOLOGY", color = colorSel) + 
  facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bub)
dev.off()

# GO circle plot
data <- GO[order(GO$pvalue), ]
datasig <- data[data$pvalue < 0.05, , drop = FALSE]
BP <- head(datasig[datasig$ONTOLOGY == "BP", , drop = FALSE], 6)
CC <- head(datasig[datasig$ONTOLOGY == "CC", , drop = FALSE], 6)
MF <- head(datasig[datasig$ONTOLOGY == "MF", , drop = FALSE], 6)
data <- rbind(BP, CC, MF)

main.col <- ontology.col[as.numeric(as.factor(data$ONTOLOGY))]
BgGene <- as.numeric(sapply(strsplit(data$BgRatio, "/"), '[', 1))
Gene <- as.numeric(sapply(strsplit(data$GeneRatio, '/'), '[', 1))
ratio <- Gene / BgGene
logpvalue <- -log(data$pvalue, 10)
logpvalue.col <- brewer.pal(n = 8, name = "Reds")
f <- colorRamp2(breaks = c(0, 2, 4, 6, 8, 10, 15, 20), colors = logpvalue.col)
BgGene.col <- f(logpvalue)

df <- data.frame(GO = data$ID, start = 1, end = max(BgGene))
rownames(df) <- df$GO
bed2 <- data.frame(GO = data$ID, start = 1, end = BgGene, BgGene = BgGene, BgGene.col = BgGene.col)
bed3 <- data.frame(GO = data$ID, start = 1, end = Gene, BgGene = Gene)
bed4 <- data.frame(GO = data$ID, start = 1, end = max(BgGene), ratio = ratio, col = main.col)
bed4$ratio <- bed4$ratio / max(bed4$ratio) * 9.5

pdf("up_gene_GO_circlize.pdf", width = 10, height = 10)
par(omi = c(0.1, 0.1, 0.1, 1.5))
circos.genomicInitialize(df, plotType = "none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, 
              facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA, bg.col = main.col)

for (si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si, track.index = 1,
              major.at = seq(0, max(BgGene), by = 100), labels.facing = "clockwise")
}

circos.genomicTrack(bed2, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, 
                                         col = value[,2], border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], 
                                         adj = 0, cex = 0.8, ...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, 
                                         col = '#BA55D3', border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[,1], 
                                         cex = 0.9, adj = 0, ...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10), track.height = 0.35, bg.border = "white", bg.col = "grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim <- get.cell.meta.data("cell.xlim")
                      cell.ylim <- get.cell.meta.data("cell.ylim")
                      for (j in 1:9) {
                        y <- cell.ylim[1] + (cell.ylim[2] - cell.ylim[1]) / 10 * j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[,1], 
                                         col = value[,2], border = NA, ...)
                    })
circos.clear()

middle.legend <- Legend(labels = c('Number of Genes', 'Number of Select', 'Rich Factor(0-1)'),
                        type = "points", pch = c(15, 15, 17), 
                        legend_gp = gpar(col = c('pink', '#BA55D3', ontology.col[1])),
                        title = "", nrow = 3, size = unit(3, "mm"))
circle_size <- unit(1, "snpc")
draw(middle.legend, x = circle_size * 0.42)

main.legend <- Legend(labels = c("Biological Process", "Cellular Component", "Molecular Function"),
                      type = "points", pch = 15, legend_gp = gpar(col = ontology.col),
                      title_position = "topcenter", title = "ONTOLOGY", nrow = 3,
                      size = unit(3, "mm"), grid_height = unit(5, "mm"), grid_width = unit(5, "mm"))
logp.legend <- Legend(labels = c('(0,2]', '(2,4]', '(4,6]', '(6,8]', '(8,10]', '(10,15]', '(15,20]', '>=20'),
                      type = "points", pch = 16, legend_gp = gpar(col = logpvalue.col),
                      title = "-log10(Pvalue)", title_position = "topcenter",
                      grid_height = unit(5, "mm"), grid_width = unit(5, "mm"), size = unit(3, "mm"))
lgd <- packLegend(main.legend, logp.legend)
draw(lgd, x = circle_size * 0.85, y = circle_size * 0.55, just = "left")
dev.off()

# KEGG enrichment for up-regulated genes
options(timeout = 300)
genes <- unique(as.vector(rt[,1]))
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezIDs)
rt_entrez <- data.frame(genes, entrezID = entrezIDs)
gene <- entrezIDs[entrezIDs != "NA"]
gene <- gsub("c\\(\"(\\d+)\".*", "\\1", gene)

kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)
KEGG <- as.data.frame(kk)
KEGG$geneID <- as.character(sapply(KEGG$geneID, function(x) 
  paste(rt_entrez$genes[match(strsplit(x, "/")[[1]], as.character(rt_entrez$entrezID))], collapse = "/")))
KEGG <- KEGG[(KEGG$pvalue < pvalueFilter & KEGG$qvalue < qvalueFilter), ]
write.csv(KEGG, file = "up_protein_KEGG.csv")

showNum <- min(30, nrow(KEGG))
pdf(file = "KEGGbarplot.pdf", width = 9, height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, label_format = 130, color = colorSel)
dev.off()

pdf(file = "KEGGbubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", label_format = 130, color = colorSel)
dev.off()

# ===========================================================================
# 9. Functional enrichment analysis - Down-regulated genes
# ===========================================================================
rt <- read.csv("protein_gene.csv")
rownames(rt) <- rt$gene
rt <- rt[rt$logFC < 0, ]

genes <- unique(as.vector(rt[,1]))
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezIDs)
gene <- entrezIDs[entrezIDs != "NA"]
gene <- gsub("c\\(\"(\\d+)\".*", "\\1", gene)

if (length(gene) == 0) stop("Gene list is empty, please check input")

# GO enrichment for down-regulated genes
kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, 
               qvalueCutoff = 0.05, ont = "all", readable = TRUE)

GO <- as.data.frame(kk)
GO <- GO[(GO$pvalue < pvalueFilter & GO$qvalue < qvalueFilter), ]
write.csv(GO, file = "down_protein_GO.csv")

showNum <- min(10, nrow(GO))
pdf(file = "down_GObarplot.pdf", width = 10, height = 7)
bar <- barplot(kk, drop = TRUE, showCategory = showNum, label_format = 130, 
               split = "ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bar)
dev.off()

pdf(file = "down_GObubble.pdf", width = 10, height = 7)
bub <- dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", 
               label_format = 130, split = "ONTOLOGY", color = colorSel) + 
  facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bub)
dev.off()

# KEGG enrichment for down-regulated genes
options(timeout = 300)
genes <- unique(as.vector(rt[,1]))
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezIDs)
rt_entrez2 <- data.frame(genes, entrezID = entrezIDs)
gene <- entrezIDs[entrezIDs != "NA"]
gene <- gsub("c\\(\"(\\d+)\".*", "\\1", gene)

kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)
KEGG <- as.data.frame(kk)
KEGG$geneID <- as.character(sapply(KEGG$geneID, function(x) 
  paste(rt_entrez2$genes[match(strsplit(x, "/")[[1]], as.character(rt_entrez2$entrezID))], collapse = "/")))
KEGG <- KEGG[(KEGG$pvalue < pvalueFilter & KEGG$qvalue < qvalueFilter), ]
write.csv(KEGG, file = "down_protein_KEGG.csv")

showNum <- min(30, nrow(KEGG))
pdf(file = "down_KEGGbarplot.pdf", width = 9, height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, label_format = 130, color = colorSel)
dev.off()

pdf(file = "down_KEGGbubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio", label_format = 130, color = colorSel)
dev.off()

cat("\n=== Analysis Complete ===", Sys.time(), "\n")
