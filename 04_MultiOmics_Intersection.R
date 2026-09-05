# ===========================================================================
# Script: 04_MultiOmics_Intersection.R
# Project: Thymoma Multi-Omics Biomarker Discovery
# Description: Three-omics integration by intersection
#              - Global three-way Venn diagram (RNA-seq x Methylation x Protein)
#              - Direction-aware intersection, keeping only genes whose
#                expression change is consistent with an epigenetic driver:
#                  Up   set: RNA up   + protein up   + methylation hypomethylated
#                  Down set: RNA down + protein down + methylation hypermethylated
# Date: 2026-09-05
# ===========================================================================
# Load packages
suppressPackageStartupMessages({
  library(VennDiagram)
  library(dplyr)
  library(grid)
})

# Set working directory
setwd("D:/Thymoma_MultiOmics/results/multiomics")
outdir <- "D:/Thymoma_MultiOmics/results/multiomics"

# Setup logging
log_file <- file.path(outdir, "analysis.log")
sink(log_file, split = TRUE)
cat("=== Analysis Started ===", Sys.time(), "\n")

# Main analysis function
main <- function() {
  cat("\nReading data...\n")
  
  # Load data files
  RNA_seq <- read.csv("D:/Thymoma_MultiOmics/results/multiomics/RNAseq_DEGs.csv")
  methylation <- read.csv("D:/Thymoma_MultiOmics/results/multiomics/DNA_methylation_DMPs_gene.csv", row.names = 1)
  protein <- read.csv("D:/Thymoma_MultiOmics/results/multiomics/proteomics_DEGs.csv")
  
  # Filter differentially expressed genes
  cat("Filtering DEGs...\n")
  RNA_seq <- filter(RNA_seq, abs(logFC) > 1 & adj.P.Val < 0.05)
  
  methylation <- methylation[methylation$gene != "", ]
  methylation <- filter(methylation, abs(deltaBeta) > 0.2 & adj.P.Val < 0.05)
  
  protein <- filter(protein, abs(logFC) > 1 & adj.P.Val < 0.05)
  
  # Extract gene names
  genes_RNA <- as.character(RNA_seq[, 1])
  genes_methylation <- as.character(methylation$gene)
  genes_protein <- as.character(protein$X)
  
  cat("Gene counts: RNA=", length(genes_RNA), 
      ", Methylation=", length(genes_methylation), 
      ", Protein=", length(genes_protein), "\n")
  
  # Triple Venn diagram
  cat("Drawing Venn diagram...\n")
  venn.plot <- draw.triple.venn(
    area1 = length(genes_RNA),
    area2 = length(genes_methylation),
    area3 = length(genes_protein),
    n12 = length(intersect(genes_RNA, genes_methylation)),
    n23 = length(intersect(genes_methylation, genes_protein)),
    n13 = length(intersect(genes_RNA, genes_protein)),
    n123 = length(Reduce(intersect, list(genes_RNA, genes_methylation, genes_protein))),
    category = c("RNA-seq", "Methylation", "Protein"),
    fill = c("skyblue", "pink", "lightgreen")
  )
  
  png(file.path(outdir, "venn_three_omics.png"), width = 800, height = 600)
  grid.draw(venn.plot)
  dev.off()
  
  # Export three-way intersection
  intersect_genes <- Reduce(intersect, list(genes_RNA, genes_methylation, genes_protein))
  write.csv(intersect_genes, file.path(outdir, "three_omics_intersect_genes.csv"), row.names = FALSE)
  
  # Directional grouping analysis
  cat("Directional analysis...\n")
  protein_up <- filter(protein, logFC > 1 & adj.P.Val < 0.05)
  protein_down <- filter(protein, logFC < -1 & adj.P.Val < 0.05)
  
  RNA_seq_up <- filter(RNA_seq, logFC > 1 & adj.P.Val < 0.05)
  RNA_seq_down <- filter(RNA_seq, logFC < -1 & adj.P.Val < 0.05)
  
  methylation_up <- filter(methylation, deltaBeta > 0.2 & adj.P.Val < 0.05)
  methylation_down <- filter(methylation, deltaBeta < -0.2 & adj.P.Val < 0.05)
  
  # Up-regulated intersection (Protein/RNA up + Methylation down)
  protein_up_genes <- protein_up$X
  rna_up_genes <- RNA_seq_up$X
  methylation_down_genes <- methylation_down$gene
  
  intersection1 <- intersect(intersect(protein_up_genes, rna_up_genes), methylation_down_genes)
  
  # Down-regulated intersection (Protein/RNA down + Methylation up)
  protein_down_genes <- protein_down$X
  rna_down_genes <- RNA_seq_down$X
  methylation_up_genes <- methylation_up$gene
  
  intersection2 <- Reduce(intersect, list(protein_down_genes, rna_down_genes, methylation_up_genes))
  
  # Merge final gene list
  intersect_final <- c(intersection1, intersection2)
  write.csv(intersect_final, file.path(outdir, "final_26gene.csv"), row.names = FALSE)
  
  # Up-regulated gene Venn diagram
  cat("Drawing up-regulated Venn diagram...\n")
  up_gene_sets <- list(
    Protein = protein_up_genes,
    RNA = rna_up_genes,
    Methylation = methylation_down_genes
  )
  
  venn_up <- venn.diagram(
    x = up_gene_sets,
    category.names = c("Protein", "RNA_seq", "Methylation"),
    filename = NULL,
    output = TRUE,
    height = 2000, width = 2000, resolution = 300,
    compression = "lzw",
    fill = c("#FF6B6B", "#4ECDC4", "#45B7D1"),
    alpha = 0.5,
    cex = 1.2,
    cat.cex = 1.2,
    cat.fontface = "bold",
    fontfamily = "sans",
    main = "Common upgene",
    main.cex = 1.5,
    main.fontface = "bold"
  )
  
  png(file.path(outdir, "venn_upregulated.png"), width = 800, height = 600)
  grid.draw(venn_up)
  dev.off()
  
  # Down-regulated gene Venn diagram
  cat("Drawing down-regulated Venn diagram...\n")
  down_gene_sets <- list(
    Protein = protein_down_genes,
    RNA = rna_down_genes,
    Methylation = methylation_up_genes
  )
  
  venn_down <- venn.diagram(
    x = down_gene_sets,
    category.names = c("Protein", "RNA_seq", "Methylation"),
    filename = NULL,
    output = TRUE,
    height = 2000, width = 2000, resolution = 300,
    compression = "lzw",
    fill = c("#E9A6A6", "#98DDCA", "#87A2FF"),
    alpha = 0.5,
    cex = 1.2,
    cat.cex = 1.2,
    cat.fontface = "bold",
    fontfamily = "sans",
    main = "Common downgene",
    main.cex = 1.5,
    main.fontface = "bold"
  )
  
  png(file.path(outdir, "venn_downregulated.png"), width = 800, height = 600)
  grid.draw(venn_down)
  dev.off()
  
  cat("\nAnalysis complete! Results saved to:", outdir, "\n")
  cat("Output files:\n")
  cat("  - venn_three_omics.png\n")
  cat("  - three_omics_intersect_genes.csv\n")
  cat("  - final_26gene.csv\n")
  cat("  - venn_upregulated.png\n")
  cat("  - venn_downregulated.png\n")
  cat("  - analysis.log\n")
  cat("=== Analysis Finished ===", Sys.time(), "\n")
}

# Execute with error handling
tryCatch({
  main()
}, error = function(e) {
  cat("\nError:", e$message, "\n")
  quit(status = 1)
})

sink()