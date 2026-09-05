# ===========================================================================
# Script: 06_Machine_Learning_Analysis.R
# Project: Thymoma Multi-Omics Biomarker Discovery
# Description: End-to-end machine learning pipeline for thymoma biomarkers
#              - ComBat-seq batch correction across GSE177522 + GSE79978
#              - Builds every <feature selection> + <classifier> combination
#                listed in methods.xlsx over the 26 multi-omics candidates
#              - AUC benchmarking on the training set and the TCGA test set
#              - Differential boxplots, per-gene ROC and multi-gene ROC for
#                the best model
#              - Nested CV
# Date: 2026-09-05
# ===========================================================================
# Data Preparation and Batch Effect Correction
rm(list = ls())
setwd("D:/Thymoma_MultiOmics/results/machine_learning/")

# Load expression data from two GEO datasets
data22 <- read.csv("GSE177522counts.csv", row.names = 1)  # GSE177522
data78 <- read.csv("GSE79978counts.csv", row.names = 1)    # GSE79978

# Keep common genes
samegene <- intersect(rownames(data22), rownames(data78))
data22 <- data22[samegene, ]
data78 <- data78[samegene, ]

# Merge datasets
combined_counts <- cbind(data78, data22)

# Define batch and biological groups
batch <- factor(c(rep("data78", ncol(data78)), rep("data22", ncol(data22))))

group_78 <- factor(ifelse(
  colnames(data78) %in% c("GSM2109588", "GSM2109589", "GSM2109590"),
  "Normal", "Tumor"
))

group_22 <- factor(ifelse(
  colnames(data22) %in% c("GSM2109588", "GSM2109589", "GSM2109590",
                          "SRR14794296", "SRR14794297", "SRR14794298",
                          "SRR14794299", "SRR14794301"),
  "Normal", "Tumor"
))

condition <- factor(c(as.character(group_78), as.character(group_22)),
                    levels = c("Normal", "Tumor"))

# Batch effect correction using ComBat-seq
library(limma)
library(sva)
library(ggplot2)
library(RColorBrewer)

stopifnot(all(combined_counts == floor(combined_counts)))  # Check integer counts

combat_seq_result <- ComBat_seq(
  counts = as.matrix(combined_counts),
  batch = batch,
  group = condition,
  full_mod = TRUE,
  shrink = FALSE,
  shrink.disp = FALSE
)

# PCA visualization after correction
log_corrected <- log2(combat_seq_result + 1)
pca_result <- prcomp(t(log_corrected))
pca_data <- as.data.frame(pca_result$x[, 1:2])

batch_colors <- brewer.pal(n = length(levels(batch)), name = "Set1")
condition_shapes <- c(16, 17, 18)[1:length(levels(condition))]

ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = batch, shape = condition), size = 4, alpha = 0.8) +
  scale_color_manual(values = batch_colors, name = "Batch") +
  scale_shape_manual(values = condition_shapes, name = "Condition") +
  labs(title = "PCA Plot After ComBat-seq Correction",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "gray80", fill = NA)
  )

# Load candidate genes (26 multi-omics intersecting genes)
newsamegene <- read.csv("common_gene.csv", row.names = 1)
newsamegene <- newsamegene$x
Train_expr <- combat_seq_result[newsamegene, ]

# Load training labels
Train_class <- read.table("Training_class.txt", header = TRUE)

# Load TCGA test data
Testing_expr <- read.csv("TCGA_expression.csv", row.names = 1)
Testing_expr <- Testing_expr[, -1]
Test_class <- read.table("Testing_class.txt", header = TRUE)
Test_class$Cohort <- "TCGA"
rownames(Test_class) <- Test_class$sample

Testing_expr <- Testing_expr[newsamegene, ]
Testing_expr <- as.matrix(Testing_expr)

# Normalize training and test sets using training set parameters
library(caret)

train_data <- Train_expr
test_data <- Testing_expr

preproc_params <- preProcess(t(train_data), method = c("center", "scale"))
train_data_norm <- predict(preproc_params, t(train_data))
test_data_norm <- predict(preproc_params, t(test_data))

train_data_norm <- t(train_data_norm)
test_data_norm <- t(test_data_norm)

Train_set <- t(train_data_norm)
Test_set <- t(test_data_norm)

identical(colnames(Train_set), colnames(Test_set))

# ============================================================================
# Machine Learning Model Training and Evaluation
# ============================================================================
library(openxlsx)
library(glmnet)
library(plsRglm)
library(gbm)
library(mboost)
library(e1071)
library(MASS)
library(xgboost)
library(ComplexHeatmap)
library(pROC)
library(ggpubr)

source("./05_ML_Algorithm_Wrappers.R")

rownames(Train_class) <- Train_class$sample
Train_class <- read.table("Training_class.txt", header = TRUE)

# Load model list
methods <- read.xlsx("methods.xlsx", startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)

classVar <- "outcome"
min.selected.var <- 3

# ============================================================================
# Pre-training: Feature selection on full training set
# ============================================================================

preTrain.var <- list()
set.seed(777)

Variable <- colnames(Train_set)
preTrain.method <- strsplit(methods, "\\+")
preTrain.method <- lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method <- unique(unlist(preTrain.method))

for (method in preTrain.method) {
  preTrain.var[[method]] <- RunML(
    method = method,
    Train_set = Train_set,
    Train_label = Train_class,
    mode = "Variable",
    classVar = classVar
  )
}

preTrain.var[["simple"]] <- colnames(Train_set)

# ============================================================================
# Model training
# ============================================================================
model <- list()
set.seed(777)
Train_set_bk <- Train_set

for (method in methods) {
  cat(match(method, methods), ":", method, "\n")
  method_name <- method
  method_split <- strsplit(method, "\\+")[[1]]
  if (length(method_split) == 1) method_split <- c("simple", method_split)
  
  Variable <- preTrain.var[[method_split[1]]]
  Train_set <- Train_set_bk[, Variable, drop = FALSE]
  Train_label <- Train_class
  
  # Skip XGBoost due to compatibility issues
  if (grepl("XGBoost", method_name, ignore.case = TRUE)) {
    cat("  Skipping XGBoost\n")
    next
  }
  
  tryCatch({
    model[[method_name]] <- RunML(
      method = method_split[2],
      Train_set = Train_set,
      Train_label = Train_label,
      mode = "Model",
      classVar = classVar
    )
    
    if (!is.null(model[[method_name]]) &&
        length(ExtractVar(model[[method_name]])) <= min.selected.var) {
      model[[method_name]] <- NULL
      cat("  Removed: too few variables\n")
    } else {
      cat("  Model trained successfully\n")
    }
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
  })
}

Train_set <- Train_set_bk
methodsValid <- names(model)

# ============================================================================
# Calculate risk scores and class predictions
# ============================================================================

RS_list <- list()
for (method in methodsValid) {
  RS_list[[method]] <- CalPredictScore(
    fit = model[[method]],
    new_data = rbind.data.frame(Train_set, Test_set)
  )
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))

Class_list <- list()
for (method in methodsValid) {
  Class_list[[method]] <- PredictClass(
    fit = model[[method]],
    new_data = rbind.data.frame(Train_set, Test_set)
  )
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))


# ============================================================================
# Extract feature genes from each model
# ============================================================================

fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}

fea_df <- lapply(model, function(fit) {
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"

# ============================================================================
# Calculate AUC for each model 
# ============================================================================
AUC_list <- list()
for (method in methodsValid) {
  AUC_list[[method]] <- RunEval(
    fit = model[[method]],
    Test_set = Test_set,
    Test_label = Test_class,
    Train_name = "Training",
    Train_set = Train_set,
    Train_label = Train_class,
    cohortVar = "Cohort",
    classVar = classVar
  )
}
AUC_mat <- do.call(rbind, AUC_list)
write.table(AUC_mat,file = "26gene_AUC_mat.txt")
# ============================================================================
# Generate heatmap 
# ============================================================================

AUC_mat <- read.table("26gene_AUC_mat.txt", sep = "\t", row.names = 1,
                      header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Sort by TCGA column in descending order
AUC_mat <- AUC_mat[order(AUC_mat$TCGA, decreasing = TRUE), ]

# Define cohort colors
if (ncol(AUC_mat) < 3) {
  CohortCol <- c("#F55A5A", "#65B2FA")
} else {
  CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired")
}
names(CohortCol) <- colnames(AUC_mat)

cellwidth <- 1
cellheight <- 0.5

# Custom heatmap function without average bar plot
SimpleHeatmap2 <- function(Cindex_mat, CohortCol,
                           heatmap_col = c("#4195C1", "#FFFFFF", "#CB5746"),
                           cellwidth = 1, cellheight = 0.5,
                           cluster_columns, cluster_rows) {
  
  col_ha <- columnAnnotation(
    "Cohort" = colnames(Cindex_mat),
    col = list("Cohort" = CohortCol),
    show_annotation_name = FALSE
  )
  
  Heatmap(
    as.matrix(Cindex_mat),
    name = "AUC",
    top_annotation = col_ha,
    col = heatmap_col,
    rect_gp = gpar(col = "black", lwd = 1),
    cluster_columns = cluster_columns,
    cluster_rows = cluster_rows,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
    height = unit(cellheight * nrow(Cindex_mat), "cm"),
    column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)),
    column_title = NULL,
    cell_fun = function(j, i, x, y, w, h, col) {
      grid.text(
        label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
        x, y, gp = gpar(fontsize = 10)
      )
    }
  )
}

hm <- SimpleHeatmap2(
  AUC_mat,
  CohortCol,
  heatmap_col = c("#65B2FA", "white", "#F55A5A"),
  cellwidth = cellwidth,
  cellheight = cellheight,
  cluster_columns = FALSE,
  cluster_rows = FALSE
)

pdf("new_color_26_new_AUC.pdf",
    width = cellwidth * ncol(AUC_mat) + 7,
    height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm)
dev.off()

# ============================================================================
# Define the optimal model
# ============================================================================

best.model <- "glmBoost+Enet[alpha=0.1]"
best.model.gene <- subset(fea_df, algorithm == best.model)[, 1]

# ============================================================================
# Differential expression analysis for the selected genes
# ============================================================================

# Training set
Train <- t(Train_expr)
Train <- t(Train[, best.model.gene])

con <- subset(Train_class, outcome == 0)
treat <- subset(Train_class, outcome == 1)
rownames(con) <- con$sample
rownames(treat) <- treat$sample

conData <- Train[, rownames(con)]
treatData <- Train[, rownames(treat)]
Train <- cbind(conData, treatData)

Type <- c(rep("Normal", ncol(conData)), rep("Disease", ncol(treatData)))
my_comparisons <- list()
my_comparisons[[1]] <- levels(factor(Type))

dir.create("diff", showWarnings = FALSE)

for (i in rownames(Train)) {
  rt1 <- data.frame(expression = Train[i, ], Type = Type)
  boxplot <- ggboxplot(rt1, x = "Type", y = "expression", color = "Type",
                       xlab = "", ylab = paste(i, "expression"),
                       legend.title = "",
                       palette = c("#00AF50", "#F5B700"),
                       add = "jitter") +
    stat_compare_means(comparisons = my_comparisons, method = "t.test")
  pdf(file = paste0("diff/Training.diff.", i, ".pdf"), width = 5, height = 4.5)
  print(boxplot)
  dev.off()
}

# TCGA test set
Test_expr <- read.csv("merged_TCGA_expression.csv", row.names = 1)
Test_expr <- Test_expr[, -1]
Test <- t(Test_expr)
Test <- t(Test[, best.model.gene])

rownames(Test_class) <- Test_class$sample
con <- subset(Test_class, outcome == 0 & Cohort == "TCGA")[, 2, drop = FALSE]
treat <- subset(Test_class, outcome == 1 & Cohort == "TCGA")[, 2, drop = FALSE]

conData <- Test[, rownames(con)]
treatData <- Test[, rownames(treat)]
Test <- cbind(conData, treatData)

Type <- c(rep("Normal", ncol(conData)), rep("Disease", ncol(treatData)))
my_comparisons <- list()
my_comparisons[[1]] <- levels(factor(Type))

for (i in rownames(Test)) {
  rt1 <- data.frame(expression = Test[i, ], Type = Type)
  boxplot <- ggboxplot(rt1, x = "Type", y = "expression", color = "Type",
                       xlab = "", ylab = paste(i, "expression"),
                       legend.title = "",
                       palette = c("#00AF50", "#F5B700"),
                       add = "jitter") +
    stat_compare_means(comparisons = my_comparisons, method = "t.test")
  pdf(file = paste0("diff/TCGA.diff.", i, ".pdf"), width = 5, height = 4.5)
  print(boxplot)
  dev.off()
}

# ============================================================================
# Individual gene ROC analysis
# ============================================================================

dir.create("ROC", showWarnings = FALSE)

# Training set
Train <- Train_set[, best.model.gene]
rownames(Train_class) <- Train_class$sample

con <- Train_class[Train_class$outcome == 0, "sample"]
treat <- Train_class[Train_class$outcome == 1, "sample"]

Train_sub <- Train[c(con, treat), ]
y_train <- c(rep(0, length(con)), rep(1, length(treat)))

for (gene in best.model.gene) {
  roc_result <- roc(y_train, Train_sub[, gene], quiet = TRUE)
  auc_value <- auc(roc_result)
  ci_result <- ci.auc(roc_result, method = "bootstrap")
  
  pdf(file = paste0("ROC/Training_", gene, "_ROC.pdf"), width = 5, height = 5)
  plot(roc_result, print.auc = TRUE, auc.polygon = TRUE, grid = TRUE,
       print.thres = "best", print.thres.best.method = "youden",
       col = "red", legacy.axes = TRUE,
       main = paste("Training -", gene),
       xlab = "1 - Specificity", ylab = "Sensitivity")
  dev.off()
}

# TCGA test set
Test <- t(Test_set[, best.model.gene])

con_samples <- Test_class$sample[Test_class$outcome == 0 & Test_class$Cohort == "TCGA"]
treat_samples <- Test_class$sample[Test_class$outcome == 1 & Test_class$Cohort == "TCGA"]

con_samples <- con_samples[con_samples %in% colnames(Test)]
treat_samples <- treat_samples[treat_samples %in% colnames(Test)]

Test_sub <- Test[, c(con_samples, treat_samples)]
y_test <- c(rep(0, length(con_samples)), rep(1, length(treat_samples)))

for (gene in rownames(Test_sub)) {
  if (gene %in% best.model.gene) {
    roc_obj <- roc(y_test, Test_sub[gene, ], quiet = TRUE)
    auc_value <- auc(roc_obj)
    ci_result <- ci.auc(roc_obj, method = "bootstrap")
    
    pdf(paste0("ROC/TCGA_", gene, "_ROC.pdf"), width = 5, height = 5)
    plot(roc_obj, print.auc = TRUE, auc.polygon = TRUE, grid = TRUE,
         print.thres = "best", print.thres.best.method = "youden",
         col = "blue", legacy.axes = TRUE,
         main = paste("TCGA -", gene),
         xlab = "1 - Specificity", ylab = "Sensitivity")
    dev.off()
  }
}

# ============================================================================
# Generate AUC summary table
# ============================================================================

auc_summary <- data.frame()

# Training set
for (gene in best.model.gene) {
  roc_obj <- roc(y_train, Train_sub[, gene], quiet = TRUE)
  ci_result <- ci.auc(roc_obj, method = "bootstrap")
  
  auc_summary <- rbind(auc_summary, data.frame(
    Gene = gene,
    Dataset = "Training",
    AUC = auc(roc_obj),
    AUC_CI = paste0(sprintf("%.3f", ci_result[1]), "-", sprintf("%.3f", ci_result[3])),
    Sensitivity = sprintf("%.1f%%", roc_obj$sensitivities[which.max(roc_obj$sensitivities + roc_obj$specificities - 1)] * 100),
    Specificity = sprintf("%.1f%%", roc_obj$specificities[which.max(roc_obj$sensitivities + roc_obj$specificities - 1)] * 100),
    stringsAsFactors = FALSE
  ))
}

# TCGA test set
for (gene in rownames(Test_sub)) {
  roc_obj <- roc(y_test, Test_sub[gene, ], quiet = TRUE)
  ci_result <- ci.auc(roc_obj, method = "bootstrap")
  
  auc_summary <- rbind(auc_summary, data.frame(
    Gene = gene,
    Dataset = "TCGA",
    AUC = auc(roc_obj),
    AUC_CI = paste0(sprintf("%.3f", ci_result[1]), "-", sprintf("%.3f", ci_result[3])),
    Sensitivity = sprintf("%.1f%%", roc_obj$sensitivities[which.max(roc_obj$sensitivities + roc_obj$specificities - 1)] * 100),
    Specificity = sprintf("%.1f%%", roc_obj$specificities[which.max(roc_obj$sensitivities + roc_obj$specificities - 1)] * 100),
    stringsAsFactors = FALSE
  ))
}

write.csv(auc_summary, "AUC_summary_long.csv", row.names = FALSE)

# ============================================================================
# Multi-gene ROC plot
# ============================================================================

plot_multi_gene_roc <- function(data_set, class_set, dataset_name, filename, colors) {
  rownames(class_set) <- class_set$sample
  con <- subset(class_set, outcome == 0)
  treat <- subset(class_set, outcome == 1)
  
  conData <- data_set[rownames(con), , drop = FALSE]
  treatData <- data_set[rownames(treat), , drop = FALSE]
  combined <- rbind(conData, treatData)
  
  y <- c(rep(0, nrow(conData)), rep(1, nrow(treatData)))
  
  roc_list <- list()
  auc_values <- c()
  
  for (gene in best.model.gene) {
    if (!gene %in% colnames(combined)) next
    expr <- combined[, gene]
    if (any(is.na(expr))) next
    roc_obj <- roc(y, expr, quiet = TRUE)
    roc_list[[gene]] <- roc_obj
    auc_values[gene] <- auc(roc_obj)
  }
  
  gene_order <- names(sort(auc_values, decreasing = TRUE))
  
  pdf(filename, width = 8, height = 8)
  par(mar = c(5, 4, 4, 6) + 0.1)
  
  first_gene <- gene_order[1]
  plot(roc_list[[first_gene]], col = colors[first_gene], lwd = 2,
       legacy.axes = TRUE,
       main = paste("ROC Curves -", dataset_name),
       xlab = "1 - Specificity", ylab = "Sensitivity",
       grid = TRUE, print.auc = FALSE, print.thres = FALSE)
  
  if (length(gene_order) > 1) {
    for (i in 2:length(gene_order)) {
      gene <- gene_order[i]
      lines(roc_list[[gene]], col = colors[gene], lwd = 2)
    }
  }
  
  legend_text <- sapply(gene_order, function(g) {
    sprintf("%s (AUC=%.3f)", g, auc_values[g])
  })
  
  legend("bottomright", legend = legend_text,
         col = colors[gene_order], lty = 1, lwd = 2,
         bty = "n", cex = 0.8, title = "Genes")
  dev.off()
  
  return(list(roc_list = roc_list, auc_values = auc_values))
}

# Set colors
n_genes <- length(best.model.gene)
if (n_genes <= 8) {
  colors <- brewer.pal(n_genes, "Set1")
} else {
  colors <- rainbow(n_genes, alpha = 0.8)
}
names(colors) <- best.model.gene

# Generate multi-gene ROC plots
train_result <- plot_multi_gene_roc(
  data_set = Train_set[, best.model.gene],
  class_set = Train_class,
  dataset_name = "Training Set",
  filename = "ROC/Training_MultiGene_ROC.pdf",
  colors = colors
)

test_result <- plot_multi_gene_roc(
  data_set = Test_set[, best.model.gene],
  class_set = Test_class,
  dataset_name = "TCGA Test Set",
  filename = "ROC/TCGA_MultiGene_ROC.pdf",
  colors = colors
)

# ==================== Nested Cross-Validation  ====================
library(glmnet)
library(caret)
library(pROC)
library(ggplot2)

set.seed(123)
outer_folds <- createFolds(y, k = 5, list = TRUE)

auc_outer <- c()
all_pred <- c()
all_true <- c()

alpha_candidates <- c(0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 1)

for (i in seq_along(outer_folds)) {
  outer_test_idx <- outer_folds[[i]]
  outer_train_idx <- setdiff(1:length(y), outer_test_idx)
  
  train_x_outer <- Train_data[outer_train_idx, , drop = FALSE]
  train_y_outer <- y[outer_train_idx]
  test_x_outer  <- Train_data[outer_test_idx, , drop = FALSE]
  test_y_outer  <- y[outer_test_idx]
  
  nfolds_outer <- min(5, length(train_y_outer))
  if (nfolds_outer < 4 || length(unique(train_y_outer)) < 2) {
    cat("Fold", i, ": insufficient outer training samples or classes, skipped\n")
    next
  }
  
  # ------ Inner 5-fold: select alpha ------
  inner_folds <- createFolds(train_y_outer, k = 5, list = TRUE)
  best_alpha <- 0.1
  best_inner_auc <- 0
  
  for (alpha_val in alpha_candidates) {
    inner_aucs <- c()
    
    for (j in 1:5) {
      inner_test_idx <- inner_folds[[j]]
      inner_train_idx <- setdiff(1:length(train_y_outer), inner_test_idx)
      
      train_x_inner <- train_x_outer[inner_train_idx, , drop = FALSE]
      train_y_inner <- train_y_outer[inner_train_idx]
      test_x_inner  <- train_x_outer[inner_test_idx, , drop = FALSE]
      test_y_inner  <- train_y_outer[inner_test_idx]
      
      if (length(unique(train_y_inner)) < 2) {
        inner_aucs <- c(inner_aucs, NA)
        next
      }
      
      nfolds_inner <- min(5, length(train_y_inner))
      if (nfolds_inner < 4) {
        inner_aucs <- c(inner_aucs, NA)
        next
      }
      
      cv_fit_inner <- tryCatch(
        cv.glmnet(
          x = train_x_inner,
          y = train_y_inner,
          family = "binomial",
          alpha = alpha_val,
          nfolds = nfolds_inner
        ),
        error = function(e) NULL
      )
      
      if (is.null(cv_fit_inner)) {
        inner_aucs <- c(inner_aucs, NA)
        next
      }
      
      pred_inner <- as.vector(
        predict(cv_fit_inner,
                newx = test_x_inner,
                type = "response",
                s = "lambda.min")
      )
      
      if (!any(is.na(pred_inner)) && length(unique(test_y_inner)) == 2) {
        roc_obj <- tryCatch(roc(test_y_inner, pred_inner, quiet = TRUE),
                            error = function(e) NULL)
        if (!is.null(roc_obj)) {
          inner_aucs <- c(inner_aucs, auc(roc_obj))
        } else {
          inner_aucs <- c(inner_aucs, NA)
        }
      } else {
        inner_aucs <- c(inner_aucs, NA)
      }
    }
    
    mean_inner <- mean(inner_aucs, na.rm = TRUE)
    if (!is.na(mean_inner) && mean_inner > best_inner_auc) {
      best_inner_auc <- mean_inner
      best_alpha <- alpha_val
    }
  }
  
  cat("Fold", i, "best alpha:", best_alpha,
      "(inner mean AUC =", round(best_inner_auc, 4), ")\n")
  
  # ------ Outer evaluation: fit with selected alpha ------
  cv_fit_outer <- tryCatch(
    cv.glmnet(
      x = train_x_outer,
      y = train_y_outer,
      family = "binomial",
      alpha = best_alpha,
      nfolds = nfolds_outer
    ),
    error = function(e) NULL
  )
  
  if (is.null(cv_fit_outer)) {
    cat("Fold", i, ": outer training failed\n")
    next
  }
  
  pred_outer <- as.vector(
    predict(cv_fit_outer,
            newx = test_x_outer,
            type = "response",
            s = "lambda.min")
  )
  
  if (!any(is.na(pred_outer)) && length(unique(test_y_outer)) == 2) {
    roc_obj <- tryCatch(roc(test_y_outer, pred_outer, quiet = TRUE),
                        error = function(e) NULL)
    if (!is.null(roc_obj)) {
      outer_auc <- auc(roc_obj)
      auc_outer <- c(auc_outer, outer_auc)
      all_pred <- c(all_pred, pred_outer)
      all_true <- c(all_true, test_y_outer)
      cat("  Outer AUC =", round(outer_auc, 4), "\n")
    } else {
      cat("Fold", i, ": outer ROC computation failed\n")
    }
  } else {
    cat("Fold", i, ": outer evaluation failed / insufficient classes\n")
  }
}

# ==================== Nested CV Summary ====================
cat("\n========== Nested CV Results ==========\n")
cat("Valid outer folds:", length(auc_outer), "/ 5\n")

if (length(auc_outer) > 0) {
  cat("Nested CV AUC mean:", round(mean(auc_outer), 4), "\n")
  cat("Nested CV AUC SD:", round(sd(auc_outer), 4), "\n")
  cat("Nested CV AUC 95% CI:",
      round(quantile(auc_outer, 0.025), 4), "-",
      round(quantile(auc_outer, 0.975), 4), "\n")
} else {
  cat("No valid nested CV AUC values\n")
}

# ==================== Pooled Outer Predictions: AUC + Bootstrap CI + ROC ====================
if (length(all_true) > 0 && length(unique(all_true)) == 2) {
  
  roc_all <- roc(all_true, all_pred, quiet = TRUE)
  auc_all <- auc(roc_all)
  
  set.seed(777)
  n_boot <- 1000
  auc_boot <- c()
  
  for (b in 1:n_boot) {
    idx <- sample(1:length(all_true), size = length(all_true), replace = TRUE)
    boot_true <- all_true[idx]
    boot_pred <- all_pred[idx]
    
    if (length(unique(boot_true)) == 2) {
      roc_boot <- tryCatch(roc(boot_true, boot_pred, quiet = TRUE),
                           error = function(e) NULL)
      if (!is.null(roc_boot)) {
        auc_boot <- c(auc_boot, auc(roc_boot))
      }
    }
  }
  
  ci_boot <- quantile(auc_boot, c(0.025, 0.975), na.rm = TRUE)
  
  cat("\n========== Pooled AUC (all outer folds) & Bootstrap ==========\n")
  cat("Pooled AUC:", round(auc_all, 4), "\n")
  cat("Bootstrap AUC mean:", round(mean(auc_boot), 4), "\n")
  cat("Bootstrap AUC SD:", round(sd(auc_boot), 4), "\n")
  cat("Bootstrap 95% CI:", round(ci_boot[1], 4), "-", round(ci_boot[2], 4), "\n")
  
  pdf("Nested_CV_ROC.pdf", width = 6, height = 6)
  plot(roc_all,
       main = paste("Nested CV ROC (AUC =", round(auc_all, 3), ")"),
       col = "#D95F02",
       lwd = 3,
       legacy.axes = TRUE,
       xlab = "1 - Specificity",
       ylab = "Sensitivity")
  abline(a = 0, b = 1, col = "gray50", lty = 2, lwd = 1.5)
  legend("bottomright",
         legend = c(paste("AUC =", round(auc_all, 3)),
                    paste("95% CI:", round(ci_boot[1], 3), "-", round(ci_boot[2], 3))),
         col = c("#D95F02", NA),
         lty = c(1, NA),
         lwd = c(3, NA),
         bty = "n")
  dev.off()
  
  cat("Nested CV ROC saved as Nested_CV_ROC.pdf\n")
}

# ==================== TCGA  ====================
Test_data <- as.matrix(Test_set[, best.model.gene])
pred_prob <- CalPredictScore(fit = model[[best.model]], new_data = Test_data)

Test_y <- Test_class[rownames(Test_data), "outcome"]

cat("TCGA rowname match:",
    identical(rownames(Test_class), rownames(Test_data)), "\n")
cat("TCGA class distribution:\n")
print(table(Test_y))

if (length(unique(Test_y)) == 2) {
  tcga_roc <- roc(Test_y, pred_prob, quiet = TRUE)
  tcga_auc <- auc(tcga_roc)
  tcga_ci  <- ci.auc(tcga_roc, method = "bootstrap")
  
  cat("\n========== TCGA Results ==========\n")
  cat("Samples:", length(Test_y), "\n")
  cat("Normal:", sum(Test_y == 0), ", Tumor:", sum(Test_y == 1), "\n")
  cat("AUC:", round(tcga_auc, 4), "\n")
  cat("95% CI:", round(tcga_ci[1], 4), "-", round(tcga_ci[2], 4), "\n")
  
  pdf("TCGA_ROC.pdf", width = 6, height = 6)
  plot(tcga_roc,
       main = paste("TCGA  (AUC =", round(tcga_auc, 3), ")"),
       col = "#2E9FDF",
       lwd = 3,
       legacy.axes = TRUE,
       xlab = "1 - Specificity",
       ylab = "Sensitivity")
  abline(a = 0, b = 1, col = "gray50", lty = 2, lwd = 1.5)
  legend("bottomright",
         legend = c(paste("AUC =", round(tcga_auc, 3)),
                    paste("95% CI:", round(tcga_ci[1], 3), "-", round(tcga_ci[2], 3))),
         col = c("#2E9FDF", NA),
         lty = c(1, NA),
         lwd = c(3, NA),
         bty = "n")
  dev.off()
  
  cat("ROC saved as TCGA_ROC.pdf\n")
} else {
  cat("TCGA validation: insufficient classes for AUC/ROC.\n")
}                          
