# Results

Final result tables produced by the pipeline. All files here are tracked in
version control; intermediate artefacts (`diff/`, `ROC/`, plots `.rds`) are
excluded — see [`../.gitignore`](../.gitignore).

## Overview

| File                            | Rows   | Produced by              | Description                            |
| ------------------------------- | ------ | ------------------------ | -------------------------------------- |
| `RNAseq_DEGs.csv`               | 1,833  | `01` (exported manually) | Transcriptome DEGs, limma-voom         |
| `proteomics_DEGs.csv`           | 1,501  | `03` (exported manually) | Proteome DEPs, limma                   |
| `DNA_methylation_DMPs_gene.csv` | 7,181  | `02`                     | Differentially methylated genes, ChAMP |
| `common_gene.csv`               | 26     | `04`                     | Direction-aware three-omics signature  |
| `AUC_mat..xlsx`                 | 56 × 2 | `06`                     | Per-model AUC (Training, TCGA)         |
| `AUC_summary_long.csv`          | 12     | `06`                     | Per-gene AUC for the 6 feature genes   |
| `scTenifoldKnk_results.csv`     | 86     | `07`                     | PKP1 virtual-knockout regulation table |

## Differential-expression tables

`RNAseq_DEGs.csv`, `proteomics_DEGs.csv` and `DNA_methylation_DMPs_gene.csv`
share the limma `topTable` schema:

`logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`, plus a gene-symbol column
(`Gene` for RNA-seq, `gene` for proteomics and methylation).

The methylation table additionally carries `Tumor_AVG`, `Normal_AVG`,
`deltaBeta` and probe annotation (`CHR`, `MAPINFO`, `Type`, `feature`, `cgi`,
…).

Thresholds used throughout: **|log2FC| > 1** and **adj.P.Val < 0.05** for
RNA and protein; **|deltaBeta| > 0.2** and **adj.P.Val < 0.05** for methylation.



## `common_gene.csv` — the 26-gene signature

Column `x` holds the gene symbols consumed by script `06`
(`newsamegene <- read.csv("common_gene.csv"); newsamegene$x`).

Direction-aware composition: **25** genes from the up-set (RNA up + protein up

+ hypomethylated) and **1** gene from the down-set (RNA down + protein down +
  hypermethylated).

```
IRF6   GGA2    S100A14 TRIP6   SH3PXD2B PKP1   HSPA4L  KIAA1217 FBLIM1
S100A16 MYO1E  SQSTM1  CTSA    HGS      PKP3   CDH1    PFN2     ENAH
PBX1   EPN3    LTBP2   HRAS    EPS8L2   AGRN   JUP     KRT1
```

## `AUC_mat.xlsx` — model benchmark

56 model combinations × 2 cohorts (`ID`, `Training`, `TCGA`), sorted by TCGA
AUC. Best model: **`glmBoost+Enet[alpha=0.1]`**, training AUC 0.983,
validation AUC 0.983.

The file name has a double dot (`AUC_mat..xlsx`). It can be renamed to
`AUC_mat.xlsx` without affecting any script — nothing reads it back.

## `AUC_summary_long.csv` — per-gene diagnostic performance

Six feature genes × 2 cohorts. These are the numbers behind Figure 3C/D and
Table S3.

| Gene | Training AUC | 95% CI | Sens | Spec | TCGA AUC | 95% CI | Sens | Spec |
|------|--------------|--------|------|------|----------|--------|------|
| IRF6   | 0.775 | 0.558–0.942 | 66.7% | 100% | 0.796 | 0.725–0.867 | 79.2% | 100% |
| TRIP6  | 0.800 | 0.608–0.967 | 80.0% | 87.5% | 0.900 | 0.825–0.967 | 86.7% | 100% |
| **PKP1** | **0.833** | 0.633–0.983 | 80.0% | 87.5% | **0.942** | 0.892–0.983 | 93.3% | 100% |
| CDH1   | 0.875 | 0.692–1.000 | 80.0% | 100% | 0.742 | 0.625–0.858 | 67.5% | 100% |
| PBX1   | 0.817 | 0.617–0.983 | 73.3% | 100% | 0.546 | 0.292–0.800 | 34.2% | 100% |
| KRT1   | 0.983 | 0.925–1.000 | 100%  | 87.5% | 0.958 | 0.883–1.000 | 91.7% | 100% |

## `scTenifoldKnk_results.csv` — PKP1 virtual knockout

86 rows: the knocked-out gene `PKP1` itself plus **85 perturbation-responsive
genes** (adj. p < 0.05), matching the manuscript's count of 85.

Columns: `gene`, `distance`, `Z`, `FC`, `p.value`, `p.adj`, `log2FC`.

