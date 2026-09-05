# Thymoma Multi-Omics Biomarker Discovery

Reproducible analysis code for:

> **Integrative Multi-Omics and Machine Learning Unveil PKP1 as a Robust Diagnostic Biomarker for Thymoma**
> Yunhan Li¹'³†, Shen Zhang¹'³†, Qian Gao¹'³, Shanshan Ding¹'³, Mengling Sun¹'³, Haiying Xue¹'², Fan Yang¹'², Yangyang Liu¹'³, Xiaokang An¹'², Xiaotao Dong¹'²'³\*, Guoyu Zhang¹'²\*
> ¹ Department of Thoracic Surgery, The First Affiliated Hospital, Henan University, Kaifeng 475004, China
> ² Kaifeng Key Laboratory of Translational Medicine for Thoracic Diseases, Kaifeng 475004, China
> ³ Laboratory of Receptor Gene Regulation and Drug Discovery, School of Basic Medical Sciences, Henan University, Kaifeng 475004, China

This integrated multi-omics pipeline consolidates bulk transcriptomic, DNA methylation, and proteomic profiling to distill approximately 20,000 genes into a robust **26-gene signature**. We systematically benchmarked **56 machine learning algorithmic combinations** on this signature, and subsequently validated the leading candidate, **PKP1**, through single-cell resolution mapping and immunohistochemical staining

---

## Headline results

| Step                                     | Result                                                       |
| ---------------------------------------- | ------------------------------------------------------------ |
| Transcriptome DEGs                       | **1,833** (941 up / 892 down)                                |
| Proteome DEPs                            | **1,501** (1,122 up / 379 down)                              |
| Differentially methylated genes          | **7,181** (1,354 hyper / 5,827 hypo)                         |
| Direction-aware three-omics intersection | **26 genes** (25 up-set + 1 down-set)                        |
| ML models benchmarked                    | **56**; best = `glmBoost + Enet[alpha=0.1]`                  |
| Best model performance                   | Training AUC **0.983** (95% CI 0.944–1.000); TCGA validation AUC **0.983** (0.947–1.000) |
| Core feature genes                       | **IRF6, TRIP6, PKP1, CDH1, PBX1, KRT1**                      |
| Lead biomarker PKP1                      | Training AUC **0.833** (0.633–0.983); validation AUC **0.942** (0.892–0.983) |
| Single-cell                              | 33,725 cells → 6 cell types                                  |
| scTenifoldKnk virtual knockout           | **85** perturbation-responsive genes (adj. p < 0.05)         |
| IHC validation                           | 25 thymoma vs 8 adjacent normal, ***P* < 0.001**             |

---

## Workflow

```
        ┌───────────────────────┐  ┌───────────────────────┐  ┌───────────────────────┐
        │ 01 Transcriptome      │  │ 02 Methylation        │  │ 03 Proteome           │
        │ GSE79978 + GSE177522  │  │ GSE218549 (EPIC)      │  │ PXD039744             │
        │ ComBat-seq + voom     │  │ ChAMP DMP             │  │ limma + GO/KEGG       │
        │ 1,833 DEGs            │  │ 7,181 DMGs            │  │ 1,501 DEPs            │
        └───────────┬───────────┘  └───────────┬───────────┘  └───────────┬───────────┘
                    └──────────────────────────┼──────────────────────────┘
                                               v
                    ┌──────────────────────────────────────────────────┐
                    │ 04 Direction-aware three-omics intersection      │
                    │    up-set   : RNA up   + protein up   + hypo-    │
                    │    down-set : RNA down + protein down + hyper-   │
                    │    -> 26 candidate genes                          │
                    └───────────────────────┬──────────────────────────┘
                                            v
              ┌─────────────────────────────────────────────────────┐
              │ 05 + 06 Machine learning                            │
              │ 56 model combinations, TCGA as independent test set  │
              │ -> 6 feature genes, PKP1 ranked best                 │
              └──────────────────────┬──────────────────────────────┘
                                     v
              ┌─────────────────────────────────────────────────────┐
              │ 07 Single-cell validation (SCP1532)                 │
              │ 6 cell types + scTenifoldKnk PKP1 knockout          │
              └─────────────────────────────────────────────────────┘
```

The design decision that matters is in **script 04**: a plain three-way
intersection is not enough. A gene survives only if its expression change is
*compatible with an epigenetic driver* — DNA methylation generally represses
expression, so an up-regulated gene must be **hypo**methylated and a
down-regulated gene must be **hyper**methylated.

---

## Datasets

| Layer         | Source             | Accession     | Tumor                               | Normal             |
| ------------- | ------------------ | ------------- | ----------------------------------- | ------------------ |
| Transcriptome | GEO                | **GSE79978**  | 13                                  | 3                  |
| Transcriptome | GEO                | **GSE177522** | 2                                   | 5                  |
| Transcriptome | TCGA               | **TCGA-THYM** | 120                                 | 2                  |
| Methylation   | GEO                | **GSE218549** | 113 → **28** analysed               | 9 → **8** analysed |
| Proteome      | ProteomeXchange    | **PXD039744** | 15                                  | 11                 |
| Single cell   | Single Cell Portal | **SCP1532**   | 4 MG-thymoma patients, 33,725 cells | —                  |

GSE79978 + GSE177522 are merged and batch-corrected into a single **training
set of 23 samples (15 tumour / 8 normal)**; TCGA provides an **independent
validation set of 122 samples (120 tumour / 2 normal)**.

Methylation sample counts drop from 122 to 36 because script `02` removes
outliers — tumour samples falling inside the normal 95% PCA ellipsoid, plus 28
manually flagged arrays.


---

## Repository structure

```
Thymoma_MultiOmics/
├── R/                                   Analysis scripts, numbered in run order
│   ├── 01_RNA_seq_Analysis.R            Transcriptome DE + visualisation
│   ├── 02_DNA_Methylation_Analysis.R    ChAMP normalisation, PCA outlier removal, DMP
│   ├── 03_Proteomics_Analysis.R         Proteome DE + GO/KEGG enrichment
│   ├── 04_MultiOmics_Intersection.R     Direction-aware intersection -> 26 genes
│   ├── 05_ML_Algorithm_Wrappers.R       Function library (sourced by 06)
│   ├── 06_Machine_Learning_Analysis.R   56-model benchmark + TCGA validation
│   └── 07_scRNA_seq_Analysis.R          Seurat QC/clustering + scTenifoldKnk
├── data/                                Input data (see data/README.md)
├── results/                             Final result tables (see results/README.md)
├── docs/
│   └── DEPENDENCIES.md                  Package list + install commands
├── .gitignore
├── .gitattributes
├── LICENSE                              MIT
└── README.md
```

Scripts are numbered in execution order. `01`–`03` are mutually independent;
`04` needs all three outputs; `05` is a function library sourced by `06`.

---

## Requirements

- **R 4.5.0** with **Bioconductor 3.21** for scripts `01`–`06`
- **R 4.4.2** for script `07` (developed on a Linux HPC via Paratera; Seurat v4.4.0)

Key package versions used in the study: `sva` 3.58.0, `limma` 3.64.1, `ChAMP`
2.38.0, `clusterProfiler` 4.18.2, `DESeq2` 1.50.2, `ggplot2` 3.5.2, `Seurat`
4.4.0, `scTenifoldKnk` 1.0.3.

See [`docs/DEPENDENCIES.md`](docs/DEPENDENCIES.md) for the full list and
install commands.

---

## Running the analysis

Each script opens with a hardcoded `setwd()` pointing at
`D:/Thymoma_MultiOmics/...`. **Edit these to match your machine** — and note
that several currently point at sub-directories that do not match where the
data actually lives (see [Known gaps](#known-gaps)).

```r
source("R/01_RNA_seq_Analysis.R")          # -> transcriptome DEGs
source("R/02_DNA_Methylation_Analysis.R")  # -> methylation DMPs
source("R/03_Proteomics_Analysis.R")       # -> proteome DEPs + enrichment
source("R/04_MultiOmics_Intersection.R")   # -> 26 candidate genes
source("R/06_Machine_Learning_Analysis.R") # sources 05 internally, then benchmarks
source("R/07_scRNA_seq_Analysis.R")        # single-cell validation (needs HPC-scale RAM)
```

Script `02` expects an EPIC beta matrix as `myLoad$beta` — the raw IDAT files
are not distributed, so section 1 (ChAMP normalisation) must be run against
your own copy of GSE218549.

---

## Input and output files

Detailed per-file documentation lives in [`data/README.md`](data/README.md)
and [`results/README.md`](results/README.md).

**Inputs** (`data/`): `GSE177522counts.csv`, `GSE79978counts.csv`,
`TCGA_expression.csv`, `Training_class.txt`, `Testing_class.txt`,
`proteomics.xlsx`, `methods.xlsx`.

**Outputs** (`results/`): `RNAseq_DEGs.csv`, `proteomics_DEGs.csv`,
`DNA_methylation_DMPs_gene.csv`, `common_gene.csv`, `AUC_mat..xlsx`,
`AUC_summary_long.csv`, `scTenifoldKnk_results.csv`.

---

## Data availability

- **GSE79978**, **GSE177522**, **GSE218549** — [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/)
- **TCGA-THYM** — [GDC Portal](https://portal.gdc.cancer.gov/)
- **PXD039744** — [ProteomeXchange](https://www.proteomexchange.org/)
- **SCP1532** — [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP1532)

---

## Known gaps

These are documented rather than silently patched, so the repository stays
faithful to what was actually executed. **Please read before running.**

## Citation

If you use this code, please cite the manuscript above. *(Add the journal
reference and this repository's DOI once available.)*

## License

Released under the [MIT License](LICENSE).
