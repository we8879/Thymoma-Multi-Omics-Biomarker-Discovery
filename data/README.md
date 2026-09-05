# Input data

All files below are **shipped with the repository** (25 MB total), so scripts
`01`, `04` and `06` can be re-run without downloading anything. Raw methylation
IDAT files and the single-cell expression matrix are **not** included — see the
root [README](../README.md#data-availability).

## Files

| File                  | Size   | Shape                         | Used by    | Description                                                  |
| --------------------- | ------ | ----------------------------- | ---------- | ------------------------------------------------------------ |
| `GSE177522counts.csv` | 0.8 MB | 24,992 genes × **7** samples  | `01`, `06` | Raw counts, GSE177522 (2 thymoma / 5 normal). Columns `SRR14794294`–`SRR14794301` |
| `GSE79978counts.csv`  | 2.1 MB | 39,374 genes × **16** samples | `01`, `06` | Raw counts, GSE79978 (13 tumour / 3 normal). Columns `GSM2109575`–`GSM2109590` |
| `TCGA_expression.csv` | 21 MB  | 59,427 rows × **122** samples | `06`       | TCGA-THYM expression. First two columns are `gene_symbol` and `gene_type`; sample columns are prefixed `Normal_` or `FullData_` |
| `Training_class.txt`  | 0.6 KB | **23** rows                   | `06`       | Training labels: `sample`, `outcome` (0 = normal n=8, 1 = tumour n=15), `batch` |
| `Testing_class.txt`   | 5.6 KB | **122** rows                  | `06`       | TCGA labels: `sample`, `outcome` (0 = normal n=2, 1 = tumour n=120) |
| `proteomics.xlsx`     | 0.9 MB | Sheet1, 31 sample columns     | `03`       | PXD039744 protein matrix. Columns: `Protein Accession`, `geneSymbol`, then samples |
| `methods.xlsx`        | 14 KB  | 107 model names               | `06`       | Model combination list (see caveat below)                    |

## Sources

- **GSE79978**, **GSE177522**, **GSE218549** — [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/)
- **TCGA-THYM** — [GDC Portal](https://portal.gdc.cancer.gov/)
- **PXD039744** — [ProteomeXchange](https://www.proteomexchange.org/)
- **SCP1532** — [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP1532)
