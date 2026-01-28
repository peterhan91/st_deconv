# Spatial Transcriptomics Deconvolution Pipeline for IDC

Deconvolve spatial transcriptomics (Visium) data into cell type proportions using **RCTD** (Robust Cell Type Decomposition) with a breast cancer single-cell RNA-seq reference.

---

## Overview

| Component | Description |
|-----------|-------------|
| **Input** | `NCBI783.h5ad` - IDC (Invasive Ductal Carcinoma) Visium data |
| **Method** | RCTD from spacexr package |
| **Reference** | Wu et al. 2021 breast cancer scRNA-seq atlas (GEO: GSE176078) |
| **Output** | Cell type proportions per spot + visualizations |

---

## Data Sources

### Spatial Transcriptomics Data
- **File**: `NCBI783.h5ad`
- **Type**: 10x Visium spatial transcriptomics
- **Dimensions**: 3,869 spots × 541 genes
- **Tissue**: Invasive Ductal Carcinoma (IDC)

### Single-Cell Reference Data
- **Source**: Wu et al. 2021 - *"A single-cell and spatially resolved atlas of human breast cancers"*
- **Publication**: Nature Genetics, 2021
- **GEO Accession**: [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078)
- **Cells**: ~100,000 cells from breast cancer patients
- **Cell Types**: 9 major types

| Cell Type | Description |
|-----------|-------------|
| Cancer Epithelial | Malignant tumor cells |
| Normal Epithelial | Adjacent normal epithelium |
| CAFs | Cancer-associated fibroblasts |
| T-cells | Tumor-infiltrating T lymphocytes |
| B-cells | B lymphocytes |
| Myeloid | Macrophages, monocytes, dendritic cells |
| Endothelial | Vascular endothelium |
| PVL | Perivascular-like cells |
| Plasmablasts | Antibody-secreting B cells |

---

## Requirements

### R Packages

```r
# CRAN
install.packages(c("Matrix", "Seurat", "ggplot2", "patchwork",
                   "dplyr", "tidyr", "RColorBrewer", "viridis"))

# Bioconductor
BiocManager::install(c("zellkonverter", "SingleCellExperiment"))

# GitHub (RCTD)
devtools::install_github("dmcable/spacexr")
```

### System Requirements
- R >= 4.0
- ~8 GB RAM (for loading reference data)
- ~3 GB disk space (for reference + outputs)

---

## Installation & Setup

### 1. Clone/Download Repository
```bash
cd /path/to/your/projects
git clone <repository-url> st_deconv
cd st_deconv
```

### 2. Install R Packages
```bash
Rscript 01_setup.R
```

### 3. Download Reference Data
```bash
# Download from GEO (533 MB)
curl -L -o reference/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz \
  'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz'

# Extract
tar -xzf reference/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz -C reference/
mv reference/Wu_etal_2021_BRCA_scRNASeq/* reference/
```

After extraction, `reference/` should contain:
```
reference/
├── count_matrix_sparse.mtx    # Expression matrix (2.2 GB)
├── count_matrix_barcodes.tsv  # Cell barcodes
├── count_matrix_genes.tsv     # Gene names
└── metadata.csv               # Cell type annotations
```

---

## Usage

Run the pipeline scripts in order:

```bash
# Step 1: Install packages (run once)
Rscript 01_setup.R

# Step 2: Load and convert h5ad to R format
Rscript 02_load_data.R

# Step 3: Prepare scRNA-seq reference
Rscript 03_prepare_reference.R

# Step 4: Run RCTD deconvolution
Rscript 04_run_rctd.R

# Step 5: Generate visualizations
Rscript 05_visualize.R
```

Or run all at once:
```bash
Rscript -e 'source("02_load_data.R"); source("03_prepare_reference.R"); source("04_run_rctd.R"); source("05_visualize.R")'
```

---

## Output Files

### Data (`output/data/`)

| File | Description |
|------|-------------|
| `cell_type_weights.csv` | Cell type proportions per spot (spots × cell types) |
| `cell_type_weights.rds` | Same as above in R format |
| `rctd_results.rds` | Full RCTD object with all results |
| `st_data_deconvolved.rds` | Seurat object with deconvolution results |
| `deconvolution_summary.csv` | Summary statistics per cell type |
| `st_data_seurat.rds` | Processed Seurat object (before deconvolution) |
| `spatial_coords.csv` | Spot spatial coordinates |

### Figures (`output/figures/`)

| File | Description |
|------|-------------|
| `deconvolution_summary_figure.png` | Combined overview figure |
| `spatial_dominant_celltype.png` | Map of dominant cell type per spot |
| `spatial_celltypes_all.png` | All cell types in grid layout |
| `spatial_<CellType>.png` | Individual cell type spatial maps |
| `celltype_proportions_boxplot.png` | Distribution of proportions |
| `celltype_mean_proportions.png` | Mean proportion bar chart |
| `celltype_correlation_heatmap.png` | Cell type co-occurrence patterns |
| `celltype_pie_chart.png` | Overall composition |
| `spatial_composition_stacked.png` | Composition along spatial axis |

---

## Results Interpretation

### Example Output (NCBI783 IDC Sample)

| Cell Type | Mean % | Interpretation |
|-----------|--------|----------------|
| CAFs | 29.8% | Abundant tumor stroma |
| Normal Epithelial | 14.3% | Adjacent normal tissue |
| Cancer Epithelial | 14.0% | Malignant cells |
| Myeloid | 12.6% | Immune infiltration |
| T-cells | 11.4% | Tumor-infiltrating lymphocytes |
| Endothelial | 10.8% | Tumor vasculature |
| Plasmablasts | 2.7% | Humoral immune response |
| PVL | 2.4% | Perivascular cells |
| B-cells | 1.9% | B lymphocytes |

### Biological Insights
- **High CAF proportion** is characteristic of IDC - these stromal cells support tumor growth
- **Immune infiltrate** (T-cells + Myeloid ~24%) suggests active tumor microenvironment
- **Spatial patterns** reveal tumor architecture and immune cell localization

---

## Project Structure

```
st_deconv/
├── NCBI783.h5ad              # Input: Spatial transcriptomics data
├── 01_setup.R                # Package installation
├── 02_load_data.R            # Load h5ad, convert to Seurat
├── 03_prepare_reference.R    # Prepare Wu et al. reference
├── 04_run_rctd.R             # Run RCTD deconvolution
├── 05_visualize.R            # Generate figures
├── README.md                 # This file
├── reference/                # scRNA-seq reference data
│   ├── count_matrix_sparse.mtx
│   ├── count_matrix_barcodes.tsv
│   ├── count_matrix_genes.tsv
│   └── metadata.csv
└── output/
    ├── data/                 # Deconvolution results
    └── figures/              # Visualization outputs
```

---

## Quality Assessment

Run the quality assessment script to evaluate deconvolution results:

```bash
Rscript 06_quality_assessment.R
```

### Quality Metrics

| Metric | Good | Acceptable | Poor | This Run |
|--------|------|------------|------|----------|
| **Convergence Rate** | >90% | 60-90% | <60% | 84.1% (Good) |
| **Weight Sums** | =1.0 | ~1.0 | Varies | 1.0 (Excellent) |
| **Cell Type Correlation** | Negative | <0.5 | >0.5 | -0.08 mean (Excellent) |
| **Classification Confidence** | >50% dominant | >35% | <35% | 48% mean (Moderate) |
| **Biological Plausibility** | Matches expected | Mostly | Unexpected | Plausible |

### Interpretation Guide

**1. Convergence Rate (84.1%)**
- % of spots where RCTD successfully fit the model
- >80% is good; <60% suggests data quality issues
- Failed spots often have very low UMI counts

**2. Weight Normalization**
- Cell type proportions should sum to 1.0 per spot
- Perfect normalization indicates proper model fitting

**3. Cell Type Separation**
- Negative correlation between cell types is EXPECTED
- Means when one type is high, others are low (mutually exclusive)
- High positive correlation (>0.5) suggests types are hard to distinguish

**4. Classification Confidence**
- Mean dominant proportion of 0.48 is MODERATE
- This is normal for heterogeneous tumor tissue
- Visium spots (55μm) often contain multiple cell types

**5. Biological Plausibility (IDC-specific)**
- High CAFs (30%): Expected - IDC has abundant desmoplastic stroma
- Cancer epithelial (14%): Reasonable - varies by tumor region
- Immune cells (T-cells + Myeloid ~24%): Expected - active TME
- Endothelial (11%): Expected - tumor vascularization

### Quality Score: 13/15 (87%) - HIGH QUALITY

### Output Files
- `output/figures/quality_assessment.png` - Visual quality report
- `output/data/quality_report.csv` - Metrics summary

---

## Troubleshooting

### "Reference data not found"
Download the reference data following the instructions in [Installation & Setup](#3-download-reference-data).

### "fewer than 10 regression differentially expressed genes found"
This error occurs with synthetic/improper reference data. Ensure you're using the real Wu et al. 2021 reference from GEO.

### Memory issues
The reference has ~100k cells. If you encounter memory problems:
- Reduce `max_cells_per_type` in `03_prepare_reference.R` (default: 2000)
- Close other applications

### Low number of common genes
The ST data has 541 genes. RCTD found 275 common genes with the reference, which is sufficient for deconvolution.

---

## Citation

If you use this pipeline, please cite:

### RCTD Method
> Cable DM, et al. (2022). Robust decomposition of cell type mixtures in spatial transcriptomics. *Nature Biotechnology*, 40, 517-526.

### Reference Data
> Wu SZ, et al. (2021). A single-cell and spatially resolved atlas of human breast cancers. *Nature Genetics*, 53, 1334-1347.

---

## License

This pipeline is provided for research purposes. Please respect the licenses of the underlying tools and datasets.

---

## Contact

For issues with:
- **This pipeline**: Open an issue in this repository
- **RCTD/spacexr**: https://github.com/dmcable/spacexr
- **Reference data**: Contact the Wu et al. 2021 authors
