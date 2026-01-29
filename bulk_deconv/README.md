# Bulk RNA-seq Deconvolution Pipeline using MuSiC

Deconvolve bulk RNA-seq data (CPTAC BRCA) into cell type proportions using **MuSiC** with the Wu et al. 2021 breast cancer single-cell reference.

---

## Overview

| Component | Description |
|-----------|-------------|
| **Method** | MuSiC (Multi-subject Single Cell deconvolution) |
| **Bulk Data** | CPTAC BRCA RNA-seq |
| **Reference** | Wu et al. 2021 breast cancer scRNA-seq atlas (GEO: GSE176078) |
| **Output** | Cell type proportions per sample |

---

## Why MuSiC?

MuSiC was chosen over other bulk deconvolution methods because:

1. **Designed for scRNA-seq references** - Unlike CIBERSORT which uses signature matrices
2. **Handles cross-subject variability** - Uses multiple subjects from reference to improve robustness
3. **Gene weighting** - Weights genes by cross-subject consistency
4. **No registration required** - Unlike CIBERSORTx

---

## Requirements

### R Packages

```r
# Run 01_setup_bulk.R to install all packages
source("01_setup_bulk.R")
```

Key packages:
- MuSiC (from GitHub)
- Biobase
- GenomicDataCommons (for CPTAC download)

---

## Usage

### Step 1: Install Packages
```bash
cd bulk_deconv
Rscript 01_setup_bulk.R
```

### Step 2: Download Bulk RNA-seq Data

**Option A: TCGA-BRCA (Recommended - Largest dataset)**
```r
# Using TCGAbiolinks in R
library(TCGAbiolinks)
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)
GDCdownload(query)
data <- GDCprepare(query)
saveRDS(data, "bulk_deconv/data/tcga_brca_se.rds")
```

**Option B: cBioPortal (Easiest)**
1. Go to: https://www.cbioportal.org/study/summary?id=brca_cptac_2020
2. Click "Download" -> "Download All Data"
3. Extract and copy `data_mrna_seq_v2_rsem.txt` to `bulk_deconv/data/`

**Option C: GDC Data Portal**
1. Go to: https://portal.gdc.cancer.gov/
2. Cases -> Project: CPTAC-3 or TCGA-BRCA
3. Filter: Primary Site = Breast
4. Files -> Data Category = Transcriptome Profiling
5. Add to Cart -> Download

### Step 3: Prepare Data
```bash
Rscript 03_prepare_bulk_data.R
```

### Step 4: Run MuSiC Deconvolution
```bash
Rscript 04_run_music.R
```

### Step 5: Visualize Results
```bash
Rscript 05_visualize_bulk.R
```

---

## Data Sources

### CPTAC BRCA or TCGA-BRCA (Bulk RNA-seq)
- **CPTAC Source**: Clinical Proteomic Tumor Analysis Consortium
- **TCGA Source**: The Cancer Genome Atlas
- **Portal**: GDC Data Portal or cBioPortal
- **Study ID**: brca_cptac_2020 (cBioPortal) or TCGA-BRCA (GDC)

**Note**: Both CPTAC-BRCA and TCGA-BRCA are suitable for bulk deconvolution with MuSiC.

### scRNA-seq Reference (Wu et al. 2021)
- **GEO**: GSE176078
- **Paper**: "A single-cell and spatially resolved atlas of human breast cancers"
- **Cells**: ~100,000 breast cancer cells
- **Cell Types**: 9 major types

The same reference used for the spatial transcriptomics RCTD pipeline is used here.

---

## Output Files

### Data (`output/data/`)

| File | Description |
|------|-------------|
| `music_proportions.csv` | Cell type proportions (samples × cell types) |
| `music_proportions.rds` | Same in R format |
| `deconvolution_summary.csv` | Summary statistics |
| `model_fit_rsquared.csv` | R² values for each sample |

### Figures (`output/figures/`)

| File | Description |
|------|-------------|
| `deconvolution_summary.png` | Combined overview figure |
| `mean_proportions_barplot.png` | Mean cell type proportions |
| `proportions_boxplot.png` | Distribution across samples |
| `composition_heatmap.png` | Sample × cell type heatmap |
| `celltype_correlation.png` | Cell type correlation matrix |

---

## Expected Results

For CPTAC BRCA, expected cell type proportions:

| Cell Type | Expected Range |
|-----------|----------------|
| CAFs | 15-30% |
| Cancer Epithelial | 10-25% |
| T-cells | 10-20% |
| Myeloid | 8-15% |
| Endothelial | 5-15% |
| Normal Epithelial | 5-15% |
| B-cells | 2-8% |
| PVL | 2-8% |
| Plasmablasts | 1-5% |

---

## Quality Assessment

MuSiC provides R² values indicating how well the estimated mixture fits the observed bulk expression:

- **R² > 0.7**: Good fit
- **R² 0.5-0.7**: Acceptable
- **R² < 0.5**: Poor fit (may indicate missing cell types or batch effects)

---

## Comparison with ST Deconvolution

| Aspect | Bulk (MuSiC) | Spatial (RCTD) |
|--------|--------------|----------------|
| Input | Bulk RNA-seq | Visium spots |
| Resolution | Whole sample | ~10-100 cells/spot |
| Spatial info | No | Yes |
| Sample size | Many samples | Single tissue section |
| Use case | Cohort analysis | Spatial architecture |

Both pipelines use the same Wu et al. 2021 reference, enabling direct comparison.

---

## References

### MuSiC Method
> Wang X, et al. (2019). Bulk tissue cell type deconvolution with multi-subject single-cell expression reference. *Nature Communications*, 10, 380.

### Reference Data
> Wu SZ, et al. (2021). A single-cell and spatially resolved atlas of human breast cancers. *Nature Genetics*, 53, 1334-1347.

---

## Troubleshooting

### "No bulk expression data found"
Download CPTAC data following instructions in Step 2, or run `02_download_cptac.R`.

### "Reference data not found"
The scRNA-seq reference is shared with the ST pipeline. Download from GEO GSE176078:
```bash
curl -L -o ../reference/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz \
  'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz'
tar -xzf ../reference/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz -C ../reference/
mv ../reference/Wu_etal_2021_BRCA_scRNASeq/* ../reference/
```

### Low number of common genes
Ensure gene names match between bulk and reference (typically gene symbols).
