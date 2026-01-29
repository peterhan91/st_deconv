# =============================================================================
# 03_prepare_bulk_data.R - Prepare bulk and reference data for MuSiC
# Bulk Deconvolution Pipeline using MuSiC
# =============================================================================

library(Biobase)
library(Matrix)
library(dplyr)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
DATA_DIR <- "data"
REF_DIR <- "../reference"
OUTPUT_DIR <- "output"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(file.path(OUTPUT_DIR, "data"))) dir.create(file.path(OUTPUT_DIR, "data"))

# -----------------------------------------------------------------------------
# Load Wu et al. scRNA-seq Reference (same as ST pipeline)
# -----------------------------------------------------------------------------
message("=== Loading Wu et al. scRNA-seq reference ===")

# Check if reference files exist
ref_files <- c(
  matrix = file.path(REF_DIR, "count_matrix_sparse.mtx"),
  barcodes = file.path(REF_DIR, "count_matrix_barcodes.tsv"),
  genes = file.path(REF_DIR, "count_matrix_genes.tsv"),
  metadata = file.path(REF_DIR, "metadata.csv")
)

if (!all(file.exists(ref_files))) {
  stop("
=============================================================================
ERROR: Reference data not found!

Please download Wu et al. 2021 breast cancer scRNA-seq reference:

  curl -L -o ../reference/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz \\
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz'

  tar -xzf ../reference/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz -C ../reference/
  mv ../reference/Wu_etal_2021_BRCA_scRNASeq/* ../reference/
=============================================================================
")
}

# Load sparse matrix
message("Loading count matrix...")
counts <- Matrix::readMM(ref_files["matrix"])
counts <- as(counts, "CsparseMatrix")

# Load barcodes and genes
barcodes <- read.table(ref_files["barcodes"], header = FALSE, stringsAsFactors = FALSE)$V1
genes <- read.table(ref_files["genes"], header = FALSE, stringsAsFactors = FALSE)$V1

rownames(counts) <- genes
colnames(counts) <- barcodes

# Load metadata
metadata <- read.csv(ref_files["metadata"], row.names = 1)
message(paste("Loaded:", nrow(counts), "genes x", ncol(counts), "cells"))

# Align metadata with counts
common_cells <- intersect(colnames(counts), rownames(metadata))
counts <- counts[, common_cells]
metadata <- metadata[common_cells, ]

message(paste("After alignment:", ncol(counts), "cells"))
message("\nCell type distribution:")
print(table(metadata$celltype_major))

# -----------------------------------------------------------------------------
# Subsample reference for efficiency (MuSiC can handle more cells than RCTD)
# -----------------------------------------------------------------------------
message("\n=== Subsampling reference ===")

set.seed(42)
max_cells_per_type <- 1000  # Reduced for memory efficiency

cell_types <- unique(metadata$celltype_major)
sampled_cells <- c()

for (ct in cell_types) {
  ct_cells <- rownames(metadata)[metadata$celltype_major == ct]
  n_sample <- min(length(ct_cells), max_cells_per_type)
  sampled_cells <- c(sampled_cells, sample(ct_cells, n_sample))
}

counts_sub <- counts[, sampled_cells]
metadata_sub <- metadata[sampled_cells, ]

message(paste("Subsampled to:", ncol(counts_sub), "cells"))
message("\nSubsampled cell type distribution:")
print(table(metadata_sub$celltype_major))

# -----------------------------------------------------------------------------
# Create ExpressionSet for MuSiC (scRNA-seq reference)
# -----------------------------------------------------------------------------
message("\n=== Creating ExpressionSet for reference ===")

# MuSiC requires ExpressionSet objects
# Convert sparse to dense (MuSiC works better with dense for reference)
counts_dense <- as.matrix(counts_sub)

# Create phenoData
# Use orig.ident as subject/patient identifier (patient_id doesn't exist in this dataset)
pdata <- new("AnnotatedDataFrame", data = data.frame(
  cellType = metadata_sub$celltype_major,
  SubjectName = metadata_sub$orig.ident,
  row.names = colnames(counts_dense)
))

# Create ExpressionSet
sc_eset <- ExpressionSet(
  assayData = counts_dense,
  phenoData = pdata
)

message(paste("Created scRNA-seq ExpressionSet:",
              nrow(sc_eset), "genes x", ncol(sc_eset), "cells"))

# Save reference ExpressionSet
saveRDS(sc_eset, file.path(OUTPUT_DIR, "data", "sc_reference_eset.rds"))
message("Saved: output/data/sc_reference_eset.rds")

# -----------------------------------------------------------------------------
# Load and Prepare Bulk Data
# -----------------------------------------------------------------------------
message("\n=== Loading bulk expression data ===")

# Check for different possible bulk data sources (NO demo/simulated data)
# Supports CPTAC or TCGA-BRCA data
bulk_files <- c(
  cptac_rds = file.path(DATA_DIR, "cptac_brca_expression.rds"),
  cptac_csv = file.path(DATA_DIR, "cptac_brca_expression.csv"),
  cptac_txt = file.path(DATA_DIR, "data_mrna_seq_v2_rsem.txt"),
  tcga_rds = file.path(DATA_DIR, "tcga_brca_expression.rds"),
  tcga_se = file.path(DATA_DIR, "tcga_brca_se.rds")
)

bulk_data <- NULL
bulk_source <- NULL

# Try loading in order of preference (real data only)
if (file.exists(bulk_files["cptac_rds"])) {
  message("Loading CPTAC data from RDS...")
  bulk_data <- readRDS(bulk_files["cptac_rds"])
  bulk_source <- "CPTAC"

} else if (file.exists(bulk_files["cptac_csv"])) {
  message("Loading CPTAC data from CSV...")
  bulk_data <- read.csv(bulk_files["cptac_csv"], row.names = 1)
  bulk_data <- as.matrix(bulk_data)
  bulk_source <- "CPTAC"

} else if (file.exists(bulk_files["cptac_txt"])) {
  message("Loading cBioPortal CPTAC data...")
  bulk_raw <- read.delim(bulk_files["cptac_txt"], row.names = 1)
  # Remove Hugo_Symbol column if present
  if ("Hugo_Symbol" %in% colnames(bulk_raw)) {
    bulk_raw <- bulk_raw[, -which(colnames(bulk_raw) == "Hugo_Symbol")]
  }
  # Remove Entrez_Gene_Id if present
  if ("Entrez_Gene_Id" %in% colnames(bulk_raw)) {
    bulk_raw <- bulk_raw[, -which(colnames(bulk_raw) == "Entrez_Gene_Id")]
  }
  bulk_data <- as.matrix(bulk_raw)
  bulk_source <- "cBioPortal"

} else if (file.exists(bulk_files["tcga_rds"])) {
  message("Loading TCGA-BRCA data from RDS...")
  bulk_data <- readRDS(bulk_files["tcga_rds"])
  if (inherits(bulk_data, "SummarizedExperiment") || inherits(bulk_data, "RangedSummarizedExperiment")) {
    # Extract expression matrix from SummarizedExperiment
    bulk_data <- SummarizedExperiment::assay(bulk_data, "unstranded")
    # Convert Ensembl IDs to gene symbols if rowData contains gene names
  }
  bulk_data <- as.matrix(bulk_data)
  bulk_source <- "TCGA-BRCA"

} else if (file.exists(bulk_files["tcga_se"])) {
  message("Loading TCGA-BRCA SummarizedExperiment...")
  library(SummarizedExperiment)
  se <- readRDS(bulk_files["tcga_se"])
  # Extract expression matrix
  bulk_data <- assay(se, "unstranded")
  bulk_data <- as.matrix(bulk_data)
  bulk_source <- "TCGA-BRCA"

} else {
  stop("
=============================================================================
ERROR: No bulk expression data found!

Please download CPTAC BRCA data:

Option 1: Run 02_download_cptac.R first

Option 2: Manual download from cBioPortal:
  1. Go to: https://www.cbioportal.org/study/summary?id=brca_cptac_2020
  2. Click 'Download' -> 'All Data'
  3. Extract data_mrna_seq_v2_rsem.txt to bulk_deconv/data/

Option 3: Manual download from GDC:
  1. Go to: https://portal.gdc.cancer.gov/
  2. Search for CPTAC-3 breast cancer RNA-seq
  3. Download and process the data
=============================================================================
")
}

message(paste("Loaded bulk data:", nrow(bulk_data), "genes x", ncol(bulk_data), "samples"))
message(paste("Data source:", bulk_source))

# -----------------------------------------------------------------------------
# Preprocess Bulk Data
# -----------------------------------------------------------------------------
message("\n=== Preprocessing bulk data ===")

# Remove genes with all zeros
nonzero_genes <- rowSums(bulk_data) > 0
bulk_data <- bulk_data[nonzero_genes, ]
message(paste("After removing zero genes:", nrow(bulk_data), "genes"))

# Find common genes between bulk and reference
common_genes <- intersect(rownames(bulk_data), rownames(sc_eset))
message(paste("Common genes with reference:", length(common_genes)))

if (length(common_genes) < 500) {
  warning("Low number of common genes. Results may be less reliable.")
}

# Subset both to common genes
bulk_data <- bulk_data[common_genes, ]
sc_eset <- sc_eset[common_genes, ]

message(paste("Final dimensions - Bulk:", nrow(bulk_data), "x", ncol(bulk_data)))
message(paste("Final dimensions - Reference:", nrow(sc_eset), "x", ncol(sc_eset)))

# -----------------------------------------------------------------------------
# Create ExpressionSet for Bulk Data
# -----------------------------------------------------------------------------
message("\n=== Creating ExpressionSet for bulk data ===")

# Create simple phenoData for bulk samples
bulk_pdata <- new("AnnotatedDataFrame", data = data.frame(
  SampleID = colnames(bulk_data),
  row.names = colnames(bulk_data)
))

# Create ExpressionSet
bulk_eset <- ExpressionSet(
  assayData = as.matrix(bulk_data),
  phenoData = bulk_pdata
)

message(paste("Created bulk ExpressionSet:",
              nrow(bulk_eset), "genes x", ncol(bulk_eset), "samples"))

# Save bulk ExpressionSet
saveRDS(bulk_eset, file.path(OUTPUT_DIR, "data", "bulk_eset.rds"))
message("Saved: output/data/bulk_eset.rds")

# Also update and save reference with common genes
saveRDS(sc_eset, file.path(OUTPUT_DIR, "data", "sc_reference_eset_filtered.rds"))
message("Saved: output/data/sc_reference_eset_filtered.rds")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
message("\n", strrep("=", 60))
message("DATA PREPARATION COMPLETE")
message(strrep("=", 60))
message(paste("\nBulk data source:", bulk_source))
message(paste("Bulk samples:", ncol(bulk_eset)))
message(paste("Reference cells:", ncol(sc_eset)))
message(paste("Cell types:", length(unique(pData(sc_eset)$cellType))))
message(paste("Common genes:", nrow(sc_eset)))
message("\nCell types in reference:")
print(table(pData(sc_eset)$cellType))

message("\n=== Next: Run 04_run_music.R ===")
