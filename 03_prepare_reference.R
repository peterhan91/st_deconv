# =============================================================================
# 03_prepare_reference.R - Download and Prepare scRNA-seq Reference
# ST Deconvolution Pipeline for IDC Data using RCTD
# Reference: Wu et al. 2021 Breast Cancer Atlas (GSE176078)
# =============================================================================

library(Seurat)
library(Matrix)
library(dplyr)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
REFERENCE_DIR <- "reference"
OUTPUT_DIR <- "output/data"

# Create directories if they don't exist
if (!dir.exists(REFERENCE_DIR)) dir.create(REFERENCE_DIR, recursive = TRUE)
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# -----------------------------------------------------------------------------
# Download Reference Data
# -----------------------------------------------------------------------------
# Wu et al. 2021 Breast Cancer scRNA-seq data
# Paper: "A single-cell and spatially resolved atlas of human breast cancers"
# GEO: GSE176078

message("=== Reference Data: Wu et al. 2021 Breast Cancer Atlas ===")
message("
Reference: Wu et al. 2021 'A single-cell and spatially resolved atlas of human breast cancers'
Source: GEO GSE176078
URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078

Download command:
  curl -L -o reference/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz \\
    'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz'

Extract command:
  tar -xzf reference/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz -C reference/
  mv reference/Wu_etal_2021_BRCA_scRNASeq/* reference/
")

# -----------------------------------------------------------------------------
# Check if reference file exists
# -----------------------------------------------------------------------------
ref_file <- file.path(REFERENCE_DIR, "Wu_etal_2021_BRCA_scRNASeq.rds")
ref_file_alt <- file.path(REFERENCE_DIR, "reference_seurat.rds")

# Check for MTX format from GEO
ref_mtx_file <- file.path(REFERENCE_DIR, "count_matrix_sparse.mtx")
ref_metadata_file <- file.path(REFERENCE_DIR, "metadata.csv")

# Function to prepare reference from Seurat object
prepare_reference_from_seurat <- function(ref_seurat,
                                          celltype_col = "celltype_major",
                                          max_cells_per_type = 1000,
                                          min_cells_per_type = 25) {

  message("\n=== Preparing reference data ===")

  # Get cell type annotations
  if (!celltype_col %in% colnames(ref_seurat@meta.data)) {
    # Try alternative column names
    alt_cols <- c("cell_type", "CellType", "celltype", "cluster", "Cluster",
                  "celltype_minor", "annotation", "cell_annotation")
    for (col in alt_cols) {
      if (col %in% colnames(ref_seurat@meta.data)) {
        celltype_col <- col
        message(paste("Using cell type column:", col))
        break
      }
    }
  }

  message(paste("Cell type column:", celltype_col))

  # Get cell types
  celltypes <- ref_seurat@meta.data[[celltype_col]]
  names(celltypes) <- colnames(ref_seurat)

  # Remove cells with NA or empty cell types
  valid_cells <- !is.na(celltypes) & celltypes != "" & celltypes != "NA"
  celltypes <- celltypes[valid_cells]

  message(paste("\nTotal cells after removing NA:", length(celltypes)))
  message("\nCell type distribution:")
  print(table(celltypes))

  # Filter cell types with too few cells
  celltype_counts <- table(celltypes)
  valid_types <- names(celltype_counts)[celltype_counts >= min_cells_per_type]

  message(paste("\nCell types with >=", min_cells_per_type, "cells:", length(valid_types)))

  # Subset to valid cell types
  cells_to_keep <- names(celltypes)[celltypes %in% valid_types]

  # Downsample to max_cells_per_type per cell type for computational efficiency
  set.seed(42)
  downsampled_cells <- c()
  for (ct in valid_types) {
    ct_cells <- names(celltypes)[celltypes == ct]
    if (length(ct_cells) > max_cells_per_type) {
      ct_cells <- sample(ct_cells, max_cells_per_type)
    }
    downsampled_cells <- c(downsampled_cells, ct_cells)
  }

  message(paste("\nCells after downsampling (max", max_cells_per_type, "per type):",
                length(downsampled_cells)))

  # Subset Seurat object
  ref_subset <- subset(ref_seurat, cells = downsampled_cells)

  # Get raw counts
  counts <- GetAssayData(ref_subset, layer = "counts")

  # Ensure counts are integers
  if (!all(counts@x == floor(counts@x))) {
    message("Converting to integer counts...")
    counts@x <- round(counts@x)
  }

  # Get final cell types
  final_celltypes <- factor(ref_subset@meta.data[[celltype_col]])
  names(final_celltypes) <- colnames(ref_subset)

  message("\nFinal cell type distribution:")
  print(table(final_celltypes))

  # Prepare output list
  reference_data <- list(
    counts = counts,
    cell_types = final_celltypes,
    nUMI = colSums(counts)
  )

  return(reference_data)
}

# -----------------------------------------------------------------------------
# Try to load and prepare reference
# -----------------------------------------------------------------------------
reference_data <- NULL

if (file.exists(ref_file)) {
  message(paste("\n=== Loading reference from:", ref_file, "==="))
  ref_seurat <- readRDS(ref_file)
  reference_data <- prepare_reference_from_seurat(ref_seurat)

} else if (file.exists(ref_file_alt)) {
  message(paste("\n=== Loading reference from:", ref_file_alt, "==="))
  ref_seurat <- readRDS(ref_file_alt)
  reference_data <- prepare_reference_from_seurat(ref_seurat)

} else if (file.exists(ref_mtx_file) && file.exists(ref_metadata_file)) {
  # Load MTX format reference from GEO (Wu et al. 2021)
  message("\n=== Loading reference from GEO MTX format ===")
  message("This is the Wu et al. 2021 breast cancer scRNA-seq atlas")

  # Load count matrix
  message("Loading count matrix (this may take a few minutes)...")
  counts <- Matrix::readMM(ref_mtx_file)
  counts <- as(counts, "dgCMatrix")

  # Load barcodes (cell IDs)
  barcodes_file <- file.path(REFERENCE_DIR, "count_matrix_barcodes.tsv")
  barcodes <- read.table(barcodes_file, header = FALSE, stringsAsFactors = FALSE)$V1

  # Load genes
  genes_file <- file.path(REFERENCE_DIR, "count_matrix_genes.tsv")
  genes <- read.table(genes_file, header = FALSE, stringsAsFactors = FALSE)$V1

  # Set row and column names
  rownames(counts) <- genes
  colnames(counts) <- barcodes

  message(paste("Loaded count matrix:", nrow(counts), "genes x", ncol(counts), "cells"))

  # Load metadata
  metadata <- read.csv(ref_metadata_file, row.names = 1, stringsAsFactors = FALSE)
  message(paste("Loaded metadata for", nrow(metadata), "cells"))

  # Use celltype_major for deconvolution
  celltype_col <- "celltype_major"
  celltypes <- metadata[[celltype_col]]
  names(celltypes) <- rownames(metadata)

  message("\nCell type distribution (celltype_major):")
  print(table(celltypes))

  # Subsample for computational efficiency (max cells per type)
  max_cells_per_type <- 2000
  min_cells_per_type <- 25

  set.seed(42)
  celltype_counts <- table(celltypes)
  valid_types <- names(celltype_counts)[celltype_counts >= min_cells_per_type]

  downsampled_cells <- c()
  for (ct in valid_types) {
    ct_cells <- names(celltypes)[celltypes == ct]
    if (length(ct_cells) > max_cells_per_type) {
      ct_cells <- sample(ct_cells, max_cells_per_type)
    }
    downsampled_cells <- c(downsampled_cells, ct_cells)
  }

  message(paste("\nSubsampling to max", max_cells_per_type, "cells per type..."))
  message(paste("Total cells after subsampling:", length(downsampled_cells)))

  # Subset counts and celltypes
  counts <- counts[, downsampled_cells]
  final_celltypes <- factor(celltypes[downsampled_cells])
  names(final_celltypes) <- downsampled_cells

  message("\nFinal cell type distribution:")
  print(table(final_celltypes))

  # Create reference data
  reference_data <- list(
    counts = counts,
    cell_types = final_celltypes,
    nUMI = colSums(counts)
  )

  message("\n=== Wu et al. 2021 reference loaded successfully ===")

} else {
  # No reference data found - stop with instructions
  stop("
=============================================================================
ERROR: Reference data not found!
=============================================================================

Please download the Wu et al. 2021 breast cancer scRNA-seq reference:

  1. Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078

  2. Download: GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz

  3. Extract to the 'reference/' directory:
     tar -xzf GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz -C reference/

  4. Ensure these files exist in reference/:
     - count_matrix_sparse.mtx
     - count_matrix_barcodes.tsv
     - count_matrix_genes.tsv
     - metadata.csv

Then re-run this script.
=============================================================================
")
}

# -----------------------------------------------------------------------------
# Save processed reference
# -----------------------------------------------------------------------------
if (!is.null(reference_data)) {
  ref_output <- file.path(OUTPUT_DIR, "reference_for_rctd.rds")
  saveRDS(reference_data, ref_output)
  message(paste("\n=== Saved processed reference to:", ref_output, "==="))

  # Print summary
  message("\n=== Reference Summary ===")
  message(paste("Total cells:", ncol(reference_data$counts)))
  message(paste("Total genes:", nrow(reference_data$counts)))
  message(paste("Cell types:", length(levels(reference_data$cell_types))))
  message("\nCell types included:")
  for (ct in levels(reference_data$cell_types)) {
    n <- sum(reference_data$cell_types == ct)
    message(paste("  -", ct, ":", n, "cells"))
  }
}

message("\n=== Reference preparation complete! ===")
message("Next step: Run 04_run_rctd.R")

# Return reference data for interactive use
reference_data
