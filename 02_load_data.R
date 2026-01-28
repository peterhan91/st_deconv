# =============================================================================
# 02_load_data.R - Load and Convert h5ad to R Format
# ST Deconvolution Pipeline for IDC Data using RCTD
# =============================================================================

library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)
library(Matrix)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
INPUT_FILE <- "NCBI783.h5ad"
OUTPUT_DIR <- "output/data"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Step 1: Read h5ad file using zellkonverter
# -----------------------------------------------------------------------------
message("=== Loading h5ad file ===")
message(paste("Input file:", INPUT_FILE))

# Read the h5ad file into a SingleCellExperiment object
sce <- readH5AD(INPUT_FILE)

message(paste("Loaded SCE object with", ncol(sce), "spots and", nrow(sce), "genes"))

# -----------------------------------------------------------------------------
# Step 2: Explore the data structure
# -----------------------------------------------------------------------------
message("\n=== Data structure ===")

# Print basic info
message(paste("Number of spots:", ncol(sce)))
message(paste("Number of genes:", nrow(sce)))

# Check available assays
message("\nAvailable assays:")
print(assayNames(sce))

# Check available metadata
message("\nColumn data (spot metadata):")
print(head(colData(sce)))

message("\nRow data (gene metadata):")
print(head(rowData(sce)))

# Check for spatial coordinates
if ("spatial" %in% reducedDimNames(sce)) {
  message("\nSpatial coordinates found in reducedDims")
  coords <- reducedDim(sce, "spatial")
  print(head(coords))
}

# -----------------------------------------------------------------------------
# Step 3: Convert to Seurat object
# -----------------------------------------------------------------------------
message("\n=== Converting to Seurat object ===")

# Get the count matrix
counts <- assay(sce, "X")

# Ensure it's a sparse matrix with proper format
if (!inherits(counts, "dgCMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

# Transpose if genes are in columns (h5ad stores cells x genes)
# Check dimensions to determine orientation
if (ncol(counts) == nrow(sce)) {
  # Genes are in rows already
  message("Count matrix is in gene x spot format")
} else if (nrow(counts) == ncol(sce)) {
  # Need to transpose (spots x genes -> genes x spots)
  message("Transposing count matrix from spot x gene to gene x spot format")
  counts <- t(counts)
}

# Set row and column names
rownames(counts) <- rownames(sce)
colnames(counts) <- colnames(sce)

# Create Seurat object
st_data <- CreateSeuratObject(
  counts = counts,
  assay = "Spatial",
  meta.data = as.data.frame(colData(sce))
)

message(paste("Created Seurat object with", ncol(st_data), "spots and", nrow(st_data), "genes"))

# -----------------------------------------------------------------------------
# Step 4: Add spatial coordinates
# -----------------------------------------------------------------------------
message("\n=== Adding spatial coordinates ===")

# Try to get spatial coordinates from different possible locations
coords <- NULL

# Option 1: From reducedDims
if ("spatial" %in% reducedDimNames(sce)) {
  coords <- reducedDim(sce, "spatial")
  message("Found coordinates in reducedDims$spatial")
}

# Option 2: From obsm (stored in metadata)
if (is.null(coords) && "obsm" %in% names(metadata(sce))) {
  if ("spatial" %in% names(metadata(sce)$obsm)) {
    coords <- metadata(sce)$obsm$spatial
    message("Found coordinates in obsm$spatial")
  }
}

# Option 3: From colData
if (is.null(coords)) {
  coord_cols <- c("x", "y", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
  available_cols <- intersect(coord_cols, colnames(colData(sce)))
  if (length(available_cols) >= 2) {
    coords <- as.matrix(colData(sce)[, available_cols[1:2]])
    message(paste("Found coordinates in colData:", paste(available_cols[1:2], collapse = ", ")))
  }
}

if (!is.null(coords)) {
  # Ensure coords is a matrix with proper format
  coords <- as.data.frame(coords)
  colnames(coords) <- c("x", "y")[1:ncol(coords)]
  rownames(coords) <- colnames(st_data)

  # Add to Seurat metadata
  st_data$x_coord <- coords[, 1]
  st_data$y_coord <- coords[, 2]

  # Create a basic spatial image slot for visualization
  # This is a simplified version - full Visium integration would need image data
  centroids <- CreateCentroids(coords)
  fov <- CreateFOV(
    coords = centroids,
    type = "centroids",
    assay = "Spatial"
  )
  st_data[["fov"]] <- fov

  message("Spatial coordinates added to Seurat object")
  message(paste("X range:", round(min(coords[,1]), 2), "-", round(max(coords[,1]), 2)))
  message(paste("Y range:", round(min(coords[,2]), 2), "-", round(max(coords[,2]), 2)))
} else {
  warning("No spatial coordinates found! Visualization may be limited.")
}

# -----------------------------------------------------------------------------
# Step 5: Basic QC
# -----------------------------------------------------------------------------
message("\n=== Basic QC metrics ===")

# Calculate QC metrics
st_data[["nCount_Spatial"]] <- colSums(GetAssayData(st_data, layer = "counts"))
st_data[["nFeature_Spatial"]] <- colSums(GetAssayData(st_data, layer = "counts") > 0)

message(paste("Total UMI range:", min(st_data$nCount_Spatial), "-", max(st_data$nCount_Spatial)))
message(paste("Genes detected range:", min(st_data$nFeature_Spatial), "-", max(st_data$nFeature_Spatial)))
message(paste("Median UMI per spot:", median(st_data$nCount_Spatial)))
message(paste("Median genes per spot:", median(st_data$nFeature_Spatial)))

# -----------------------------------------------------------------------------
# Step 6: Save processed data
# -----------------------------------------------------------------------------
message("\n=== Saving processed data ===")

# Save Seurat object
seurat_output <- file.path(OUTPUT_DIR, "st_data_seurat.rds")
saveRDS(st_data, seurat_output)
message(paste("Saved Seurat object to:", seurat_output))

# Save coordinates separately for RCTD
if (!is.null(coords)) {
  coords_output <- file.path(OUTPUT_DIR, "spatial_coords.csv")
  write.csv(coords, coords_output)
  message(paste("Saved spatial coordinates to:", coords_output))
}

# Save count matrix for RCTD
counts_output <- file.path(OUTPUT_DIR, "counts_matrix.rds")
saveRDS(GetAssayData(st_data, layer = "counts"), counts_output)
message(paste("Saved count matrix to:", counts_output))

message("\n=== Data loading complete! ===")
message("Next step: Run 03_prepare_reference.R")

# Return the Seurat object for interactive use
st_data
