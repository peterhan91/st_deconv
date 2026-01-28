# =============================================================================
# 04_run_rctd.R - Run RCTD Deconvolution
# ST Deconvolution Pipeline for IDC Data using RCTD
# =============================================================================

library(spacexr)
library(Matrix)
library(Seurat)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
OUTPUT_DIR <- "output/data"
N_CORES <- 4  # Adjust based on your system

# RCTD parameters
DOUBLET_MODE <- "full"  # "full" for multiple cell types per spot (recommended for Visium)
                        # "doublet" for max 2 cell types per spot
                        # "multi" for multiple cell types with different algorithm

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
message("=== Loading data ===")

# Load ST data
st_seurat_file <- file.path(OUTPUT_DIR, "st_data_seurat.rds")
st_counts_file <- file.path(OUTPUT_DIR, "counts_matrix.rds")
coords_file <- file.path(OUTPUT_DIR, "spatial_coords.csv")

if (!file.exists(st_counts_file)) {
  stop("ST count matrix not found. Please run 02_load_data.R first.")
}

st_counts <- readRDS(st_counts_file)
message(paste("Loaded ST counts:", nrow(st_counts), "genes x", ncol(st_counts), "spots"))

# Load coordinates
if (file.exists(coords_file)) {
  coords <- read.csv(coords_file, row.names = 1)
  message(paste("Loaded spatial coordinates for", nrow(coords), "spots"))
} else {
  # Create dummy coordinates if not available
  message("No coordinates file found. Creating dummy coordinates...")
  coords <- data.frame(
    x = seq_len(ncol(st_counts)),
    y = rep(1, ncol(st_counts)),
    row.names = colnames(st_counts)
  )
}

# Ensure coordinate column names
if (!all(c("x", "y") %in% colnames(coords))) {
  colnames(coords)[1:2] <- c("x", "y")
}

# Load reference
ref_file <- file.path(OUTPUT_DIR, "reference_for_rctd.rds")
if (!file.exists(ref_file)) {
  stop("Reference data not found. Please run 03_prepare_reference.R first.")
}

reference_data <- readRDS(ref_file)
message(paste("Loaded reference:", nrow(reference_data$counts), "genes x",
              ncol(reference_data$counts), "cells"))
message(paste("Cell types:", length(levels(reference_data$cell_types))))

# -----------------------------------------------------------------------------
# Find common genes
# -----------------------------------------------------------------------------
message("\n=== Finding common genes ===")

common_genes <- intersect(rownames(st_counts), rownames(reference_data$counts))
message(paste("Common genes between ST and reference:", length(common_genes)))

if (length(common_genes) < 100) {
  warning("Very few common genes found! Check gene name formats.")
  message("ST gene examples:", paste(head(rownames(st_counts), 5), collapse = ", "))
  message("Reference gene examples:", paste(head(rownames(reference_data$counts), 5), collapse = ", "))
}

# Subset to common genes
st_counts_common <- st_counts[common_genes, ]
ref_counts_common <- reference_data$counts[common_genes, ]

# -----------------------------------------------------------------------------
# Create RCTD Reference object
# -----------------------------------------------------------------------------
message("\n=== Creating RCTD Reference object ===")

# Ensure counts are in correct format (integer)
ref_counts_common <- as(ref_counts_common, "dgCMatrix")

# Verify no negative values
if (any(ref_counts_common@x < 0)) {
  warning("Negative values found in reference counts. Setting to 0.")
  ref_counts_common@x[ref_counts_common@x < 0] <- 0
}

# Create Reference object
reference <- Reference(
  counts = ref_counts_common,
  cell_types = reference_data$cell_types,
  nUMI = reference_data$nUMI
)

message("Reference object created successfully")
message(paste("  Cell types:", paste(levels(reference@cell_types), collapse = ", ")))

# -----------------------------------------------------------------------------
# Create RCTD SpatialRNA object
# -----------------------------------------------------------------------------
message("\n=== Creating SpatialRNA object ===")

# Ensure ST counts are in correct format
st_counts_common <- as(st_counts_common, "dgCMatrix")

# Verify no negative values
if (any(st_counts_common@x < 0)) {
  warning("Negative values found in ST counts. Setting to 0.")
  st_counts_common@x[st_counts_common@x < 0] <- 0
}

# Ensure coords are data.frame with numeric columns
coords_df <- data.frame(
  x = as.numeric(coords[, 1]),
  y = as.numeric(coords[, 2]),
  row.names = rownames(coords)
)

# Ensure spot names match
common_spots <- intersect(colnames(st_counts_common), rownames(coords_df))
st_counts_common <- st_counts_common[, common_spots]
coords_df <- coords_df[common_spots, , drop = FALSE]

# Create SpatialRNA object
puck <- SpatialRNA(
  coords = coords_df,
  counts = st_counts_common,
  nUMI = colSums(st_counts_common)
)

message("SpatialRNA object created successfully")
message(paste("  Spots:", ncol(puck@counts)))
message(paste("  Genes:", nrow(puck@counts)))

# -----------------------------------------------------------------------------
# Run RCTD
# -----------------------------------------------------------------------------
message("\n=== Running RCTD deconvolution ===")
message(paste("Mode:", DOUBLET_MODE))
message(paste("Cores:", N_CORES))

# Create RCTD object
myRCTD <- create.RCTD(
  spatialRNA = puck,
  reference = reference,
  max_cores = N_CORES,
  CELL_MIN_INSTANCE = 10  # Minimum cells per type in reference
)

message("\nStarting RCTD fitting...")
start_time <- Sys.time()

# Run RCTD
myRCTD <- run.RCTD(
  myRCTD,
  doublet_mode = DOUBLET_MODE
)

end_time <- Sys.time()
message(paste("\nRCTD completed in:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))

# -----------------------------------------------------------------------------
# Extract Results
# -----------------------------------------------------------------------------
message("\n=== Extracting results ===")

# Get results based on mode
if (DOUBLET_MODE == "full") {
  # Full mode: get cell type weights for each spot
  weights <- myRCTD@results$weights

  # Normalize weights to sum to 1
  weights_normalized <- normalize_weights(weights)

  message("Extracted cell type weights (full mode)")
  message(paste("  Spots:", nrow(weights_normalized)))
  message(paste("  Cell types:", ncol(weights_normalized)))

} else if (DOUBLET_MODE == "doublet") {
  # Doublet mode: get spot classifications
  results_df <- myRCTD@results$results_df

  # Extract first and second cell types with weights
  weights_normalized <- data.frame(
    spot = rownames(results_df),
    first_type = results_df$first_type,
    second_type = results_df$second_type,
    first_weight = 1 - results_df$second_type_weight,
    second_weight = results_df$second_type_weight,
    spot_class = results_df$spot_class
  )

  message("Extracted cell type classifications (doublet mode)")

} else {
  # Multi mode
  weights <- myRCTD@results$weights
  weights_normalized <- normalize_weights(weights)
  message("Extracted cell type weights (multi mode)")
}

# Print summary statistics
message("\n=== Results Summary ===")

if (DOUBLET_MODE == "full") {
  message("\nMean cell type proportions across all spots:")
  mean_props <- colMeans(weights_normalized)
  mean_props <- sort(mean_props, decreasing = TRUE)
  for (ct in names(mean_props)) {
    message(sprintf("  %s: %.1f%%", ct, mean_props[ct] * 100))
  }

  # Check that weights sum to ~1
  weight_sums <- rowSums(weights_normalized)
  message(paste("\nWeight sum range:", round(min(weight_sums), 3), "-", round(max(weight_sums), 3)))
}

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------
message("\n=== Saving results ===")

# Save full RCTD object
rctd_output <- file.path(OUTPUT_DIR, "rctd_results.rds")
saveRDS(myRCTD, rctd_output)
message(paste("Saved RCTD object to:", rctd_output))

# Save normalized weights
weights_output <- file.path(OUTPUT_DIR, "cell_type_weights.csv")
write.csv(weights_normalized, weights_output)
message(paste("Saved cell type weights to:", weights_output))

# Save weights as RDS for easier R loading
weights_rds <- file.path(OUTPUT_DIR, "cell_type_weights.rds")
saveRDS(weights_normalized, weights_rds)
message(paste("Saved cell type weights (RDS) to:", weights_rds))

# -----------------------------------------------------------------------------
# Add results to Seurat object
# -----------------------------------------------------------------------------
if (file.exists(st_seurat_file)) {
  message("\n=== Adding results to Seurat object ===")

  st_data <- readRDS(st_seurat_file)

  if (DOUBLET_MODE == "full") {
    # Find common spots between Seurat and weights
    common_spots_seurat <- intersect(colnames(st_data), rownames(weights_normalized))
    message(paste("Adding results for", length(common_spots_seurat), "spots"))

    # Add each cell type as metadata (only for spots with results)
    for (ct in colnames(weights_normalized)) {
      ct_clean <- gsub("[^a-zA-Z0-9]", "_", ct)  # Clean column name
      # Initialize with NA
      ct_values <- rep(NA, ncol(st_data))
      names(ct_values) <- colnames(st_data)
      # Fill in values for spots with results
      ct_values[common_spots_seurat] <- weights_normalized[common_spots_seurat, ct]
      st_data[[ct_clean]] <- ct_values
    }
  }

  # Save updated Seurat object
  seurat_output <- file.path(OUTPUT_DIR, "st_data_deconvolved.rds")
  saveRDS(st_data, seurat_output)
  message(paste("Saved deconvolved Seurat object to:", seurat_output))
}

message("\n=== RCTD deconvolution complete! ===")
message("Next step: Run 05_visualize.R")

# Return weights for interactive use
weights_normalized
