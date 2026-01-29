# =============================================================================
# 08_run_bayesprism.R - Run BayesPrism deconvolution
# =============================================================================

library(BayesPrism)
library(Biobase)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
OUTPUT_DIR <- "output"
DATA_DIR <- file.path(OUTPUT_DIR, "data")

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
message("=== Loading data ===")

bulk_eset <- readRDS(file.path(DATA_DIR, "bulk_eset.rds"))
sc_eset <- readRDS(file.path(DATA_DIR, "sc_reference_eset_filtered.rds"))

bulk_mtx <- exprs(bulk_eset)
sc_mtx <- exprs(sc_eset)
sc_meta <- pData(sc_eset)

message(paste("Bulk:", nrow(bulk_mtx), "genes x", ncol(bulk_mtx), "samples"))
message(paste("Reference:", nrow(sc_mtx), "genes x", ncol(sc_mtx), "cells"))

# -----------------------------------------------------------------------------
# Prepare data for BayesPrism
# -----------------------------------------------------------------------------
message("\n=== Preparing data for BayesPrism ===")

# BayesPrism expects:
# - Reference: genes x cells matrix (raw counts)
# - Bulk: genes x samples matrix (raw counts)
# - Cell type labels as a character vector

# Get common genes
common_genes <- intersect(rownames(bulk_mtx), rownames(sc_mtx))
message(paste("Common genes:", length(common_genes)))

# Subset to common genes
bulk_sub <- bulk_mtx[common_genes, ]
sc_sub <- sc_mtx[common_genes, ]

# BayesPrism needs cells x genes format for reference
# And samples x genes for bulk
sc_t <- t(sc_sub)  # cells x genes
bulk_t <- t(bulk_sub)  # samples x genes

# Cell type labels
cell_types <- as.character(sc_meta$cellType)

message(paste("Reference matrix:", nrow(sc_t), "cells x", ncol(sc_t), "genes"))
message(paste("Bulk matrix:", nrow(bulk_t), "samples x", ncol(bulk_t), "genes"))
message(paste("Cell types:", length(unique(cell_types))))

# -----------------------------------------------------------------------------
# Create BayesPrism reference
# -----------------------------------------------------------------------------
message("\n=== Creating BayesPrism reference ===")

# BayesPrism reference object
# Note: BayesPrism expects raw counts, not normalized data
# If data is normalized, we may need to round or use as-is

# Check if data looks like counts
if (max(sc_t) < 100) {
  message("Data appears to be normalized - converting to pseudo-counts")
  sc_t <- round(sc_t * 100)  # Scale up for pseudo-counts
  bulk_t <- round(bulk_t * 100)
}

# Filter genes with low expression
gene_sums <- colSums(sc_t)
keep_genes <- gene_sums > ncol(sc_t) * 0.01  # Keep genes expressed in >1% cells
message(paste("Genes after filtering:", sum(keep_genes)))

sc_t_filt <- sc_t[, keep_genes]
bulk_t_filt <- bulk_t[, keep_genes]

# Create prism reference
message("Creating prism reference...")
myPrism <- new.prism(
  reference = sc_t_filt,
  mixture = bulk_t_filt,
  input.type = "count.matrix",
  cell.type.labels = cell_types,
  cell.state.labels = cell_types,  # Use same as cell types for simplicity
  key = NULL,
  outlier.cut = 0.01,
  outlier.fraction = 0.1
)

# -----------------------------------------------------------------------------
# Run BayesPrism
# -----------------------------------------------------------------------------
message("\n=== Running BayesPrism deconvolution ===")
message("This may take a while...")

# Run with default settings
bp_result <- run.prism(
  prism = myPrism,
  n.cores = 2
)

# -----------------------------------------------------------------------------
# Extract Results
# -----------------------------------------------------------------------------
message("\n=== Extracting results ===")

# Get cell type fractions
theta <- get.fraction(
  bp = bp_result,
  which.theta = "final",
  state.or.type = "type"
)

# theta is samples x cell types
props_bp <- as.data.frame(theta)
props_bp$Sample <- rownames(props_bp)

message("BayesPrism proportions (first 5 samples):")
print(head(props_bp[, 1:min(5, ncol(props_bp)-1)]))

# Calculate mean proportions
mean_props <- colMeans(theta) * 100
message("\nMean proportions across all samples:")
for (ct in names(sort(mean_props, decreasing = TRUE))) {
  message(sprintf("  %s: %.1f%%", ct, mean_props[ct]))
}

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------
message("\n=== Saving results ===")

saveRDS(bp_result, file.path(DATA_DIR, "bayesprism_result.rds"))
saveRDS(props_bp, file.path(DATA_DIR, "props_bayesprism.rds"))

# Save as CSV for easy viewing
write.csv(props_bp, file.path(DATA_DIR, "props_bayesprism.csv"), row.names = FALSE)

message("Results saved to:")
message(paste("  -", file.path(DATA_DIR, "bayesprism_result.rds")))
message(paste("  -", file.path(DATA_DIR, "props_bayesprism.rds")))
message(paste("  -", file.path(DATA_DIR, "props_bayesprism.csv")))

message("\n=== BayesPrism deconvolution complete ===")
