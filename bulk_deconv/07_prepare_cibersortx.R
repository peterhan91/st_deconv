# =============================================================================
# 07_prepare_cibersortx.R - Prepare files for CIBERSORTx
# Can be used for web upload OR Docker
# =============================================================================

library(Biobase)
library(dplyr)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
OUTPUT_DIR <- "output"
CIBERSORT_DIR <- file.path(OUTPUT_DIR, "cibersortx")
if (!dir.exists(CIBERSORT_DIR)) dir.create(CIBERSORT_DIR, recursive = TRUE)

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
message("=== Loading data ===")

bulk_eset <- readRDS(file.path(OUTPUT_DIR, "data", "bulk_eset.rds"))
sc_eset <- readRDS(file.path(OUTPUT_DIR, "data", "sc_reference_eset_filtered.rds"))

bulk_mtx <- exprs(bulk_eset)
sc_mtx <- exprs(sc_eset)
sc_meta <- pData(sc_eset)

message(paste("Bulk:", nrow(bulk_mtx), "genes x", ncol(bulk_mtx), "samples"))
message(paste("Reference:", nrow(sc_mtx), "genes x", ncol(sc_mtx), "cells"))

# -----------------------------------------------------------------------------
# Create Signature Matrix for CIBERSORTx
# -----------------------------------------------------------------------------
message("\n=== Creating signature matrix ===")

# CIBERSORTx needs a signature matrix (genes x cell types)
# Calculate mean expression per cell type
cell_types <- unique(sc_meta$cellType)

signature_matrix <- sapply(cell_types, function(ct) {
  cells <- which(sc_meta$cellType == ct)
  rowMeans(sc_mtx[, cells, drop = FALSE])
})

# Add gene names as first column
sig_df <- data.frame(
  GeneSymbol = rownames(signature_matrix),
  signature_matrix,
  check.names = FALSE
)

# Save signature matrix
sig_file <- file.path(CIBERSORT_DIR, "signature_matrix.txt")
write.table(sig_df, sig_file, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste("Saved:", sig_file))
message(paste("Dimensions:", nrow(sig_df), "genes x", length(cell_types), "cell types"))

# -----------------------------------------------------------------------------
# Create Mixture File (Bulk Expression)
# -----------------------------------------------------------------------------
message("\n=== Creating mixture file ===")

# CIBERSORTx needs mixture file (genes x samples)
mixture_df <- data.frame(
  GeneSymbol = rownames(bulk_mtx),
  bulk_mtx,
  check.names = FALSE
)

# Save mixture file
mix_file <- file.path(CIBERSORT_DIR, "mixture_file.txt")
write.table(mixture_df, mix_file, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste("Saved:", mix_file))
message(paste("Dimensions:", nrow(mixture_df), "genes x", ncol(bulk_mtx), "samples"))

# -----------------------------------------------------------------------------
# Create Single-Cell Reference (for CIBERSORTx scRNA-seq mode)
# -----------------------------------------------------------------------------
message("\n=== Creating single-cell reference files ===")

# Sample 500 cells per type for manageable size
set.seed(42)
max_cells <- 500
sampled_cells <- c()

for (ct in cell_types) {
  ct_cells <- which(sc_meta$cellType == ct)
  n_sample <- min(length(ct_cells), max_cells)
  sampled_cells <- c(sampled_cells, sample(ct_cells, n_sample))
}

sc_mtx_sub <- sc_mtx[, sampled_cells]
sc_meta_sub <- sc_meta[sampled_cells, ]

# Single-cell reference matrix
sc_ref_df <- data.frame(
  GeneSymbol = rownames(sc_mtx_sub),
  sc_mtx_sub,
  check.names = FALSE
)

sc_ref_file <- file.path(CIBERSORT_DIR, "sc_reference_matrix.txt")
write.table(sc_ref_df, sc_ref_file, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste("Saved:", sc_ref_file))

# Phenotype/label file for single cells
phenotype_df <- data.frame(
  SampleID = colnames(sc_mtx_sub),
  CellType = sc_meta_sub$cellType
)

pheno_file <- file.path(CIBERSORT_DIR, "sc_phenotypes.txt")
write.table(phenotype_df, pheno_file, sep = "\t", quote = FALSE, row.names = FALSE)
message(paste("Saved:", pheno_file))

# -----------------------------------------------------------------------------
# Instructions
# -----------------------------------------------------------------------------
message("\n", strrep("=", 60))
message("FILES READY FOR CIBERSORTx")
message(strrep("=", 60))

message("
Output files in: ", CIBERSORT_DIR, "

FILES CREATED:
  1. signature_matrix.txt    - Pre-computed signature (genes x cell types)
  2. mixture_file.txt        - Bulk expression to deconvolve
  3. sc_reference_matrix.txt - Single-cell reference (for S-mode)
  4. sc_phenotypes.txt       - Cell type labels for sc reference

=== OPTION A: CIBERSORTx WEB INTERFACE ===
1. Go to: https://cibersortx.stanford.edu/
2. Register/Login
3. Click 'Run CIBERSORTx Fractions'
4. Upload:
   - Signature matrix: signature_matrix.txt
   - Mixture file: mixture_file.txt
5. Run and download results

=== OPTION B: CIBERSORTx DOCKER ===
# Install Docker Desktop first, then:

# 1. Get token from CIBERSORTx website after registration

# 2. Pull image
docker pull cibersortx/fractions

# 3. Run (replace YOUR_EMAIL and YOUR_TOKEN)
docker run -v ", normalizePath(CIBERSORT_DIR), ":/src/data \\
  -v ", normalizePath(CIBERSORT_DIR), ":/src/outdir \\
  cibersortx/fractions \\
  --username YOUR_EMAIL \\
  --token YOUR_TOKEN \\
  --single_cell TRUE \\
  --refsample /src/data/sc_reference_matrix.txt \\
  --phenoclasses /src/data/sc_phenotypes.txt \\
  --mixture /src/data/mixture_file.txt \\
  --perm 100 \\
  --verbose TRUE

# Results will be saved to: ", CIBERSORT_DIR, "/CIBERSORTx_Results.txt
")

message("\n=== Preparation complete ===")
