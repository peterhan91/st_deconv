# =============================================================================
# 01_setup_bulk.R - Install packages for bulk RNA-seq deconvolution
# Bulk Deconvolution Pipeline using MuSiC
# =============================================================================

message("=== Installing packages for bulk RNA-seq deconvolution ===")

# -----------------------------------------------------------------------------
# CRAN packages
# -----------------------------------------------------------------------------
cran_packages <- c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "tidyr",
  "pheatmap",
  "RColorBrewer",
  "devtools"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing", pkg, "..."))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    message(paste(pkg, "already installed"))
  }
}

# -----------------------------------------------------------------------------
# Bioconductor packages
# -----------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "Biobase",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "GenomicDataCommons",    # For downloading GDC/CPTAC data
  "TCGAbiolinks",          # Alternative for TCGA/GDC data
  "ExperimentHub",
  "scater"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing", pkg, "from Bioconductor..."))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    message(paste(pkg, "already installed"))
  }
}

# -----------------------------------------------------------------------------
# MuSiC - from GitHub
# -----------------------------------------------------------------------------
message("\n=== Installing MuSiC ===")

if (!requireNamespace("MuSiC", quietly = TRUE)) {
  message("Installing MuSiC from GitHub...")
  devtools::install_github("xuranw/MuSiC", build_vignettes = FALSE)
} else {
  message("MuSiC already installed")
}

# -----------------------------------------------------------------------------
# Verify installations
# -----------------------------------------------------------------------------
message("\n=== Verifying installations ===")

all_packages <- c(cran_packages, bioc_packages, "MuSiC")

check_results <- sapply(all_packages, function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
})

if (all(check_results)) {
  message("\nAll packages installed successfully!")
} else {
  failed <- names(check_results)[!check_results]
  warning(paste("\nFailed to install:", paste(failed, collapse = ", ")))
}

message("\n=== Setup complete ===")
message("Next: Run 02_download_cptac.R to download CPTAC BRCA data")
