# =============================================================================
# 01_setup.R - Package Installation and Environment Setup
# ST Deconvolution Pipeline for IDC Data using RCTD
# =============================================================================

# Function to check and install packages
install_if_missing <- function(packages, source = "CRAN") {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing", pkg, "from", source, "..."))
      if (source == "CRAN") {
        install.packages(pkg)
      } else if (source == "Bioconductor") {
        BiocManager::install(pkg)
      }
    } else {
      message(paste(pkg, "is already installed."))
    }
  }
}

# -----------------------------------------------------------------------------
# Step 1: Install CRAN packages
# -----------------------------------------------------------------------------
message("=== Installing CRAN packages ===")

cran_packages <- c(
  "Matrix",       # Sparse matrix operations
  "Seurat",       # Single-cell/spatial analysis
  "ggplot2",      # Plotting
  "patchwork",    # Combining plots
  "dplyr",        # Data manipulation
  "tidyr",        # Data tidying
  "RColorBrewer", # Color palettes
  "viridis",      # Color scales
  "devtools",     # For GitHub installations
  "remotes"       # Alternative for GitHub installations
)

install_if_missing(cran_packages, "CRAN")

# -----------------------------------------------------------------------------
# Step 2: Install BiocManager if needed
# -----------------------------------------------------------------------------
message("\n=== Setting up Bioconductor ===")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# -----------------------------------------------------------------------------
# Step 3: Install Bioconductor packages
# -----------------------------------------------------------------------------
message("\n=== Installing Bioconductor packages ===")

bioc_packages <- c(
  "zellkonverter",          # Read h5ad files
  "SingleCellExperiment",   # SCE object class
  "SummarizedExperiment",   # Required for SCE
  "basilisk"                # Python environment management for zellkonverter
)

install_if_missing(bioc_packages, "Bioconductor")

# -----------------------------------------------------------------------------
# Step 4: Install spacexr (RCTD)
# -----------------------------------------------------------------------------
message("\n=== Installing spacexr (RCTD) ===")

if (!requireNamespace("spacexr", quietly = TRUE)) {
  message("Installing spacexr from GitHub...")
  devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
} else {
  message("spacexr is already installed.")
}

# -----------------------------------------------------------------------------
# Step 5: Verify installations
# -----------------------------------------------------------------------------
message("\n=== Verifying installations ===")

all_packages <- c(cran_packages, bioc_packages, "spacexr")

check_results <- sapply(all_packages, function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
})

if (all(check_results)) {
  message("\nAll packages installed successfully!")
} else {
  failed <- names(check_results)[!check_results]
  warning(paste("\nFailed to install:", paste(failed, collapse = ", ")))
}

# Print package versions
message("\n=== Package versions ===")
for (pkg in all_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    version <- as.character(packageVersion(pkg))
    message(paste(pkg, ":", version))
  }
}

# -----------------------------------------------------------------------------
# Step 6: Create output directories
# -----------------------------------------------------------------------------
message("\n=== Creating output directories ===")

dirs <- c("output", "output/figures", "output/data", "reference")
for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    message(paste("Created directory:", d))
  }
}

message("\n=== Setup complete! ===")
message("Next step: Run 02_load_data.R")
