# =============================================================================
# 02_download_cptac.R - Download CPTAC BRCA RNA-seq data
# Bulk Deconvolution Pipeline using MuSiC
# =============================================================================

# Load GenomicDataCommons
library(GenomicDataCommons)
library(SummarizedExperiment)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
DATA_DIR <- "data"
if (!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive = TRUE)

# -----------------------------------------------------------------------------
# Check GDC API status
# -----------------------------------------------------------------------------
message("=== Checking GDC API status ===")
status <- GenomicDataCommons::status()
message(paste("GDC API Status:", status$status))

# -----------------------------------------------------------------------------
# Query CPTAC-3 BRCA RNA-seq data
# -----------------------------------------------------------------------------
message("\n=== Querying CPTAC BRCA RNA-seq data ===")

# Use explicit GenomicDataCommons namespace to avoid dplyr conflicts
# Query for gene expression quantification files from CPTAC-3
query <- GenomicDataCommons::files() |>
  GenomicDataCommons::filter(
    cases.project.project_id == "CPTAC-3" &
    data_category == "Transcriptome Profiling" &
    data_type == "Gene Expression Quantification"
  )

# Get results (returns a list structure)
results <- GenomicDataCommons::results(query, size = 500)

# Check how many files found
n_files_found <- length(results$file_id)
message(paste("Found", n_files_found, "CPTAC-3 RNA-seq files"))

if (n_files_found == 0) {
  message("\nNo CPTAC-3 data found with initial query. Trying broader query...")

  # Try broader CPTAC query
  query2 <- GenomicDataCommons::files() |>
    GenomicDataCommons::filter(
      cases.project.program.name == "CPTAC" &
      data_category == "Transcriptome Profiling"
    )

  results2 <- GenomicDataCommons::results(query2, size = 100)

  message(paste("\nFound", length(results2$file_id), "CPTAC files total"))
}

# -----------------------------------------------------------------------------
# Alternative: Use pre-processed CPTAC data from cBioPortal
# -----------------------------------------------------------------------------
message("\n=== Alternative: Download from cBioPortal ===")
message("
CPTAC BRCA data is also available from cBioPortal:
  https://www.cbioportal.org/study/summary?id=brca_cptac_2020

Manual download:
  1. Go to: https://www.cbioportal.org/study/summary?id=brca_cptac_2020
  2. Click 'Download' -> 'All Data'
  3. Extract and place RNA-seq data in bulk_deconv/data/
")

# -----------------------------------------------------------------------------
# Download files if available
# -----------------------------------------------------------------------------
if (n_files_found > 0) {
  message("\n=== Downloading sample files ===")

  # Download first 20 files for testing
  n_download <- min(20, n_files_found)
  file_ids <- results$file_id[1:n_download]
  file_names <- results$file_name[1:n_download]

  message(paste("Downloading", n_download, "files..."))

  # Create manifest
  manifest <- data.frame(
    file_id = file_ids,
    file_name = file_names,
    stringsAsFactors = FALSE
  )
  write.csv(manifest, file.path(DATA_DIR, "cptac_manifest.csv"), row.names = FALSE)
  message("Saved manifest: data/cptac_manifest.csv")

  # Download files
  download_dir <- file.path(DATA_DIR, "cptac_raw")
  if (!dir.exists(download_dir)) dir.create(download_dir)

  tryCatch({
    message("Starting download from GDC...")
    downloaded_files <- GenomicDataCommons::gdcdata(
      file_ids,
      destination_dir = download_dir,
      progress = TRUE
    )
    message("Download complete!")
    message(paste("Files saved to:", download_dir))
  }, error = function(e) {
    message(paste("Download error:", e$message))
    message("Please download manually from GDC portal or cBioPortal")
  })
}

# -----------------------------------------------------------------------------
# Instructions for real data
# -----------------------------------------------------------------------------
message("\n", strrep("=", 60))
message("INSTRUCTIONS FOR REAL CPTAC DATA")
message(strrep("=", 60))
message("
Option 1: GDC Data Portal (Recommended)
  1. Go to: https://portal.gdc.cancer.gov/
  2. Repository -> Cases -> Project: CPTAC-3
  3. Filter: Primary Site = Breast
  4. Files -> Data Category = Transcriptome Profiling
  5. Add to Cart -> Download -> Manifest or Cart

Option 2: cBioPortal (Easiest)
  1. Go to: https://www.cbioportal.org/study/summary?id=brca_cptac_2020
  2. Click 'Download' -> 'All Data'
  3. Extract brca_cptac_2020.tar.gz
  4. Copy data_mrna_seq_v2_rsem.txt to bulk_deconv/data/

Option 3: CPTAC Data Portal
  1. Go to: https://proteomic.datacommons.cancer.gov/pdc/
  2. Search: CPTAC BRCA
  3. Download RNA-seq data

After download, place the expression matrix in:
  bulk_deconv/data/cptac_brca_expression.rds (or .csv/.txt)
  OR
  bulk_deconv/data/data_mrna_seq_v2_rsem.txt (from cBioPortal)

Then run: 03_prepare_bulk_data.R
")

message("\n=== Download script complete ===")
