# =============================================================================
# 04_run_music.R - Run MuSiC Deconvolution
# Bulk Deconvolution Pipeline using MuSiC
# =============================================================================

library(MuSiC)
library(Biobase)
library(SingleCellExperiment)
library(dplyr)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
OUTPUT_DIR <- "output"

# -----------------------------------------------------------------------------
# Load Prepared Data
# -----------------------------------------------------------------------------
message("=== Loading prepared data ===")

# Load ExpressionSets
bulk_eset <- readRDS(file.path(OUTPUT_DIR, "data", "bulk_eset.rds"))
sc_eset <- readRDS(file.path(OUTPUT_DIR, "data", "sc_reference_eset_filtered.rds"))

message(paste("Bulk data:", nrow(bulk_eset), "genes x", ncol(bulk_eset), "samples"))
message(paste("Reference:", nrow(sc_eset), "genes x", ncol(sc_eset), "cells"))

# Check cell types
cell_types <- unique(pData(sc_eset)$cellType)
message(paste("\nCell types:", paste(cell_types, collapse = ", ")))

# -----------------------------------------------------------------------------
# Convert to formats required by new MuSiC API
# -----------------------------------------------------------------------------
message("\n=== Preparing data for MuSiC ===")

# Extract bulk matrix
bulk_mtx <- exprs(bulk_eset)
message(paste("Bulk matrix:", nrow(bulk_mtx), "x", ncol(bulk_mtx)))

# Convert ExpressionSet to SingleCellExperiment for new MuSiC API
sc_counts <- exprs(sc_eset)
sc_coldata <- pData(sc_eset)

sc_sce <- SingleCellExperiment(
  assays = list(counts = sc_counts),
  colData = sc_coldata
)
message(paste("SingleCellExperiment:", nrow(sc_sce), "x", ncol(sc_sce)))

# -----------------------------------------------------------------------------
# Run MuSiC Deconvolution
# -----------------------------------------------------------------------------
message("\n=== Running MuSiC deconvolution ===")
message("This may take several minutes...")

# Run MuSiC with new API
start_time <- Sys.time()

music_result <- music_prop(
  bulk.mtx = bulk_mtx,
  sc.sce = sc_sce,
  clusters = "cellType",
  samples = "SubjectName",
  verbose = TRUE
)

end_time <- Sys.time()
message(paste("\nDeconvolution completed in:", round(difftime(end_time, start_time, units = "mins"), 2), "minutes"))

# -----------------------------------------------------------------------------
# Extract and Process Results
# -----------------------------------------------------------------------------
message("\n=== Processing results ===")

# MuSiC returns estimated proportions in Est.prop.weighted
proportions <- music_result$Est.prop.weighted

message(paste("Result dimensions:", nrow(proportions), "samples x", ncol(proportions), "cell types"))

# Check that proportions sum to ~1
prop_sums <- rowSums(proportions)
message(paste("\nProportion sums - Mean:", round(mean(prop_sums), 4),
              "SD:", round(sd(prop_sums), 4)))

# -----------------------------------------------------------------------------
# Summary Statistics
# -----------------------------------------------------------------------------
message("\n=== Cell Type Proportions Summary ===")

mean_props <- colMeans(proportions) * 100
sd_props <- apply(proportions, 2, sd) * 100

summary_df <- data.frame(
  CellType = names(mean_props),
  Mean_Percent = round(mean_props, 2),
  SD_Percent = round(sd_props, 2)
) %>%
  arrange(desc(Mean_Percent))

print(summary_df)

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------
message("\n=== Saving results ===")

# Save proportions as matrix
saveRDS(proportions, file.path(OUTPUT_DIR, "data", "music_proportions.rds"))
message("Saved: music_proportions.rds")

# Save as CSV for easy viewing
write.csv(proportions, file.path(OUTPUT_DIR, "data", "music_proportions.csv"))
message("Saved: music_proportions.csv")

# Save full MuSiC result object
saveRDS(music_result, file.path(OUTPUT_DIR, "data", "music_result_full.rds"))
message("Saved: music_result_full.rds")

# Save summary statistics
write.csv(summary_df, file.path(OUTPUT_DIR, "data", "deconvolution_summary.csv"), row.names = FALSE)
message("Saved: deconvolution_summary.csv")

# -----------------------------------------------------------------------------
# Additional MuSiC Outputs
# -----------------------------------------------------------------------------
message("\n=== Additional MuSiC information ===")

# MuSiC also provides:
# - Est.prop.allgene: Proportions using all genes (no weighting)
# - Weight.gene: Gene-specific weights used in estimation
# - r.squared.full: R-squared for each sample fit

if (!is.null(music_result$r.squared.full)) {
  r2_values <- music_result$r.squared.full
  message(paste("Model fit (R-squared) - Mean:", round(mean(r2_values), 3),
                "Range:", round(min(r2_values), 3), "-", round(max(r2_values), 3)))

  # Save R-squared values
  r2_df <- data.frame(
    Sample = names(r2_values),
    R_squared = r2_values
  )
  write.csv(r2_df, file.path(OUTPUT_DIR, "data", "model_fit_rsquared.csv"), row.names = FALSE)
}

# Compare weighted vs unweighted estimates
if (!is.null(music_result$Est.prop.allgene)) {
  prop_allgene <- music_result$Est.prop.allgene

  # Correlation between methods
  cor_methods <- diag(cor(proportions, prop_allgene))
  message("\nCorrelation between weighted and unweighted estimates:")
  print(round(cor_methods, 3))

  # Save unweighted proportions
  write.csv(prop_allgene, file.path(OUTPUT_DIR, "data", "music_proportions_unweighted.csv"))
  message("Saved: music_proportions_unweighted.csv")
}

# -----------------------------------------------------------------------------
# Quality Checks
# -----------------------------------------------------------------------------
message("\n", strrep("=", 60))
message("QUALITY CHECKS")
message(strrep("=", 60))

# Check 1: Any samples with failed deconvolution?
failed_samples <- which(prop_sums < 0.9 | prop_sums > 1.1)
if (length(failed_samples) > 0) {
  message(paste("WARNING:", length(failed_samples), "samples have unusual proportion sums"))
} else {
  message("PASS: All samples have valid proportion sums")
}

# Check 2: Any cell type completely absent?
absent_types <- names(mean_props)[mean_props < 0.1]
if (length(absent_types) > 0) {
  message(paste("NOTE: Low abundance cell types (<0.1%):", paste(absent_types, collapse = ", ")))
}

# Check 3: Dominant cell types
dominant_type <- names(which.max(mean_props))
message(paste("\nMost abundant cell type:", dominant_type, "(", round(max(mean_props), 1), "%)"))

# Check 4: Expected patterns for breast cancer
message("\nBiological plausibility check:")
if ("CAFs" %in% names(mean_props) && mean_props["CAFs"] > 10) {
  message("  - CAFs abundant: Expected for breast cancer")
}
if ("Cancer Epithelial" %in% names(mean_props) && mean_props["Cancer Epithelial"] > 5) {
  message("  - Cancer cells detected: Expected")
}
immune_total <- sum(mean_props[names(mean_props) %in% c("T-cells", "B-cells", "Myeloid")])
message(paste("  - Total immune infiltration:", round(immune_total, 1), "%"))

# -----------------------------------------------------------------------------
# Final Summary
# -----------------------------------------------------------------------------
message("\n", strrep("=", 60))
message("DECONVOLUTION COMPLETE")
message(strrep("=", 60))
message(paste("\nProcessed", nrow(proportions), "bulk samples"))
message(paste("Estimated proportions for", ncol(proportions), "cell types"))
message("\nOutput files in output/data/:")
message("  - music_proportions.csv (main results)")
message("  - music_proportions.rds (R format)")
message("  - deconvolution_summary.csv (statistics)")
message("  - model_fit_rsquared.csv (quality metrics)")

message("\n=== Next: Run 05_visualize_bulk.R ===")
