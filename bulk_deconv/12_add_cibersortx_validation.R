# =============================================================================
# 12_add_cibersortx_validation.R - Add CIBERSORTx to Pseudo-bulk Validation
# =============================================================================
# Run this AFTER you've run CIBERSORTx on pseudobulk_for_cibersortx.txt
# and downloaded the results.

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(patchwork)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
VALIDATION_DIR <- "output/validation"

# UPDATE THIS PATH to your downloaded CIBERSORTx results file
CIBERSORTX_RESULTS <- file.path(VALIDATION_DIR, "CIBERSORTx_Job3_Results.txt")

# -----------------------------------------------------------------------------
# Load Previous Results
# -----------------------------------------------------------------------------
message("=== Loading previous validation results ===")

ground_truth <- readRDS(file.path(VALIDATION_DIR, "ground_truth_proportions.rds"))
results_list <- readRDS(file.path(VALIDATION_DIR, "estimated_proportions_all.rds"))
prev_summary <- read.csv(file.path(VALIDATION_DIR, "method_performance_summary.csv"))

message("Previous methods:", paste(names(results_list), collapse = ", "))
message("Ground truth samples:", nrow(ground_truth))

# -----------------------------------------------------------------------------
# Load CIBERSORTx Results
# -----------------------------------------------------------------------------
message("\n=== Loading CIBERSORTx results ===")

if (!file.exists(CIBERSORTX_RESULTS)) {
  stop(paste("CIBERSORTx results not found at:", CIBERSORTX_RESULTS,
             "\nPlease run CIBERSORTx on pseudobulk_for_cibersortx.txt and save results here."))
}

cibersortx_raw <- read.delim(CIBERSORTX_RESULTS, stringsAsFactors = FALSE)
message(paste("Loaded", nrow(cibersortx_raw), "samples"))

# Extract proportions (exclude Mixture, P-value, Correlation, RMSE columns)
cibersortx_samples <- cibersortx_raw$Mixture
cibersortx_celltypes <- colnames(cibersortx_raw)[2:(ncol(cibersortx_raw) - 3)]
cibersortx_props <- as.matrix(cibersortx_raw[, cibersortx_celltypes])
rownames(cibersortx_props) <- cibersortx_samples

# Fix cell type names (dots to dashes/spaces)
colnames(cibersortx_props) <- gsub("\\.", " ", colnames(cibersortx_props))
colnames(cibersortx_props) <- gsub("B cells", "B-cells", colnames(cibersortx_props))
colnames(cibersortx_props) <- gsub("T cells", "T-cells", colnames(cibersortx_props))

message("CIBERSORTx cell types:", paste(colnames(cibersortx_props), collapse = ", "))

# CIBERSORTx QC
message("\nCIBERSORTx QC:")
message(sprintf("  Mean P-value: %.4f", mean(cibersortx_raw$P.value)))
message(sprintf("  Mean Correlation: %.4f", mean(cibersortx_raw$Correlation)))
message(sprintf("  Mean RMSE: %.4f", mean(cibersortx_raw$RMSE)))

# -----------------------------------------------------------------------------
# Evaluate CIBERSORTx
# -----------------------------------------------------------------------------
message("\n=== Evaluating CIBERSORTx ===")

evaluate_method <- function(estimated, ground_truth) {
  common_samples <- intersect(rownames(estimated), rownames(ground_truth))
  common_types <- intersect(colnames(estimated), colnames(ground_truth))

  est <- estimated[common_samples, common_types]
  truth <- ground_truth[common_samples, common_types]

  overall_cor <- cor(as.vector(est), as.vector(truth))
  overall_rmse <- sqrt(mean((as.vector(est) - as.vector(truth))^2))
  overall_mae <- mean(abs(as.vector(est) - as.vector(truth)))

  celltype_cors <- sapply(common_types, function(ct) cor(est[, ct], truth[, ct]))
  celltype_rmse <- sapply(common_types, function(ct) sqrt(mean((est[, ct] - truth[, ct])^2)))

  return(list(
    overall_cor = overall_cor,
    overall_rmse = overall_rmse,
    overall_mae = overall_mae,
    celltype_cors = celltype_cors,
    celltype_rmse = celltype_rmse
  ))
}

cibersortx_eval <- evaluate_method(cibersortx_props, ground_truth)

message(sprintf("\nCIBERSORTx Performance:"))
message(sprintf("  Correlation: %.4f", cibersortx_eval$overall_cor))
message(sprintf("  RMSE: %.4f", cibersortx_eval$overall_rmse))
message(sprintf("  MAE: %.4f", cibersortx_eval$overall_mae))

# Add to results list
results_list$CIBERSORTx <- cibersortx_props

# -----------------------------------------------------------------------------
# Updated Summary
# -----------------------------------------------------------------------------
message("\n=== Updated Method Comparison ===")

# Re-evaluate all methods
all_evals <- lapply(results_list, function(est) evaluate_method(est, ground_truth))

summary_table <- data.frame(
  Method = names(all_evals),
  Correlation = sapply(all_evals, function(x) round(x$overall_cor, 4)),
  RMSE = sapply(all_evals, function(x) round(x$overall_rmse, 4)),
  MAE = sapply(all_evals, function(x) round(x$overall_mae, 4))
) %>%
  filter(!is.na(Correlation)) %>%
  arrange(desc(Correlation))

message("\nFinal Method Ranking:")
print(summary_table)

# Per cell type
celltype_cors_df <- do.call(rbind, lapply(names(all_evals), function(m) {
  if (!is.na(all_evals[[m]]$overall_cor)) {
    data.frame(
      Method = m,
      CellType = names(all_evals[[m]]$celltype_cors),
      Correlation = all_evals[[m]]$celltype_cors
    )
  }
}))

message("\nPer Cell Type Correlations:")
celltype_cors_wide <- celltype_cors_df %>%
  pivot_wider(names_from = Method, values_from = Correlation)
print(celltype_cors_wide)

# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------
message("\n=== Generating updated figures ===")

# Color palette
method_colors <- c("MuSiC" = "#66C2A5", "BayesPrism" = "#FC8D62",
                   "NNLS" = "#8DA0CB", "BisqueRNA" = "#E78AC3",
                   "dtangle" = "#A6D854", "CIBERSORTx" = "#FFD92F")

# Figure 1: Method accuracy bar plot
p1 <- ggplot(summary_table, aes(x = reorder(Method, Correlation), y = Correlation,
                                  fill = Method)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", Correlation)), hjust = -0.1) +
  scale_fill_manual(values = method_colors) +
  coord_flip() + ylim(0, 1) + theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Deconvolution Accuracy (Pseudo-bulk Validation)",
       subtitle = "Including CIBERSORTx",
       x = "", y = "Correlation with Ground Truth")

ggsave(file.path(VALIDATION_DIR, "method_accuracy_with_cibersortx.png"), p1,
       width = 8, height = 6, dpi = 300)

# Figure 2: Scatter plot for CIBERSORTx
common_samples <- intersect(rownames(cibersortx_props), rownames(ground_truth))
common_types <- intersect(colnames(cibersortx_props), colnames(ground_truth))

scatter_cibersortx <- data.frame(
  Sample = rep(common_samples, length(common_types)),
  CellType = rep(common_types, each = length(common_samples)),
  Estimated = as.vector(cibersortx_props[common_samples, common_types]),
  GroundTruth = as.vector(ground_truth[common_samples, common_types])
)

p2 <- ggplot(scatter_cibersortx, aes(x = GroundTruth, y = Estimated, color = CellType)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = paste0("CIBERSORTx: Estimated vs Ground Truth (r = ",
                      round(cibersortx_eval$overall_cor, 3), ")"),
       x = "Ground Truth Proportion", y = "Estimated Proportion")

ggsave(file.path(VALIDATION_DIR, "cibersortx_scatter.png"), p2,
       width = 8, height = 6, dpi = 300)

# Figure 3: Per cell type heatmap
celltype_perf_matrix <- celltype_cors_df %>%
  pivot_wider(names_from = CellType, values_from = Correlation) %>%
  tibble::column_to_rownames("Method") %>%
  as.matrix()

png(file.path(VALIDATION_DIR, "celltype_performance_all_methods.png"),
    width = 12, height = 6, units = "in", res = 300)
pheatmap(celltype_perf_matrix,
         display_numbers = TRUE,
         number_format = "%.2f",
         color = colorRampPalette(c("white", "steelblue"))(100),
         main = "Per Cell Type Correlation - All Methods",
         cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

# Figure 4: All methods scatter
scatter_all <- do.call(rbind, lapply(names(results_list), function(m) {
  est <- results_list[[m]]
  cs <- intersect(rownames(est), rownames(ground_truth))
  ct <- intersect(colnames(est), colnames(ground_truth))
  if (length(cs) > 0 && length(ct) > 0) {
    data.frame(
      Method = m,
      Sample = rep(cs, length(ct)),
      CellType = rep(ct, each = length(cs)),
      Estimated = as.vector(est[cs, ct]),
      GroundTruth = as.vector(ground_truth[cs, ct])
    )
  }
}))

p4 <- ggplot(scatter_all, aes(x = GroundTruth, y = Estimated, color = CellType)) +
  geom_point(alpha = 0.4, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~Method, ncol = 3) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "All Methods: Estimated vs Ground Truth",
       x = "Ground Truth", y = "Estimated")

ggsave(file.path(VALIDATION_DIR, "all_methods_scatter.png"), p4,
       width = 12, height = 10, dpi = 300)

# Combined summary
combined <- (p1 + p2) / p4 +
  plot_annotation(title = "Pseudo-bulk Validation: All Methods Including CIBERSORTx")

ggsave(file.path(VALIDATION_DIR, "validation_summary_all.png"), combined,
       width = 14, height = 16, dpi = 300)

# -----------------------------------------------------------------------------
# Save Updated Results
# -----------------------------------------------------------------------------
message("\n=== Saving updated results ===")

saveRDS(results_list, file.path(VALIDATION_DIR, "estimated_proportions_all.rds"))
saveRDS(all_evals, file.path(VALIDATION_DIR, "evaluation_results_all.rds"))
write.csv(summary_table, file.path(VALIDATION_DIR, "method_performance_all.csv"), row.names = FALSE)
write.csv(celltype_cors_wide, file.path(VALIDATION_DIR, "celltype_correlations_all.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Final Summary
# -----------------------------------------------------------------------------
message("\n", strrep("=", 70))
message("VALIDATION COMPLETE - ALL METHODS")
message(strrep("=", 70))

message("\nFinal Method Ranking:")
for (i in 1:nrow(summary_table)) {
  message(sprintf("  %d. %s: r = %.4f, RMSE = %.4f",
                  i, summary_table$Method[i],
                  summary_table$Correlation[i], summary_table$RMSE[i]))
}

message("\nFiles saved to output/validation/")
message("\n=== Done ===")
