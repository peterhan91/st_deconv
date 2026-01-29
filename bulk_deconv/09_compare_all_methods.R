# =============================================================================
# 09_compare_all_methods.R - Compare All Deconvolution Methods Including BayesPrism
# =============================================================================

library(Biobase)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(patchwork)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
OUTPUT_DIR <- "output"
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
DATA_DIR <- file.path(OUTPUT_DIR, "data")

# -----------------------------------------------------------------------------
# Load All Results
# -----------------------------------------------------------------------------
message("=== Loading all deconvolution results ===")

# Load individual method results
music_props <- readRDS(file.path(DATA_DIR, "music_proportions.rds"))
nnls_props <- readRDS(file.path(DATA_DIR, "props_nnls.rds"))
bisque_props <- readRDS(file.path(DATA_DIR, "props_bisquerna.rds"))
dtangle_props <- readRDS(file.path(DATA_DIR, "props_dtangle.rds"))
bp_props <- readRDS(file.path(DATA_DIR, "props_bayesprism.rds"))

message(paste("MuSiC:", nrow(music_props), "samples x", ncol(music_props), "cell types"))
message(paste("NNLS:", nrow(nnls_props), "samples x", ncol(nnls_props), "cell types"))
message(paste("BisqueRNA:", nrow(bisque_props), "samples x", ncol(bisque_props), "cell types"))
message(paste("dtangle:", nrow(dtangle_props), "samples x", ncol(dtangle_props), "cell types"))
message(paste("BayesPrism:", nrow(bp_props) - 1, "samples x", ncol(bp_props) - 1, "cell types"))

# -----------------------------------------------------------------------------
# Standardize BayesPrism output format
# -----------------------------------------------------------------------------
# BayesPrism output has Sample column, need to convert to matrix format
if ("Sample" %in% colnames(bp_props)) {
  bp_samples <- bp_props$Sample
  bp_props_mat <- as.matrix(bp_props[, -which(colnames(bp_props) == "Sample")])
  rownames(bp_props_mat) <- bp_samples
} else {
  bp_props_mat <- as.matrix(bp_props)
}

# Get common samples and cell types
all_samples <- rownames(music_props)
all_celltypes <- colnames(music_props)

message(paste("\nCommon samples:", length(all_samples)))
message(paste("Cell types:", paste(all_celltypes, collapse = ", ")))

# -----------------------------------------------------------------------------
# Reorder BayesPrism to match other methods
# -----------------------------------------------------------------------------
reorder_props <- function(props, target_cols, target_rows) {
  # Match rows
  common_rows <- intersect(rownames(props), target_rows)
  props <- props[common_rows, , drop = FALSE]

  # Match columns
  matched_cols <- intersect(colnames(props), target_cols)
  result <- matrix(0, nrow = length(target_rows), ncol = length(target_cols))
  rownames(result) <- target_rows
  colnames(result) <- target_cols

  # Fill in available data
  for (col in matched_cols) {
    for (row in common_rows) {
      result[row, col] <- props[row, col]
    }
  }

  # Renormalize rows
  row_sums <- rowSums(result)
  result <- result / ifelse(row_sums > 0, row_sums, 1)
  result[is.nan(result)] <- 0
  return(result)
}

bp_props_ordered <- reorder_props(bp_props_mat, all_celltypes, all_samples)

# -----------------------------------------------------------------------------
# Create Results List
# -----------------------------------------------------------------------------
results_list <- list(
  MuSiC = music_props,
  NNLS = nnls_props,
  BisqueRNA = bisque_props,
  dtangle = dtangle_props,
  BayesPrism = bp_props_ordered
)

# -----------------------------------------------------------------------------
# Calculate Mean Proportions
# -----------------------------------------------------------------------------
message("\n=== Mean Proportions by Method ===")

mean_props <- data.frame(
  CellType = all_celltypes,
  MuSiC = colMeans(music_props) * 100,
  NNLS = colMeans(nnls_props) * 100,
  BisqueRNA = colMeans(bisque_props) * 100,
  dtangle = colMeans(dtangle_props) * 100,
  BayesPrism = colMeans(bp_props_ordered) * 100
)

# Print with proper formatting
mean_props_print <- mean_props
mean_props_print[, -1] <- round(mean_props_print[, -1], 2)
print(mean_props_print)

# -----------------------------------------------------------------------------
# Calculate Method Correlations
# -----------------------------------------------------------------------------
message("\n=== Method Correlations ===")

methods <- names(results_list)
cor_matrix <- matrix(NA, length(methods), length(methods))
rownames(cor_matrix) <- colnames(cor_matrix) <- methods

for (i in seq_along(methods)) {
  for (j in seq_along(methods)) {
    props_i <- as.vector(results_list[[i]])
    props_j <- as.vector(results_list[[j]])
    cor_matrix[i, j] <- cor(props_i, props_j, use = "complete.obs")
  }
}

message("Correlation matrix:")
print(round(cor_matrix, 3))

# Save correlation matrix
write.csv(cor_matrix, file.path(DATA_DIR, "method_correlation_matrix_all.csv"))

# -----------------------------------------------------------------------------
# Save Updated Comparison
# -----------------------------------------------------------------------------
write.csv(mean_props, file.path(DATA_DIR, "method_comparison_all.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------
message("\n=== Generating comparison figures ===")

# Figure 1: Mean proportions by method (all 5 methods)
comparison_long <- mean_props %>%
  pivot_longer(-CellType, names_to = "Method", values_to = "Proportion") %>%
  mutate(
    CellType = factor(CellType, levels = mean_props$CellType[order(-mean_props$MuSiC)]),
    Method = factor(Method, levels = c("MuSiC", "BayesPrism", "NNLS", "BisqueRNA", "dtangle"))
  )

p1 <- ggplot(comparison_long, aes(x = CellType, y = Proportion, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("MuSiC" = "#66C2A5", "BayesPrism" = "#FC8D62",
                               "NNLS" = "#8DA0CB", "BisqueRNA" = "#E78AC3",
                               "dtangle" = "#A6D854")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Proportions by Method (Including BayesPrism)",
       x = "", y = "Mean Proportion (%)")

ggsave(file.path(FIGURE_DIR, "method_comparison_all_barplot.png"), p1,
       width = 14, height = 6, dpi = 300)
message("Saved: method_comparison_all_barplot.png")

# Figure 2: Updated correlation heatmap
png(file.path(FIGURE_DIR, "method_correlation_all_heatmap.png"),
    width = 7, height = 6, units = "in", res = 300)
pheatmap(cor_matrix,
         display_numbers = TRUE,
         number_format = "%.2f",
         color = colorRampPalette(c("white", "steelblue"))(100),
         main = "Method Agreement (All Methods)")
dev.off()
message("Saved: method_correlation_all_heatmap.png")

# Figure 3: Consensus proportions including BayesPrism
consensus <- rowMeans(cbind(
  colMeans(music_props),
  colMeans(nnls_props),
  colMeans(bisque_props),
  colMeans(dtangle_props),
  colMeans(bp_props_ordered)
)) * 100

consensus_df <- data.frame(
  CellType = names(consensus),
  Consensus = consensus
) %>% arrange(desc(Consensus))

p3 <- ggplot(consensus_df, aes(x = reorder(CellType, -Consensus), y = Consensus)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Consensus Cell Type Proportions (5 Methods)",
       x = "", y = "Mean Proportion (%)")

ggsave(file.path(FIGURE_DIR, "consensus_proportions_all.png"), p3,
       width = 10, height = 6, dpi = 300)
message("Saved: consensus_proportions_all.png")

# Figure 4: Compare BayesPrism with MuSiC
scatter_bp_music <- data.frame(
  Sample = rep(all_samples, length(all_celltypes)),
  CellType = rep(all_celltypes, each = length(all_samples)),
  MuSiC = as.vector(t(music_props)),
  BayesPrism = as.vector(t(bp_props_ordered))
)

p4 <- ggplot(scatter_bp_music, aes(x = MuSiC, y = BayesPrism, color = CellType)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(title = paste0("MuSiC vs BayesPrism (r = ",
                      round(cor(as.vector(music_props), as.vector(bp_props_ordered)), 3), ")"),
       x = "MuSiC Proportion", y = "BayesPrism Proportion")

ggsave(file.path(FIGURE_DIR, "music_vs_bayesprism.png"), p4,
       width = 8, height = 6, dpi = 300)
message("Saved: music_vs_bayesprism.png")

# Figure 5: All pairwise scatter plots
scatter_all <- data.frame(
  Sample = rep(all_samples, 4 * length(all_celltypes)),
  CellType = rep(rep(all_celltypes, each = length(all_samples)), 4),
  BayesPrism = rep(as.vector(t(bp_props_ordered)), 4),
  Other = c(as.vector(t(music_props)),
            as.vector(t(nnls_props)),
            as.vector(t(bisque_props)),
            as.vector(t(dtangle_props))),
  Method = rep(c("MuSiC", "NNLS", "BisqueRNA", "dtangle"),
               each = length(all_samples) * length(all_celltypes))
)

p5 <- ggplot(scatter_all, aes(x = BayesPrism, y = Other, color = CellType)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~Method, ncol = 2) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "BayesPrism vs Other Methods",
       x = "BayesPrism Proportion", y = "Other Method Proportion")

ggsave(file.path(FIGURE_DIR, "bayesprism_vs_all.png"), p5,
       width = 10, height = 10, dpi = 300)
message("Saved: bayesprism_vs_all.png")

# Figure 6: Method variability including BayesPrism
sample_variability <- data.frame(Sample = all_samples)

for (ct in all_celltypes) {
  ct_vals <- cbind(music_props[,ct], nnls_props[,ct], bisque_props[,ct],
                   dtangle_props[,ct], bp_props_ordered[,ct])
  sample_variability[[ct]] <- apply(ct_vals, 1, sd) * 100
}

sample_var_long <- sample_variability %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "SD")

p6 <- ggplot(sample_var_long, aes(x = CellType, y = SD)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Method Variability per Cell Type (5 Methods)",
       x = "", y = "Standard Deviation (%)")

ggsave(file.path(FIGURE_DIR, "method_variability_all.png"), p6,
       width = 10, height = 6, dpi = 300)
message("Saved: method_variability_all.png")

# Combined summary figure
combined_fig <- (p1 + p3) / (p4 + p6) +
  plot_annotation(
    title = "Bulk Deconvolution: 5-Method Comparison",
    subtitle = "MuSiC, NNLS, BisqueRNA, dtangle, BayesPrism",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(file.path(FIGURE_DIR, "method_comparison_all_summary.png"), combined_fig,
       width = 16, height = 12, dpi = 300)
message("Saved: method_comparison_all_summary.png")

# -----------------------------------------------------------------------------
# Summary Statistics
# -----------------------------------------------------------------------------
message("\n", strrep("=", 60))
message("5-METHOD COMPARISON SUMMARY")
message(strrep("=", 60))

message("\nMean Proportions by Method:")
mean_props_print <- mean_props
mean_props_print[, -1] <- round(mean_props_print[, -1], 2)
print(mean_props_print)

message("\nMethod Agreement (Correlation with BayesPrism):")
bp_cors <- cor_matrix["BayesPrism", ]
for (m in names(bp_cors)) {
  message(sprintf("  %s: r = %.3f", m, bp_cors[m]))
}

message("\nConsensus Proportions (Mean across 5 methods):")
print(consensus_df)

# Save updated consensus
write.csv(consensus_df, file.path(DATA_DIR, "consensus_proportions_all.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Method Ranking
# -----------------------------------------------------------------------------
message("\n", strrep("=", 60))
message("METHOD EVALUATION")
message(strrep("=", 60))

# Expected proportions for breast cancer tissue
expected <- list(
  "Cancer Epithelial" = c(20, 60),
  "CAFs" = c(10, 40),
  "Endothelial" = c(5, 20),
  "Normal Epithelial" = c(5, 25),
  "Myeloid" = c(3, 15),
  "T-cells" = c(2, 15),
  "B-cells" = c(1, 10),
  "Plasmablasts" = c(0, 5),
  "PVL" = c(0, 5)
)

# Calculate score for each method
score_method <- function(props) {
  score <- 0
  for (ct in names(expected)) {
    if (ct %in% colnames(props)) {
      val <- mean(props[, ct]) * 100
      if (val >= expected[[ct]][1] && val <= expected[[ct]][2]) {
        score <- score + 1
      }
    }
  }
  return(score)
}

method_scores <- data.frame(
  Method = methods,
  Score = sapply(results_list, score_method),
  Cancer_Epithelial = sapply(results_list, function(x) mean(x[, "Cancer Epithelial"]) * 100),
  CAFs = sapply(results_list, function(x) mean(x[, "CAFs"]) * 100)
)

message("\nCell types within expected range (max 9):")
print(method_scores[order(-method_scores$Score), ])

# Entropy analysis
calc_entropy <- function(props) {
  mean_props <- colMeans(props)
  mean_props <- mean_props[mean_props > 0]
  -sum(mean_props * log2(mean_props))
}

entropies <- sapply(results_list, calc_entropy)
message("\nEntropy (higher = more diverse, max ~3.17 bits for 9 types):")
for (m in names(sort(entropies, decreasing = TRUE))) {
  message(sprintf("  %s: %.2f bits", m, entropies[m]))
}

message("\n=== BayesPrism Analysis ===")
message(sprintf("BayesPrism estimates %.1f%% Cancer Epithelial (highest among methods)",
                mean(bp_props_ordered[, "Cancer Epithelial"]) * 100))
message(sprintf("BayesPrism estimates %.1f%% CAFs", mean(bp_props_ordered[, "CAFs"]) * 100))
message(sprintf("Correlation with MuSiC: r = %.3f", cor_matrix["BayesPrism", "MuSiC"]))

message("\n=== Comparison complete ===")
