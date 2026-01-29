# =============================================================================
# 10_compare_with_cibersortx.R - Compare All 6 Deconvolution Methods
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
CIBERSORTX_DIR <- file.path(OUTPUT_DIR, "cibersortx")

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

# Load CIBERSORTx results
cibersortx_raw <- read.delim(file.path(CIBERSORTX_DIR, "CIBERSORTx_Job1_Results.txt"),
                              stringsAsFactors = FALSE)

message(paste("MuSiC:", nrow(music_props), "samples x", ncol(music_props), "cell types"))
message(paste("NNLS:", nrow(nnls_props), "samples x", ncol(nnls_props), "cell types"))
message(paste("BisqueRNA:", nrow(bisque_props), "samples x", ncol(bisque_props), "cell types"))
message(paste("dtangle:", nrow(dtangle_props), "samples x", ncol(dtangle_props), "cell types"))
message(paste("BayesPrism:", nrow(bp_props), "samples x", ncol(bp_props), "cell types"))
message(paste("CIBERSORTx:", nrow(cibersortx_raw), "samples x",
              ncol(cibersortx_raw) - 4, "cell types (excluding QC columns)"))

# -----------------------------------------------------------------------------
# Process CIBERSORTx Results
# -----------------------------------------------------------------------------
message("\n=== Processing CIBERSORTx results ===")

# Extract sample names and cell type proportions
cibersortx_samples <- cibersortx_raw$Mixture
cibersortx_celltypes <- colnames(cibersortx_raw)[2:(ncol(cibersortx_raw) - 3)]
cibersortx_props <- as.matrix(cibersortx_raw[, cibersortx_celltypes])
rownames(cibersortx_props) <- cibersortx_samples

# Fix cell type names (CIBERSORTx uses dots instead of dashes/spaces)
colnames(cibersortx_props) <- gsub("\\.", " ", colnames(cibersortx_props))
colnames(cibersortx_props) <- gsub("B cells", "B-cells", colnames(cibersortx_props))
colnames(cibersortx_props) <- gsub("T cells", "T-cells", colnames(cibersortx_props))

message(paste("CIBERSORTx cell types:", paste(colnames(cibersortx_props), collapse = ", ")))

# CIBERSORTx QC metrics
message("\nCIBERSORTx QC metrics:")
message(sprintf("  Mean P-value: %.4f", mean(cibersortx_raw$P.value)))
message(sprintf("  Mean Correlation: %.4f", mean(cibersortx_raw$Correlation)))
message(sprintf("  Mean RMSE: %.4f", mean(cibersortx_raw$RMSE)))

# -----------------------------------------------------------------------------
# Standardize BayesPrism output format
# -----------------------------------------------------------------------------
if ("Sample" %in% colnames(bp_props)) {
  bp_samples <- bp_props$Sample
  bp_props_mat <- as.matrix(bp_props[, -which(colnames(bp_props) == "Sample")])
  rownames(bp_props_mat) <- bp_samples
} else {
  bp_props_mat <- as.matrix(bp_props)
}

# Get common samples and cell types from reference method (MuSiC)
all_samples <- rownames(music_props)
all_celltypes <- colnames(music_props)

message(paste("\nReference samples:", length(all_samples)))
message(paste("Reference cell types:", paste(all_celltypes, collapse = ", ")))

# -----------------------------------------------------------------------------
# Reorder function to match dimensions
# -----------------------------------------------------------------------------
reorder_props <- function(props, target_cols, target_rows) {
  # Match rows (samples)
  common_rows <- intersect(rownames(props), target_rows)

  if (length(common_rows) == 0) {
    message("Warning: No common samples found. Using available samples.")
    common_rows <- rownames(props)
  }

  props <- props[common_rows, , drop = FALSE]

  # Match columns (cell types)
  matched_cols <- intersect(colnames(props), target_cols)

  result <- matrix(0, nrow = length(common_rows), ncol = length(target_cols))
  rownames(result) <- common_rows
  colnames(result) <- target_cols

  # Fill in available data
  for (col in matched_cols) {
    result[, col] <- props[, col]
  }

  # Renormalize rows to sum to 1
  row_sums <- rowSums(result)
  result <- result / ifelse(row_sums > 0, row_sums, 1)
  result[is.nan(result)] <- 0

  return(result)
}

# Reorder all methods to match reference
bp_props_ordered <- reorder_props(bp_props_mat, all_celltypes, all_samples)
cibersortx_props_ordered <- reorder_props(cibersortx_props, all_celltypes, all_samples)

message(paste("\nCIBERSORTx samples matched:", nrow(cibersortx_props_ordered)))

# -----------------------------------------------------------------------------
# Create Results List
# -----------------------------------------------------------------------------
results_list <- list(
  MuSiC = music_props,
  NNLS = nnls_props,
  BisqueRNA = bisque_props,
  dtangle = dtangle_props,
  BayesPrism = bp_props_ordered,
  CIBERSORTx = cibersortx_props_ordered
)

# Save CIBERSORTx processed results
saveRDS(cibersortx_props_ordered, file.path(DATA_DIR, "props_cibersortx.rds"))

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
  BayesPrism = colMeans(bp_props_ordered) * 100,
  CIBERSORTx = colMeans(cibersortx_props_ordered) * 100
)

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

# Save results
write.csv(cor_matrix, file.path(DATA_DIR, "method_correlation_matrix_6methods.csv"))
write.csv(mean_props, file.path(DATA_DIR, "method_comparison_6methods.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------
message("\n=== Generating comparison figures ===")

# Color palette for 6 methods
method_colors <- c("MuSiC" = "#66C2A5", "BayesPrism" = "#FC8D62",
                   "NNLS" = "#8DA0CB", "BisqueRNA" = "#E78AC3",
                   "dtangle" = "#A6D854", "CIBERSORTx" = "#FFD92F")

# Figure 1: Mean proportions by method (all 6 methods)
comparison_long <- mean_props %>%
  pivot_longer(-CellType, names_to = "Method", values_to = "Proportion") %>%
  mutate(
    CellType = factor(CellType, levels = mean_props$CellType[order(-mean_props$MuSiC)]),
    Method = factor(Method, levels = c("MuSiC", "BayesPrism", "CIBERSORTx",
                                       "NNLS", "BisqueRNA", "dtangle"))
  )

p1 <- ggplot(comparison_long, aes(x = CellType, y = Proportion, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = method_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cell Type Proportions by Method (6 Methods)",
       x = "", y = "Mean Proportion (%)")

ggsave(file.path(FIGURE_DIR, "method_comparison_6methods_barplot.png"), p1,
       width = 16, height = 6, dpi = 300)
message("Saved: method_comparison_6methods_barplot.png")

# Figure 2: Correlation heatmap
png(file.path(FIGURE_DIR, "method_correlation_6methods_heatmap.png"),
    width = 8, height = 7, units = "in", res = 300)
pheatmap(cor_matrix,
         display_numbers = TRUE,
         number_format = "%.2f",
         color = colorRampPalette(c("white", "steelblue"))(100),
         main = "Method Agreement (6 Methods)")
dev.off()
message("Saved: method_correlation_6methods_heatmap.png")

# Figure 3: CIBERSORTx vs other methods scatter plots
scatter_cibersortx <- data.frame(
  Sample = rep(rownames(cibersortx_props_ordered), 5 * length(all_celltypes)),
  CellType = rep(rep(all_celltypes, each = nrow(cibersortx_props_ordered)), 5),
  CIBERSORTx = rep(as.vector(t(cibersortx_props_ordered)), 5),
  Other = c(as.vector(t(music_props)),
            as.vector(t(bp_props_ordered)),
            as.vector(t(nnls_props)),
            as.vector(t(bisque_props)),
            as.vector(t(dtangle_props))),
  Method = rep(c("MuSiC", "BayesPrism", "NNLS", "BisqueRNA", "dtangle"),
               each = nrow(cibersortx_props_ordered) * length(all_celltypes))
)

p3 <- ggplot(scatter_cibersortx, aes(x = CIBERSORTx, y = Other, color = CellType)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~Method, ncol = 3) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "CIBERSORTx vs Other Methods",
       x = "CIBERSORTx Proportion", y = "Other Method Proportion")

ggsave(file.path(FIGURE_DIR, "cibersortx_vs_all.png"), p3,
       width = 12, height = 10, dpi = 300)
message("Saved: cibersortx_vs_all.png")

# Figure 4: Consensus proportions (6 methods)
consensus <- rowMeans(cbind(
  colMeans(music_props),
  colMeans(nnls_props),
  colMeans(bisque_props),
  colMeans(dtangle_props),
  colMeans(bp_props_ordered),
  colMeans(cibersortx_props_ordered)
)) * 100

consensus_sd <- apply(cbind(
  colMeans(music_props),
  colMeans(nnls_props),
  colMeans(bisque_props),
  colMeans(dtangle_props),
  colMeans(bp_props_ordered),
  colMeans(cibersortx_props_ordered)
), 1, sd) * 100

consensus_df <- data.frame(
  CellType = names(consensus),
  Consensus = consensus,
  SD = consensus_sd
) %>% arrange(desc(Consensus))

p4 <- ggplot(consensus_df, aes(x = reorder(CellType, -Consensus), y = Consensus)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_errorbar(aes(ymin = pmax(0, Consensus - SD), ymax = Consensus + SD),
                width = 0.3, color = "gray40") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Consensus Cell Type Proportions (6 Methods)",
       subtitle = "Error bars show SD across methods",
       x = "", y = "Mean Proportion (%)")

ggsave(file.path(FIGURE_DIR, "consensus_proportions_6methods.png"), p4,
       width = 10, height = 6, dpi = 300)
message("Saved: consensus_proportions_6methods.png")

# Figure 5: Method variability
sample_variability <- data.frame(Sample = rownames(cibersortx_props_ordered))

for (ct in all_celltypes) {
  ct_vals <- cbind(music_props[rownames(cibersortx_props_ordered), ct],
                   nnls_props[rownames(cibersortx_props_ordered), ct],
                   bisque_props[rownames(cibersortx_props_ordered), ct],
                   dtangle_props[rownames(cibersortx_props_ordered), ct],
                   bp_props_ordered[, ct],
                   cibersortx_props_ordered[, ct])
  sample_variability[[ct]] <- apply(ct_vals, 1, sd, na.rm = TRUE) * 100
}

sample_var_long <- sample_variability %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "SD")

p5 <- ggplot(sample_var_long, aes(x = CellType, y = SD)) +
  geom_boxplot(fill = "lightblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Method Variability per Cell Type (6 Methods)",
       x = "", y = "Standard Deviation (%)")

ggsave(file.path(FIGURE_DIR, "method_variability_6methods.png"), p5,
       width = 10, height = 6, dpi = 300)
message("Saved: method_variability_6methods.png")

# Figure 6: Per-sample composition comparison
# Heatmap showing CIBERSORTx vs MuSiC vs BayesPrism for each sample
top_methods <- c("MuSiC", "BayesPrism", "CIBERSORTx")
sample_comparison <- do.call(rbind, lapply(top_methods, function(m) {
  props <- results_list[[m]]
  data.frame(
    Sample = rownames(props),
    Method = m,
    as.data.frame(props)
  )
}))

sample_comp_long <- sample_comparison %>%
  pivot_longer(cols = -c(Sample, Method), names_to = "CellType", values_to = "Proportion")

p6 <- ggplot(sample_comp_long, aes(x = Sample, y = CellType, fill = Proportion)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", limits = c(0, 1)) +
  facet_wrap(~Method, ncol = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  labs(title = "Sample Composition: MuSiC vs BayesPrism vs CIBERSORTx",
       x = "", y = "")

ggsave(file.path(FIGURE_DIR, "sample_composition_3methods.png"), p6,
       width = 14, height = 10, dpi = 300)
message("Saved: sample_composition_3methods.png")

# Combined summary figure
combined_fig <- (p1 + p4) / (p3 + p5) +
  plot_annotation(
    title = "Bulk Deconvolution: 6-Method Comparison",
    subtitle = "MuSiC, NNLS, BisqueRNA, dtangle, BayesPrism, CIBERSORTx",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(file.path(FIGURE_DIR, "method_comparison_6methods_summary.png"), combined_fig,
       width = 18, height = 14, dpi = 300)
message("Saved: method_comparison_6methods_summary.png")

# -----------------------------------------------------------------------------
# Summary Statistics
# -----------------------------------------------------------------------------
message("\n", strrep("=", 70))
message("6-METHOD COMPARISON SUMMARY")
message(strrep("=", 70))

message("\nMean Proportions by Method:")
print(mean_props_print)

message("\nMethod Agreement (Correlation with CIBERSORTx):")
cibersortx_cors <- cor_matrix["CIBERSORTx", ]
for (m in names(sort(cibersortx_cors, decreasing = TRUE))) {
  message(sprintf("  %s: r = %.3f", m, cibersortx_cors[m]))
}

message("\nConsensus Proportions (Mean ± SD across 6 methods):")
consensus_print <- consensus_df
consensus_print$Summary <- sprintf("%.1f ± %.1f", consensus_print$Consensus, consensus_print$SD)
print(consensus_print[, c("CellType", "Summary")])

# Save consensus
write.csv(consensus_df, file.path(DATA_DIR, "consensus_proportions_6methods.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Method Ranking
# -----------------------------------------------------------------------------
message("\n", strrep("=", 70))
message("METHOD EVALUATION")
message(strrep("=", 70))

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
  CAFs = sapply(results_list, function(x) mean(x[, "CAFs"]) * 100),
  Endothelial = sapply(results_list, function(x) mean(x[, "Endothelial"]) * 100)
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

message("\n=== CIBERSORTx Analysis ===")
message(sprintf("CIBERSORTx estimates %.1f%% Cancer Epithelial",
                mean(cibersortx_props_ordered[, "Cancer Epithelial"]) * 100))
message(sprintf("CIBERSORTx estimates %.1f%% CAFs",
                mean(cibersortx_props_ordered[, "CAFs"]) * 100))
message(sprintf("CIBERSORTx estimates %.1f%% Endothelial",
                mean(cibersortx_props_ordered[, "Endothelial"]) * 100))
message(sprintf("Highest correlation: %s (r = %.3f)",
                names(which.max(cibersortx_cors[-6])), max(cibersortx_cors[-6])))

message("\n=== Comparison complete ===")
