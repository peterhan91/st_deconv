# =============================================================================
# 06_quality_assessment.R - Evaluate Deconvolution Quality
# ST Deconvolution Pipeline for IDC Data using RCTD
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
OUTPUT_DIR <- "output/data"
FIGURE_DIR <- "output/figures"

# -----------------------------------------------------------------------------
# Load Results
# -----------------------------------------------------------------------------
message("=== Loading deconvolution results ===")

weights <- readRDS(file.path(OUTPUT_DIR, "cell_type_weights.rds"))
weights <- as.matrix(weights)

coords <- read.csv(file.path(OUTPUT_DIR, "spatial_coords.csv"), row.names = 1)

# Original spot count
total_spots <- nrow(coords)
converged_spots <- nrow(weights)

message(paste("Total spots:", total_spots))
message(paste("Converged spots:", converged_spots))

# =============================================================================
# QUALITY METRIC 1: Convergence Rate
# =============================================================================
message("\n" , strrep("=", 60))
message("METRIC 1: CONVERGENCE RATE")
message(strrep("=", 60))

convergence_rate <- converged_spots / total_spots * 100
message(paste("Convergence rate:", round(convergence_rate, 1), "%"))

if (convergence_rate >= 90) {
  message("Assessment: EXCELLENT (>=90%)")
  convergence_score <- 3
} else if (convergence_rate >= 80) {
  message("Assessment: GOOD (80-90%)")
  convergence_score <- 2
} else if (convergence_rate >= 60) {
  message("Assessment: ACCEPTABLE (60-80%)")
  convergence_score <- 1
} else {
  message("Assessment: POOR (<60%) - Many spots failed to converge")
  message("  Possible causes: low UMI counts, poor gene overlap, noisy data")
  convergence_score <- 0
}

# =============================================================================
# QUALITY METRIC 2: Weight Normalization
# =============================================================================
message("\n", strrep("=", 60))
message("METRIC 2: WEIGHT NORMALIZATION")
message(strrep("=", 60))

weight_sums <- rowSums(weights)
message(paste("Mean weight sum:", round(mean(weight_sums), 4)))
message(paste("Std dev:", round(sd(weight_sums), 4)))
message(paste("Range:", round(min(weight_sums), 4), "-", round(max(weight_sums), 4)))

if (abs(mean(weight_sums) - 1) < 0.001 && sd(weight_sums) < 0.01) {
  message("Assessment: EXCELLENT - Weights properly normalized")
  normalization_score <- 3
} else if (abs(mean(weight_sums) - 1) < 0.01) {
  message("Assessment: GOOD - Minor variation acceptable")
  normalization_score <- 2
} else {
  message("Assessment: CHECK - Weights may not be properly normalized")
  normalization_score <- 1
}

# =============================================================================
# QUALITY METRIC 3: Cell Type Separation
# =============================================================================
message("\n", strrep("=", 60))
message("METRIC 3: CELL TYPE SEPARATION")
message(strrep("=", 60))

# Calculate correlation between cell types (negative = good separation)
cor_matrix <- cor(weights)
off_diag <- cor_matrix[upper.tri(cor_matrix)]

message(paste("Mean pairwise correlation:", round(mean(off_diag), 3)))
message(paste("Max pairwise correlation:", round(max(off_diag), 3)))

# Check for highly correlated pairs (problematic if >0.5)
high_cor_pairs <- which(cor_matrix > 0.5 & upper.tri(cor_matrix), arr.ind = TRUE)
if (nrow(high_cor_pairs) > 0) {
  message("\nHighly correlated cell type pairs (r > 0.5):")
  for (i in 1:nrow(high_cor_pairs)) {
    ct1 <- colnames(cor_matrix)[high_cor_pairs[i, 1]]
    ct2 <- colnames(cor_matrix)[high_cor_pairs[i, 2]]
    r <- cor_matrix[high_cor_pairs[i, 1], high_cor_pairs[i, 2]]
    message(paste("  ", ct1, "-", ct2, ": r =", round(r, 3)))
  }
  message("  Note: High correlation may indicate difficulty distinguishing these types")
}

if (mean(off_diag) < 0) {
  message("Assessment: EXCELLENT - Cell types show negative correlation (expected)")
  separation_score <- 3
} else if (max(off_diag) < 0.5) {
  message("Assessment: GOOD - No highly correlated pairs")
  separation_score <- 2
} else {
  message("Assessment: MODERATE - Some cell types may be hard to distinguish")
  separation_score <- 1
}

# =============================================================================
# QUALITY METRIC 4: Classification Confidence
# =============================================================================
message("\n", strrep("=", 60))
message("METRIC 4: CLASSIFICATION CONFIDENCE")
message(strrep("=", 60))

max_weights <- apply(weights, 1, max)
second_max <- apply(weights, 1, function(x) sort(x, decreasing = TRUE)[2])
confidence_margin <- max_weights - second_max

message(paste("Mean dominant cell type proportion:", round(mean(max_weights), 3)))
message(paste("Mean confidence margin (1st - 2nd):", round(mean(confidence_margin), 3)))
message(paste("Spots with clear dominant (>50%):", sum(max_weights > 0.5),
              "(", round(sum(max_weights > 0.5)/converged_spots*100, 1), "%)"))
message(paste("Spots with dominant >30%:", sum(max_weights > 0.3),
              "(", round(sum(max_weights > 0.3)/converged_spots*100, 1), "%)"))

# Mixed spots (no type >30%) may indicate transition zones or technical issues
mixed_spots <- sum(max_weights < 0.3)
message(paste("Mixed spots (<30% dominant):", mixed_spots,
              "(", round(mixed_spots/converged_spots*100, 1), "%)"))

if (mean(max_weights) > 0.5) {
  message("Assessment: HIGH CONFIDENCE - Most spots have clear dominant type")
  confidence_score <- 3
} else if (mean(max_weights) > 0.35) {
  message("Assessment: MODERATE - Expected for heterogeneous tissue like tumors")
  confidence_score <- 2
} else {
  message("Assessment: LOW - Many mixed spots, may need more marker genes")
  confidence_score <- 1
}

# =============================================================================
# QUALITY METRIC 5: Biological Plausibility
# =============================================================================
message("\n", strrep("=", 60))
message("METRIC 5: BIOLOGICAL PLAUSIBILITY (IDC-specific)")
message(strrep("=", 60))

mean_props <- colMeans(weights) * 100

message("\nCell type proportions vs. expected for IDC:")
message(sprintf("  %-18s: %5.1f%% - %s", "CAFs", mean_props["CAFs"],
                ifelse(mean_props["CAFs"] > 15, "EXPECTED (high in IDC stroma)", "Lower than typical")))
message(sprintf("  %-18s: %5.1f%% - %s", "Cancer Epithelial", mean_props["Cancer Epithelial"],
                ifelse(mean_props["Cancer Epithelial"] > 5, "EXPECTED (tumor cells)", "May be low")))
message(sprintf("  %-18s: %5.1f%% - %s", "T-cells", mean_props["T-cells"],
                ifelse(mean_props["T-cells"] > 5, "EXPECTED (TILs present)", "Low immune infiltration")))
message(sprintf("  %-18s: %5.1f%% - %s", "Myeloid", mean_props["Myeloid"],
                ifelse(mean_props["Myeloid"] > 5, "EXPECTED (TAMs present)", "Low myeloid")))
message(sprintf("  %-18s: %5.1f%% - %s", "Endothelial", mean_props["Endothelial"],
                ifelse(mean_props["Endothelial"] > 3, "EXPECTED (tumor vasculature)", "May be low")))

# Check for biologically implausible results
issues <- c()
if (mean_props["CAFs"] < 5) issues <- c(issues, "Very low CAFs unusual for IDC")
if (mean_props["Cancer Epithelial"] < 1) issues <- c(issues, "Very low cancer cells - check sample")
if (sum(mean_props[c("T-cells", "B-cells", "Myeloid")]) < 5) {
  issues <- c(issues, "Very low immune infiltration")
}

if (length(issues) == 0) {
  message("\nAssessment: BIOLOGICALLY PLAUSIBLE")
  bio_score <- 3
} else {
  message("\nPotential issues:")
  for (issue in issues) message(paste("  -", issue))
  bio_score <- 1
}

# =============================================================================
# OVERALL QUALITY SCORE
# =============================================================================
message("\n", strrep("=", 60))
message("OVERALL QUALITY SUMMARY")
message(strrep("=", 60))

total_score <- convergence_score + normalization_score + separation_score +
               confidence_score + bio_score
max_score <- 15

message(sprintf("\n  Convergence:     %d/3", convergence_score))
message(sprintf("  Normalization:   %d/3", normalization_score))
message(sprintf("  Cell Separation: %d/3", separation_score))
message(sprintf("  Confidence:      %d/3", confidence_score))
message(sprintf("  Bio Plausibility:%d/3", bio_score))
message(sprintf("\n  TOTAL SCORE: %d/%d (%.0f%%)", total_score, max_score,
                total_score/max_score*100))

if (total_score >= 12) {
  message("\n  OVERALL: HIGH QUALITY DECONVOLUTION")
} else if (total_score >= 9) {
  message("\n  OVERALL: ACCEPTABLE QUALITY")
} else {
  message("\n  OVERALL: QUALITY CONCERNS - Review results carefully")
}

# =============================================================================
# Generate Quality Assessment Figures
# =============================================================================
message("\n", strrep("=", 60))
message("GENERATING QUALITY FIGURES")
message(strrep("=", 60))

# Figure 1: Weight distribution per cell type
weights_long <- as.data.frame(weights) %>%
  mutate(spot = rownames(weights)) %>%
  pivot_longer(-spot, names_to = "cell_type", values_to = "proportion")

p1 <- ggplot(weights_long, aes(x = proportion)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  facet_wrap(~cell_type, scales = "free_y", ncol = 3) +
  theme_minimal() +
  labs(title = "Distribution of Cell Type Proportions",
       x = "Proportion", y = "Count") +
  theme(strip.text = element_text(size = 8))

# Figure 2: Confidence distribution
confidence_df <- data.frame(
  max_weight = max_weights,
  margin = confidence_margin
)

p2 <- ggplot(confidence_df, aes(x = max_weight)) +
  geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0.3, linetype = "dashed", color = "orange") +
  theme_minimal() +
  labs(title = "Classification Confidence",
       subtitle = "Red: 50% threshold, Orange: 30% threshold",
       x = "Dominant Cell Type Proportion", y = "Number of Spots") +
  annotate("text", x = 0.85, y = Inf, vjust = 2,
           label = paste0(round(sum(max_weights > 0.5)/length(max_weights)*100), "% > 50%"))

# Figure 3: Cell type correlation heatmap
cor_long <- as.data.frame(cor_matrix) %>%
  mutate(type1 = rownames(cor_matrix)) %>%
  pivot_longer(-type1, names_to = "type2", values_to = "correlation")

p3 <- ggplot(cor_long, aes(x = type1, y = type2, fill = correlation)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 2.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8)) +
  labs(title = "Cell Type Correlation", x = "", y = "")

# Figure 4: Spatial distribution of confidence
common_spots <- intersect(rownames(weights), rownames(coords))
plot_df <- data.frame(
  x = coords[common_spots, 1],
  y = coords[common_spots, 2],
  confidence = max_weights[common_spots]
)

p4 <- ggplot(plot_df, aes(x = x, y = y, color = confidence)) +
  geom_point(size = 0.8) +
  scale_color_viridis_c(option = "plasma", name = "Confidence") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.title = element_blank()) +
  coord_fixed() +
  labs(title = "Spatial Distribution of Classification Confidence")

# Combine and save
combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = "Deconvolution Quality Assessment",
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

ggsave(file.path(FIGURE_DIR, "quality_assessment.png"), combined,
       width = 14, height = 12, dpi = 300)

message("Saved: quality_assessment.png")

# =============================================================================
# Save Quality Report
# =============================================================================
quality_report <- data.frame(
  Metric = c("Convergence Rate", "Weight Normalization", "Cell Type Separation",
             "Classification Confidence", "Biological Plausibility", "TOTAL"),
  Score = c(convergence_score, normalization_score, separation_score,
            confidence_score, bio_score, total_score),
  Max_Score = c(3, 3, 3, 3, 3, 15),
  Value = c(paste0(round(convergence_rate, 1), "%"),
            round(mean(weight_sums), 4),
            round(mean(off_diag), 3),
            round(mean(max_weights), 3),
            "See details",
            paste0(round(total_score/max_score*100), "%"))
)

write.csv(quality_report, file.path(OUTPUT_DIR, "quality_report.csv"), row.names = FALSE)
message("Saved: quality_report.csv")

message("\n=== Quality assessment complete ===")
