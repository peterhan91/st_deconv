# =============================================================================
# 05_visualize_bulk.R - Visualize MuSiC Deconvolution Results
# Bulk Deconvolution Pipeline using MuSiC
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(patchwork)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
OUTPUT_DIR <- "output"
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")

if (!dir.exists(FIGURE_DIR)) dir.create(FIGURE_DIR, recursive = TRUE)

# Color palette for cell types
celltype_colors <- c(
  "B-cells" = "#E41A1C",
  "CAFs" = "#377EB8",
  "Cancer Epithelial" = "#4DAF4A",
  "Endothelial" = "#984EA3",
  "Myeloid" = "#FF7F00",
  "Normal Epithelial" = "#FFFF33",
  "Plasmablasts" = "#A65628",
  "PVL" = "#F781BF",
  "T-cells" = "#999999"
)

# -----------------------------------------------------------------------------
# Load Results
# -----------------------------------------------------------------------------
message("=== Loading MuSiC results ===")

proportions <- readRDS(file.path(OUTPUT_DIR, "data", "music_proportions.rds"))
proportions <- as.data.frame(proportions)

message(paste("Loaded:", nrow(proportions), "samples x", ncol(proportions), "cell types"))

# Load R-squared if available
r2_file <- file.path(OUTPUT_DIR, "data", "model_fit_rsquared.csv")
if (file.exists(r2_file)) {
  r2_data <- read.csv(r2_file)
}

# -----------------------------------------------------------------------------
# Figure 1: Mean Cell Type Proportions (Bar Chart)
# -----------------------------------------------------------------------------
message("\n=== Generating figures ===")

mean_props <- colMeans(proportions) * 100
se_props <- apply(proportions, 2, sd) / sqrt(nrow(proportions)) * 100

summary_df <- data.frame(
  CellType = names(mean_props),
  Mean = mean_props,
  SE = se_props
) %>%
  arrange(desc(Mean)) %>%
  mutate(CellType = factor(CellType, levels = CellType))

p1 <- ggplot(summary_df, aes(x = CellType, y = Mean, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.3) +
  scale_fill_manual(values = celltype_colors, guide = "none") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  labs(title = "Mean Cell Type Proportions (CPTAC BRCA)",
       x = "", y = "Mean Proportion (%)")

ggsave(file.path(FIGURE_DIR, "mean_proportions_barplot.png"), p1,
       width = 10, height = 6, dpi = 300)
message("Saved: mean_proportions_barplot.png")

# -----------------------------------------------------------------------------
# Figure 2: Cell Type Proportions Boxplot
# -----------------------------------------------------------------------------
props_long <- proportions %>%
  mutate(Sample = rownames(proportions)) %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "Proportion") %>%
  mutate(CellType = factor(CellType, levels = summary_df$CellType))

p2 <- ggplot(props_long, aes(x = CellType, y = Proportion * 100, fill = CellType)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = celltype_colors, guide = "none") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  labs(title = "Cell Type Proportion Distribution",
       x = "", y = "Proportion (%)")

ggsave(file.path(FIGURE_DIR, "proportions_boxplot.png"), p2,
       width = 10, height = 6, dpi = 300)
message("Saved: proportions_boxplot.png")

# -----------------------------------------------------------------------------
# Figure 3: Sample Composition Heatmap
# -----------------------------------------------------------------------------
# Prepare matrix for heatmap
prop_matrix <- as.matrix(proportions)

# Order samples by dominant cell type
dominant_type <- apply(prop_matrix, 1, which.max)
sample_order <- order(dominant_type, -apply(prop_matrix, 1, max))
prop_matrix_ordered <- prop_matrix[sample_order, ]

# Create annotation for dominant type
dominant_names <- colnames(prop_matrix)[dominant_type[sample_order]]
sample_annotation <- data.frame(
  DominantType = dominant_names,
  row.names = rownames(prop_matrix_ordered)
)

# Heatmap
png(file.path(FIGURE_DIR, "composition_heatmap.png"), width = 12, height = 10, units = "in", res = 300)
pheatmap(
  t(prop_matrix_ordered),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  color = colorRampPalette(c("white", "blue", "darkblue"))(100),
  main = "Sample Cell Type Composition",
  fontsize = 10
)
dev.off()
message("Saved: composition_heatmap.png")

# -----------------------------------------------------------------------------
# Figure 4: Cell Type Correlation Heatmap
# -----------------------------------------------------------------------------
cor_matrix <- cor(proportions)

png(file.path(FIGURE_DIR, "celltype_correlation.png"), width = 8, height = 7, units = "in", res = 300)
pheatmap(
  cor_matrix,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 8,
  main = "Cell Type Proportion Correlations"
)
dev.off()
message("Saved: celltype_correlation.png")

# -----------------------------------------------------------------------------
# Figure 5: Stacked Bar Plot (Sample Composition)
# -----------------------------------------------------------------------------
# Show first 50 samples if many
n_show <- min(50, nrow(proportions))
props_subset <- proportions[1:n_show, ]

props_stacked <- props_subset %>%
  mutate(Sample = factor(rownames(props_subset), levels = rownames(props_subset))) %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "Proportion")

p5 <- ggplot(props_stacked, aes(x = Sample, y = Proportion * 100, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    legend.position = "right"
  ) +
  labs(title = paste("Sample Composition (first", n_show, "samples)"),
       x = "Sample", y = "Proportion (%)", fill = "Cell Type")

ggsave(file.path(FIGURE_DIR, "sample_composition_stacked.png"), p5,
       width = 14, height = 6, dpi = 300)
message("Saved: sample_composition_stacked.png")

# -----------------------------------------------------------------------------
# Figure 6: Pie Chart (Overall Composition)
# -----------------------------------------------------------------------------
pie_df <- data.frame(
  CellType = names(mean_props),
  Proportion = mean_props
) %>%
  arrange(desc(Proportion)) %>%
  mutate(
    CellType = factor(CellType, levels = CellType),
    ypos = cumsum(Proportion) - 0.5 * Proportion,
    label = paste0(round(Proportion, 1), "%")
  )

p6 <- ggplot(pie_df, aes(x = "", y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = celltype_colors) +
  theme_void() +
  theme(legend.position = "right") +
  labs(title = "Overall Cell Type Composition", fill = "Cell Type")

ggsave(file.path(FIGURE_DIR, "composition_pie.png"), p6,
       width = 8, height = 6, dpi = 300)
message("Saved: composition_pie.png")

# -----------------------------------------------------------------------------
# Figure 7: Model Fit (R-squared) Distribution
# -----------------------------------------------------------------------------
if (exists("r2_data")) {
  p7 <- ggplot(r2_data, aes(x = R_squared)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = mean(r2_data$R_squared), color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = "Model Fit Quality (R-squared)",
         subtitle = paste("Mean R² =", round(mean(r2_data$R_squared), 3)),
         x = "R-squared", y = "Number of Samples")

  ggsave(file.path(FIGURE_DIR, "model_fit_rsquared.png"), p7,
         width = 8, height = 6, dpi = 300)
  message("Saved: model_fit_rsquared.png")
}

# -----------------------------------------------------------------------------
# Figure 8: Individual Cell Type Distributions
# -----------------------------------------------------------------------------
p8 <- ggplot(props_long, aes(x = Proportion * 100)) +
  geom_histogram(aes(fill = CellType), bins = 30, alpha = 0.7) +
  facet_wrap(~CellType, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = celltype_colors, guide = "none") +
  theme_minimal() +
  theme(strip.text = element_text(size = 9)) +
  labs(title = "Distribution of Cell Type Proportions",
       x = "Proportion (%)", y = "Count")

ggsave(file.path(FIGURE_DIR, "celltype_distributions.png"), p8,
       width = 12, height = 8, dpi = 300)
message("Saved: celltype_distributions.png")

# -----------------------------------------------------------------------------
# Figure 9: Summary Figure (Combined)
# -----------------------------------------------------------------------------
summary_fig <- (p1 + p2) / (p6 + p8) +
  plot_annotation(
    title = "MuSiC Bulk Deconvolution Results - CPTAC BRCA",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

ggsave(file.path(FIGURE_DIR, "deconvolution_summary.png"), summary_fig,
       width = 16, height = 12, dpi = 300)
message("Saved: deconvolution_summary.png")

# -----------------------------------------------------------------------------
# Figure 10: Comparison with Expected Values (if applicable)
# -----------------------------------------------------------------------------
# Create a comparison with literature values for breast cancer
literature_props <- data.frame(
  CellType = c("CAFs", "Cancer Epithelial", "T-cells", "Myeloid", "Endothelial",
               "B-cells", "Normal Epithelial", "PVL", "Plasmablasts"),
  Literature = c(25, 20, 15, 12, 10, 5, 5, 5, 3),  # Approximate expected values
  stringsAsFactors = FALSE
)

comparison_df <- data.frame(
  CellType = names(mean_props),
  Estimated = as.numeric(mean_props)
) %>%
  left_join(literature_props, by = "CellType") %>%
  pivot_longer(c(Estimated, Literature), names_to = "Source", values_to = "Proportion") %>%
  mutate(CellType = factor(CellType, levels = summary_df$CellType))

p10 <- ggplot(comparison_df, aes(x = CellType, y = Proportion, fill = Source)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Estimated" = "steelblue", "Literature" = "gray60")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Estimated vs Expected Cell Type Proportions",
       subtitle = "Literature values are approximate ranges for breast cancer",
       x = "", y = "Proportion (%)")

ggsave(file.path(FIGURE_DIR, "comparison_literature.png"), p10,
       width = 10, height = 6, dpi = 300)
message("Saved: comparison_literature.png")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
message("\n", strrep("=", 60))
message("VISUALIZATION COMPLETE")
message(strrep("=", 60))
message("\nGenerated figures in output/figures/:")
message("  - mean_proportions_barplot.png")
message("  - proportions_boxplot.png")
message("  - composition_heatmap.png")
message("  - celltype_correlation.png")
message("  - sample_composition_stacked.png")
message("  - composition_pie.png")
message("  - model_fit_rsquared.png")
message("  - celltype_distributions.png")
message("  - deconvolution_summary.png")
message("  - comparison_literature.png")

message("\n=== Bulk deconvolution pipeline complete! ===")
