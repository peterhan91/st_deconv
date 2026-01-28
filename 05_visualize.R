# =============================================================================
# 05_visualize.R - Visualization of Deconvolution Results
# ST Deconvolution Pipeline for IDC Data using RCTD
# =============================================================================

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(viridis)
library(Matrix)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
OUTPUT_DIR <- "output/data"
FIGURE_DIR <- "output/figures"

# Create figure directory if it doesn't exist
if (!dir.exists(FIGURE_DIR)) dir.create(FIGURE_DIR, recursive = TRUE)

# Plot parameters
POINT_SIZE <- 1.5
FIGURE_WIDTH <- 10
FIGURE_HEIGHT <- 8

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
message("=== Loading data ===")

# Load cell type weights
weights_file <- file.path(OUTPUT_DIR, "cell_type_weights.rds")
if (!file.exists(weights_file)) {
  stop("Cell type weights not found. Please run 04_run_rctd.R first.")
}
weights <- readRDS(weights_file)

# Convert to regular matrix/data.frame if it's a sparse matrix
if (inherits(weights, "Matrix") || inherits(weights, "dgeMatrix")) {
  weights <- as.matrix(weights)
}
weights <- as.data.frame(weights)

message(paste("Loaded weights for", nrow(weights), "spots and", ncol(weights), "cell types"))

# Load coordinates
coords_file <- file.path(OUTPUT_DIR, "spatial_coords.csv")
if (file.exists(coords_file)) {
  coords <- read.csv(coords_file, row.names = 1)
  message(paste("Loaded coordinates for", nrow(coords), "spots"))
} else {
  stop("Coordinates file not found.")
}

# Ensure matching spots
common_spots <- intersect(rownames(weights), rownames(coords))
weights <- weights[common_spots, ]
coords <- coords[common_spots, ]

# Combine data for plotting
plot_data <- cbind(coords, weights)

# Get cell type names
cell_types <- colnames(weights)
message(paste("Cell types:", paste(cell_types, collapse = ", ")))

# -----------------------------------------------------------------------------
# Define color palette
# -----------------------------------------------------------------------------
# Create a color palette for cell types
n_types <- length(cell_types)
if (n_types <= 8) {
  type_colors <- brewer.pal(max(3, n_types), "Set2")[1:n_types]
} else if (n_types <= 12) {
  type_colors <- brewer.pal(n_types, "Set3")
} else {
  type_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_types)
}
names(type_colors) <- cell_types

# -----------------------------------------------------------------------------
# Plot 1: Spatial distribution of each cell type
# -----------------------------------------------------------------------------
message("\n=== Creating spatial distribution plots ===")

# Create individual plots for each cell type
spatial_plots <- list()

for (ct in cell_types) {
  ct_clean <- gsub("[^a-zA-Z0-9_]", "_", ct)

  p <- ggplot(plot_data, aes_string(x = "x", y = "y", color = paste0("`", ct, "`"))) +
    geom_point(size = POINT_SIZE) +
    scale_color_viridis(name = "Proportion", limits = c(0, 1), option = "plasma") +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.position = "right"
    ) +
    coord_fixed() +
    ggtitle(ct)

  spatial_plots[[ct]] <- p
}

# Combine into a grid
n_cols <- min(3, ceiling(sqrt(n_types)))
n_rows <- ceiling(n_types / n_cols)

combined_spatial <- wrap_plots(spatial_plots, ncol = n_cols)

# Save combined plot
ggsave(
  file.path(FIGURE_DIR, "spatial_celltypes_all.png"),
  combined_spatial,
  width = FIGURE_WIDTH * n_cols / 3,
  height = FIGURE_HEIGHT * n_rows / 3,
  dpi = 300
)
message("Saved: spatial_celltypes_all.png")

# Save individual plots
for (ct in cell_types) {
  ct_clean <- gsub("[^a-zA-Z0-9_]", "_", ct)
  ggsave(
    file.path(FIGURE_DIR, paste0("spatial_", ct_clean, ".png")),
    spatial_plots[[ct]],
    width = 6, height = 5, dpi = 300
  )
}
message("Saved individual cell type spatial plots")

# -----------------------------------------------------------------------------
# Plot 2: Dominant cell type per spot
# -----------------------------------------------------------------------------
message("\n=== Creating dominant cell type plot ===")

# Find dominant cell type for each spot
dominant_type <- apply(weights, 1, function(x) {
  cell_types[which.max(x)]
})
plot_data$dominant_type <- factor(dominant_type, levels = cell_types)

# Calculate dominance score (proportion of top cell type)
plot_data$dominance_score <- apply(weights, 1, max)

# Plot dominant cell types
p_dominant <- ggplot(plot_data, aes(x = x, y = y, color = dominant_type)) +
  geom_point(size = POINT_SIZE) +
  scale_color_manual(values = type_colors, name = "Cell Type") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  coord_fixed() +
  ggtitle("Dominant Cell Type per Spot")

ggsave(
  file.path(FIGURE_DIR, "spatial_dominant_celltype.png"),
  p_dominant,
  width = 8, height = 6, dpi = 300
)
message("Saved: spatial_dominant_celltype.png")

# Plot dominance score
p_dominance <- ggplot(plot_data, aes(x = x, y = y, color = dominance_score)) +
  geom_point(size = POINT_SIZE) +
  scale_color_viridis(name = "Dominance\nScore", limits = c(0, 1), option = "magma") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  coord_fixed() +
  ggtitle("Cell Type Dominance Score")

ggsave(
  file.path(FIGURE_DIR, "spatial_dominance_score.png"),
  p_dominance,
  width = 7, height = 6, dpi = 300
)
message("Saved: spatial_dominance_score.png")

# -----------------------------------------------------------------------------
# Plot 3: Cell type proportion boxplot
# -----------------------------------------------------------------------------
message("\n=== Creating proportion boxplot ===")

# Reshape data for boxplot
weights_long <- weights %>%
  as.data.frame() %>%
  mutate(spot = rownames(.)) %>%
  pivot_longer(cols = -spot, names_to = "cell_type", values_to = "proportion")

# Calculate mean proportions for ordering
mean_props <- weights_long %>%
  group_by(cell_type) %>%
  summarize(mean_prop = mean(proportion)) %>%
  arrange(desc(mean_prop))

weights_long$cell_type <- factor(weights_long$cell_type, levels = mean_props$cell_type)

# Create boxplot
p_boxplot <- ggplot(weights_long, aes(x = cell_type, y = proportion, fill = cell_type)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = type_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    x = "Cell Type",
    y = "Proportion",
    title = "Cell Type Proportions Across All Spots"
  )

ggsave(
  file.path(FIGURE_DIR, "celltype_proportions_boxplot.png"),
  p_boxplot,
  width = max(8, n_types * 0.8), height = 6, dpi = 300
)
message("Saved: celltype_proportions_boxplot.png")

# -----------------------------------------------------------------------------
# Plot 4: Cell type proportion barplot (mean)
# -----------------------------------------------------------------------------
message("\n=== Creating mean proportion barplot ===")

p_barplot <- ggplot(mean_props, aes(x = reorder(cell_type, -mean_prop), y = mean_prop, fill = cell_type)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  scale_fill_manual(values = type_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    x = "Cell Type",
    y = "Mean Proportion",
    title = "Mean Cell Type Proportions"
  ) +
  geom_text(aes(label = sprintf("%.1f%%", mean_prop * 100)),
            vjust = -0.5, size = 3)

ggsave(
  file.path(FIGURE_DIR, "celltype_mean_proportions.png"),
  p_barplot,
  width = max(8, n_types * 0.8), height = 6, dpi = 300
)
message("Saved: celltype_mean_proportions.png")

# -----------------------------------------------------------------------------
# Plot 5: Cell type correlation heatmap
# -----------------------------------------------------------------------------
message("\n=== Creating correlation heatmap ===")

# Calculate correlation between cell type proportions
cor_matrix <- cor(weights)

# Convert to long format for ggplot
cor_long <- cor_matrix %>%
  as.data.frame() %>%
  mutate(type1 = rownames(.)) %>%
  pivot_longer(cols = -type1, names_to = "type2", values_to = "correlation")

cor_long$type1 <- factor(cor_long$type1, levels = cell_types)
cor_long$type2 <- factor(cor_long$type2, levels = rev(cell_types))

# Create heatmap
p_heatmap <- ggplot(cor_long, aes(x = type1, y = type2, fill = correlation)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_blank()
  ) +
  ggtitle("Cell Type Proportion Correlations")

ggsave(
  file.path(FIGURE_DIR, "celltype_correlation_heatmap.png"),
  p_heatmap,
  width = max(8, n_types * 0.7), height = max(7, n_types * 0.6), dpi = 300
)
message("Saved: celltype_correlation_heatmap.png")

# -----------------------------------------------------------------------------
# Plot 6: Pie chart summary
# -----------------------------------------------------------------------------
message("\n=== Creating pie chart ===")

pie_data <- mean_props %>%
  mutate(
    percentage = mean_prop * 100,
    label = sprintf("%s\n(%.1f%%)", cell_type, percentage)
  )

p_pie <- ggplot(pie_data, aes(x = "", y = mean_prop, fill = cell_type)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = type_colors, name = "Cell Type") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  ggtitle("Overall Cell Type Composition")

ggsave(
  file.path(FIGURE_DIR, "celltype_pie_chart.png"),
  p_pie,
  width = 8, height = 6, dpi = 300
)
message("Saved: celltype_pie_chart.png")

# -----------------------------------------------------------------------------
# Plot 7: Stacked composition along spatial axis
# -----------------------------------------------------------------------------
message("\n=== Creating spatial composition plot ===")

# Sort spots by x coordinate and create bins
plot_data_sorted <- plot_data %>%
  arrange(x) %>%
  mutate(x_bin = cut(x, breaks = 20, labels = FALSE))

# Calculate mean proportions per bin
composition_by_bin <- plot_data_sorted %>%
  select(x_bin, all_of(cell_types)) %>%
  group_by(x_bin) %>%
  summarize(across(everything(), mean)) %>%
  pivot_longer(cols = -x_bin, names_to = "cell_type", values_to = "proportion")

composition_by_bin$cell_type <- factor(composition_by_bin$cell_type, levels = rev(cell_types))

p_stacked <- ggplot(composition_by_bin, aes(x = x_bin, y = proportion, fill = cell_type)) +
  geom_area(position = "stack", alpha = 0.8) +
  scale_fill_manual(values = type_colors, name = "Cell Type") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "right"
  ) +
  labs(
    x = "Spatial Position (X-axis bins)",
    y = "Proportion",
    title = "Cell Type Composition Along Spatial Axis"
  )

ggsave(
  file.path(FIGURE_DIR, "spatial_composition_stacked.png"),
  p_stacked,
  width = 10, height = 6, dpi = 300
)
message("Saved: spatial_composition_stacked.png")

# -----------------------------------------------------------------------------
# Summary statistics table
# -----------------------------------------------------------------------------
message("\n=== Generating summary statistics ===")

summary_stats <- data.frame(
  cell_type = cell_types,
  mean_proportion = colMeans(weights),
  sd_proportion = apply(weights, 2, sd),
  min_proportion = apply(weights, 2, min),
  max_proportion = apply(weights, 2, max),
  median_proportion = apply(weights, 2, median),
  spots_with_any = colSums(weights > 0.01),
  spots_dominant = table(dominant_type)[cell_types]
)

summary_stats <- summary_stats %>%
  arrange(desc(mean_proportion))

# Save summary statistics
write.csv(summary_stats, file.path(OUTPUT_DIR, "deconvolution_summary.csv"), row.names = FALSE)
message(paste("Saved: deconvolution_summary.csv"))

# Print summary
message("\n=== Deconvolution Summary ===")
print(summary_stats)

# -----------------------------------------------------------------------------
# Combined figure for publication
# -----------------------------------------------------------------------------
message("\n=== Creating combined figure ===")

# Combine key plots
combined_fig <- (p_dominant + p_boxplot) / (p_heatmap + p_pie) +
  plot_annotation(
    title = "ST Deconvolution Results - IDC Sample",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )

ggsave(
  file.path(FIGURE_DIR, "deconvolution_summary_figure.png"),
  combined_fig,
  width = 14, height = 12, dpi = 300
)
message("Saved: deconvolution_summary_figure.png")

message("\n=== Visualization complete! ===")
message(paste("All figures saved to:", FIGURE_DIR))
message("\nGenerated figures:")
list.files(FIGURE_DIR, pattern = "\\.png$")
