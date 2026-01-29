# =============================================================================
# 11_pseudobulk_validation.R - Validate Deconvolution Methods with Pseudo-bulk
# =============================================================================

library(Biobase)
library(SingleCellExperiment)
library(MuSiC)
library(BisqueRNA)
library(dtangle)
library(BayesPrism)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(patchwork)

set.seed(42)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
OUTPUT_DIR <- "output"
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
DATA_DIR <- file.path(OUTPUT_DIR, "data")
VALIDATION_DIR <- file.path(OUTPUT_DIR, "validation")

dir.create(VALIDATION_DIR, showWarnings = FALSE, recursive = TRUE)

N_PSEUDOBULK <- 50
CELLS_PER_SAMPLE <- 500

# -----------------------------------------------------------------------------
# Load scRNA-seq Reference
# -----------------------------------------------------------------------------
message("=== Loading scRNA-seq reference ===")

sc_eset <- readRDS(file.path(DATA_DIR, "sc_reference_eset_filtered.rds"))
message(paste("Reference:", nrow(sc_eset), "genes x", ncol(sc_eset), "cells"))

cell_types <- unique(pData(sc_eset)$cellType)
message(paste("Cell types:", paste(cell_types, collapse = ", ")))

cells_per_type <- table(pData(sc_eset)$cellType)
message("\nCells per type:")
print(cells_per_type)

# -----------------------------------------------------------------------------
# Generate Ground Truth Proportions
# -----------------------------------------------------------------------------
message("\n=== Generating ground truth proportions ===")

generate_random_proportions <- function(n_samples, cell_types, sparsity = 0.3) {
  n_types <- length(cell_types)
  props <- matrix(0, nrow = n_samples, ncol = n_types)
  colnames(props) <- cell_types

  for (i in 1:n_samples) {
    raw_props <- runif(n_types)
    mask <- runif(n_types) > sparsity
    raw_props <- raw_props * mask
    if (sum(mask) < 2) {
      top2 <- order(runif(n_types), decreasing = TRUE)[1:2]
      raw_props[top2] <- runif(2)
    }
    props[i, ] <- raw_props / sum(raw_props)
  }
  return(props)
}

ground_truth <- generate_random_proportions(N_PSEUDOBULK, cell_types)
rownames(ground_truth) <- paste0("PseudoBulk_", 1:N_PSEUDOBULK)

message(paste("Generated", N_PSEUDOBULK, "pseudo-bulk samples"))
message("\nGround truth summary (mean %):")
print(round(colMeans(ground_truth) * 100, 2))

# -----------------------------------------------------------------------------
# Create Pseudo-bulk Samples
# -----------------------------------------------------------------------------
message("\n=== Creating pseudo-bulk samples ===")

create_pseudobulk <- function(sc_eset, proportions, cells_per_sample) {
  n_samples <- nrow(proportions)
  cell_types <- colnames(proportions)
  expr_mat <- exprs(sc_eset)
  cell_info <- pData(sc_eset)

  pseudobulk <- matrix(0, nrow = nrow(expr_mat), ncol = n_samples)
  rownames(pseudobulk) <- rownames(expr_mat)
  colnames(pseudobulk) <- rownames(proportions)

  for (i in 1:n_samples) {
    sample_cells <- c()
    for (ct in cell_types) {
      n_cells <- round(proportions[i, ct] * cells_per_sample)
      if (n_cells > 0) {
        available_cells <- which(cell_info$cellType == ct)
        sampled <- sample(available_cells, n_cells, replace = TRUE)
        sample_cells <- c(sample_cells, sampled)
      }
    }
    if (length(sample_cells) > 0) {
      pseudobulk[, i] <- rowSums(expr_mat[, sample_cells, drop = FALSE])
    }
  }
  return(pseudobulk)
}

pseudobulk_counts <- create_pseudobulk(sc_eset, ground_truth, CELLS_PER_SAMPLE)
message(paste("Pseudo-bulk matrix:", nrow(pseudobulk_counts), "x", ncol(pseudobulk_counts)))

pseudobulk_eset <- ExpressionSet(
  assayData = pseudobulk_counts,
  phenoData = AnnotatedDataFrame(data.frame(
    row.names = colnames(pseudobulk_counts),
    SampleID = colnames(pseudobulk_counts)
  ))
)

# -----------------------------------------------------------------------------
# Run Deconvolution Methods
# -----------------------------------------------------------------------------
message("\n=== Running deconvolution methods ===")

results_list <- list()

# 1. MuSiC
message("\n--- MuSiC ---")
tryCatch({
  sc_sce <- SingleCellExperiment(
    assays = list(counts = exprs(sc_eset)),
    colData = pData(sc_eset)
  )
  music_result <- music_prop(
    bulk.mtx = pseudobulk_counts,
    sc.sce = sc_sce,
    clusters = "cellType",
    samples = "SubjectName",
    verbose = FALSE
  )
  results_list$MuSiC <- music_result$Est.prop.weighted
  message("MuSiC completed")
}, error = function(e) message(paste("MuSiC failed:", e$message)))

# 2. NNLS
message("\n--- NNLS ---")
tryCatch({
  ref_matrix <- sapply(cell_types, function(ct) {
    cells <- which(pData(sc_eset)$cellType == ct)
    rowMeans(exprs(sc_eset)[, cells, drop = FALSE])
  })
  common_genes <- intersect(rownames(pseudobulk_counts), rownames(ref_matrix))
  bulk_sub <- pseudobulk_counts[common_genes, ]
  ref_sub <- ref_matrix[common_genes, ]

  nnls_props <- matrix(0, nrow = ncol(bulk_sub), ncol = length(cell_types))
  rownames(nnls_props) <- colnames(bulk_sub)
  colnames(nnls_props) <- cell_types

  for (i in 1:ncol(bulk_sub)) {
    fit <- nnls::nnls(ref_sub, bulk_sub[, i])
    props <- fit$x
    nnls_props[i, ] <- props / sum(props)
  }
  results_list$NNLS <- nnls_props
  message("NNLS completed")
}, error = function(e) message(paste("NNLS failed:", e$message)))

# 3. BisqueRNA
message("\n--- BisqueRNA ---")
tryCatch({
  bisque_result <- ReferenceBasedDecomposition(
    bulk.eset = pseudobulk_eset,
    sc.eset = sc_eset,
    cell.types = "cellType",
    subject.names = "SubjectName",
    use.overlap = FALSE
  )
  results_list$BisqueRNA <- t(bisque_result$bulk.props)
  message("BisqueRNA completed")
}, error = function(e) message(paste("BisqueRNA failed:", e$message)))

# 4. dtangle
message("\n--- dtangle ---")
tryCatch({
  ref_profiles <- sapply(cell_types, function(ct) {
    cells <- which(pData(sc_eset)$cellType == ct)
    rowMeans(exprs(sc_eset)[, cells, drop = FALSE])
  })
  common_genes <- intersect(rownames(pseudobulk_counts), rownames(ref_profiles))
  combined <- cbind(ref_profiles[common_genes, ], pseudobulk_counts[common_genes, ])
  combined_log <- log2(combined + 1)

  dtangle_result <- dtangle(
    Y = t(combined_log),
    pure_samples = lapply(1:length(cell_types), function(i) i),
    n_markers = 100
  )
  bulk_indices <- (length(cell_types) + 1):ncol(combined)
  dtangle_props <- dtangle_result$estimates[bulk_indices, ]
  rownames(dtangle_props) <- colnames(pseudobulk_counts)
  results_list$dtangle <- dtangle_props
  message("dtangle completed")
}, error = function(e) message(paste("dtangle failed:", e$message)))

# 5. BayesPrism
message("\n--- BayesPrism ---")
tryCatch({
  sc_mtx <- exprs(sc_eset)
  cell_type_labels <- as.character(pData(sc_eset)$cellType)

  common_genes <- intersect(rownames(pseudobulk_counts), rownames(sc_mtx))
  sc_sub <- sc_mtx[common_genes, ]
  bulk_sub <- pseudobulk_counts[common_genes, ]

  # BayesPrism needs cells x genes and samples x genes
  sc_t <- t(sc_sub)
  bulk_t <- t(bulk_sub)

  # Scale if normalized
  if (max(sc_t) < 100) {
    sc_t <- round(sc_t * 100)
    bulk_t <- round(bulk_t * 100)
  }

  # Filter low expression genes
  gene_sums <- colSums(sc_t)
  keep_genes <- gene_sums > ncol(sc_t) * 0.01
  sc_t_filt <- sc_t[, keep_genes]
  bulk_t_filt <- bulk_t[, keep_genes]

  # Create prism reference
  myPrism <- new.prism(
    reference = sc_t_filt,
    mixture = bulk_t_filt,
    input.type = "count.matrix",
    cell.type.labels = cell_type_labels,
    cell.state.labels = cell_type_labels,
    key = NULL,
    outlier.cut = 0.01,
    outlier.fraction = 0.1
  )

  # Run BayesPrism
  bp_result <- run.prism(prism = myPrism, n.cores = 4)

  # Extract fractions
  theta <- get.fraction(bp = bp_result, which.theta = "final", state.or.type = "type")
  results_list$BayesPrism <- theta
  message("BayesPrism completed")
}, error = function(e) message(paste("BayesPrism failed:", e$message)))

message(paste("\nCompleted:", paste(names(results_list), collapse = ", ")))

# -----------------------------------------------------------------------------
# Evaluate Performance
# -----------------------------------------------------------------------------
message("\n=== Evaluating performance ===")

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

evaluation_results <- lapply(results_list, function(est) evaluate_method(est, ground_truth))

summary_table <- data.frame(
  Method = names(evaluation_results),
  Correlation = sapply(evaluation_results, function(x) round(x$overall_cor, 4)),
  RMSE = sapply(evaluation_results, function(x) round(x$overall_rmse, 4)),
  MAE = sapply(evaluation_results, function(x) round(x$overall_mae, 4))
) %>%
  filter(!is.na(Correlation)) %>%
  arrange(desc(Correlation))

message("\n=== Method Performance Summary ===")
print(summary_table)

# Remove methods with NA results
evaluation_results <- evaluation_results[!sapply(evaluation_results, function(x) is.na(x$overall_cor))]
results_list <- results_list[names(evaluation_results)]

celltype_cors_df <- do.call(rbind, lapply(names(evaluation_results), function(m) {
  data.frame(
    Method = m,
    CellType = names(evaluation_results[[m]]$celltype_cors),
    Correlation = evaluation_results[[m]]$celltype_cors
  )
}))

message("\n=== Per Cell Type Correlations ===")
celltype_cors_wide <- celltype_cors_df %>%
  pivot_wider(names_from = Method, values_from = Correlation)
print(celltype_cors_wide)

# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------
message("\n=== Generating validation figures ===")

# Figure 1: Method accuracy
p1 <- ggplot(summary_table, aes(x = reorder(Method, Correlation), y = Correlation)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", Correlation)), hjust = -0.1) +
  coord_flip() + ylim(0, 1) + theme_minimal() +
  labs(title = "Deconvolution Accuracy (Pseudo-bulk Validation)",
       subtitle = paste("n =", N_PSEUDOBULK, "pseudo-bulk samples"),
       x = "", y = "Correlation with Ground Truth")

ggsave(file.path(VALIDATION_DIR, "method_accuracy_barplot.png"), p1, width = 8, height = 5, dpi = 300)

# Figure 2: Scatter plots
scatter_data <- do.call(rbind, lapply(names(results_list), function(m) {
  est <- results_list[[m]]
  common_samples <- intersect(rownames(est), rownames(ground_truth))
  common_types <- intersect(colnames(est), colnames(ground_truth))
  data.frame(
    Method = m,
    Sample = rep(common_samples, length(common_types)),
    CellType = rep(common_types, each = length(common_samples)),
    Estimated = as.vector(est[common_samples, common_types]),
    GroundTruth = as.vector(ground_truth[common_samples, common_types])
  )
}))

p2 <- ggplot(scatter_data, aes(x = GroundTruth, y = Estimated, color = CellType)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  facet_wrap(~Method, ncol = 2) + theme_minimal() +
  theme(legend.position = "bottom") +
  labs(title = "Estimated vs Ground Truth", x = "Ground Truth", y = "Estimated")

ggsave(file.path(VALIDATION_DIR, "scatter_all_methods.png"), p2, width = 10, height = 10, dpi = 300)

# Figure 3: Heatmap
celltype_perf_matrix <- celltype_cors_df %>%
  pivot_wider(names_from = CellType, values_from = Correlation) %>%
  tibble::column_to_rownames("Method") %>% as.matrix()

png(file.path(VALIDATION_DIR, "celltype_performance_heatmap.png"), width = 10, height = 6, units = "in", res = 300)
pheatmap(celltype_perf_matrix, display_numbers = TRUE, number_format = "%.2f",
         color = colorRampPalette(c("white", "steelblue"))(100),
         main = "Per Cell Type Correlation", cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

# Combined figure
combined <- (p1 / p2) + plot_annotation(title = "Pseudo-bulk Validation Summary")
ggsave(file.path(VALIDATION_DIR, "validation_summary.png"), combined, width = 12, height = 14, dpi = 300)

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------
message("\n=== Saving results ===")

saveRDS(ground_truth, file.path(VALIDATION_DIR, "ground_truth_proportions.rds"))
write.csv(ground_truth, file.path(VALIDATION_DIR, "ground_truth_proportions.csv"))
saveRDS(results_list, file.path(VALIDATION_DIR, "estimated_proportions_all.rds"))
saveRDS(evaluation_results, file.path(VALIDATION_DIR, "evaluation_results.rds"))
write.csv(summary_table, file.path(VALIDATION_DIR, "method_performance_summary.csv"), row.names = FALSE)

# Save pseudo-bulk for CIBERSORTx
pseudobulk_export <- data.frame(Gene = rownames(pseudobulk_counts), pseudobulk_counts)
write.table(pseudobulk_export, file.path(VALIDATION_DIR, "pseudobulk_for_cibersortx.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

message("\n=== VALIDATION COMPLETE ===")
message("\nMethod Ranking:")
for (i in 1:nrow(summary_table)) {
  message(sprintf("  %d. %s: r = %.4f, RMSE = %.4f", i, summary_table$Method[i],
                  summary_table$Correlation[i], summary_table$RMSE[i]))
}
message("\nFiles saved to output/validation/")
