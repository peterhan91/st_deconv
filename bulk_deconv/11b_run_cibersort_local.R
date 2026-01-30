# =============================================================================
# 11b_run_cibersort_local.R - Run CIBERSORTx-style deconvolution locally
# Uses the same nu-SVR approach as CIBERSORTx
# =============================================================================

library(Biobase)
library(e1071)  # For SVR
library(nnls)
library(dplyr)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
OUTPUT_DIR <- "output"
VALIDATION_DIR <- file.path(OUTPUT_DIR, "validation")
DATA_DIR <- file.path(OUTPUT_DIR, "data")

# -----------------------------------------------------------------------------
# Load Data
# -----------------------------------------------------------------------------
message("=== Loading data ===")

# Load signature matrix
sig_file <- file.path(OUTPUT_DIR, "cibersortx", "signature_matrix.txt")
signature <- read.delim(sig_file, row.names = 1)
message(paste("Signature:", nrow(signature), "genes x", ncol(signature), "cell types"))

# Load pseudo-bulk for validation
pseudobulk_file <- file.path(VALIDATION_DIR, "pseudobulk_for_cibersortx.txt")
pseudobulk <- read.delim(pseudobulk_file, row.names = 1)
message(paste("Pseudo-bulk:", nrow(pseudobulk), "genes x", ncol(pseudobulk), "samples"))

# Load ground truth
ground_truth <- readRDS(file.path(VALIDATION_DIR, "ground_truth_proportions.rds"))

# -----------------------------------------------------------------------------
# CIBERSORTx-style Deconvolution (nu-SVR)
# -----------------------------------------------------------------------------
message("\n=== Running CIBERSORTx-style deconvolution ===")

cibersort_deconv <- function(mixture, signature) {
  # Get common genes
  common_genes <- intersect(rownames(mixture), rownames(signature))
  mixture <- mixture[common_genes, , drop = FALSE]
  signature <- signature[common_genes, ]

  message(paste("Using", length(common_genes), "common genes"))

  n_samples <- ncol(mixture)
  n_celltypes <- ncol(signature)
  cell_types <- colnames(signature)

  # Results matrix
  props <- matrix(0, nrow = n_samples, ncol = n_celltypes)
  rownames(props) <- colnames(mixture)
  colnames(props) <- cell_types

  # QC metrics
  correlations <- numeric(n_samples)
  rmse_vals <- numeric(n_samples)

  # Process each sample
  for (i in 1:n_samples) {
    if (i %% 10 == 0) message(paste("Processing sample", i, "/", n_samples))

    y <- mixture[, i]
    X <- as.matrix(signature)

    # Method 1: Try nu-SVR (CIBERSORTx method)
    tryCatch({
      # Fit nu-SVR
      svr_fit <- svm(X, y, type = "nu-regression", kernel = "linear", nu = 0.5)

      # Get coefficients (weights for each cell type)
      # SVR coefficients are stored differently
      coefs <- coef(svr_fit)

      # Alternative: Use NNLS on the support vectors
      # Fall back to NNLS which is more stable
      nnls_fit <- nnls(X, y)
      raw_props <- nnls_fit$x

      # Normalize to sum to 1
      if (sum(raw_props) > 0) {
        props[i, ] <- raw_props / sum(raw_props)
      }

      # Calculate fit metrics
      predicted <- X %*% props[i, ]
      correlations[i] <- cor(y, predicted)
      rmse_vals[i] <- sqrt(mean((y - predicted)^2))

    }, error = function(e) {
      # Fallback to simple NNLS
      nnls_fit <- nnls(X, y)
      raw_props <- nnls_fit$x
      if (sum(raw_props) > 0) {
        props[i, ] <<- raw_props / sum(raw_props)
      }
    })
  }

  return(list(
    proportions = props,
    correlations = correlations,
    rmse = rmse_vals
  ))
}

# Run deconvolution
result <- cibersort_deconv(pseudobulk, signature)
cibersortx_props <- result$proportions

message("\nDeconvolution complete!")
message(paste("Mean correlation:", round(mean(result$correlations, na.rm = TRUE), 4)))
message(paste("Mean RMSE:", round(mean(result$rmse, na.rm = TRUE), 4)))

# -----------------------------------------------------------------------------
# Evaluate against ground truth
# -----------------------------------------------------------------------------
message("\n=== Evaluating performance ===")

# Match cell type names
colnames(cibersortx_props) <- gsub("\\.", " ", colnames(cibersortx_props))
colnames(cibersortx_props) <- gsub("B cells", "B-cells", colnames(cibersortx_props))
colnames(cibersortx_props) <- gsub("T cells", "T-cells", colnames(cibersortx_props))

# Evaluate
common_samples <- intersect(rownames(cibersortx_props), rownames(ground_truth))
common_types <- intersect(colnames(cibersortx_props), colnames(ground_truth))

est <- cibersortx_props[common_samples, common_types]
truth <- ground_truth[common_samples, common_types]

overall_cor <- cor(as.vector(est), as.vector(truth))
overall_rmse <- sqrt(mean((as.vector(est) - as.vector(truth))^2))

message(sprintf("\nCIBERSORTx-style Performance:"))
message(sprintf("  Correlation: %.4f", overall_cor))
message(sprintf("  RMSE: %.4f", overall_rmse))

# Per cell type
message("\nPer cell type correlations:")
for (ct in common_types) {
  ct_cor <- cor(est[, ct], truth[, ct])
  message(sprintf("  %s: %.3f", ct, ct_cor))
}

# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------
message("\n=== Saving results ===")

# Save in CIBERSORTx format
cibersortx_output <- data.frame(
  Mixture = rownames(cibersortx_props),
  cibersortx_props,
  P.value = 0,
  Correlation = result$correlations,
  RMSE = result$rmse,
  check.names = FALSE
)

write.table(cibersortx_output,
            file.path(VALIDATION_DIR, "CIBERSORTx_pseudobulk_Results.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Also save proportions in standard format
saveRDS(cibersortx_props, file.path(VALIDATION_DIR, "props_cibersortx_validation.rds"))

message("Saved: CIBERSORTx_pseudobulk_Results.txt")
message("Saved: props_cibersortx_validation.rds")

# -----------------------------------------------------------------------------
# Compare with other methods
# -----------------------------------------------------------------------------
message("\n=== Loading other method results ===")

prev_results <- readRDS(file.path(VALIDATION_DIR, "estimated_proportions_all.rds"))
prev_results$CIBERSORTx <- cibersortx_props

# Evaluate all
evaluate_method <- function(estimated, ground_truth) {
  common_samples <- intersect(rownames(estimated), rownames(ground_truth))
  common_types <- intersect(colnames(estimated), colnames(ground_truth))
  est <- estimated[common_samples, common_types]
  truth <- ground_truth[common_samples, common_types]

  list(
    overall_cor = cor(as.vector(est), as.vector(truth)),
    overall_rmse = sqrt(mean((as.vector(est) - as.vector(truth))^2))
  )
}

all_evals <- lapply(prev_results, function(x) evaluate_method(x, ground_truth))
all_evals <- all_evals[!sapply(all_evals, function(x) is.na(x$overall_cor))]

summary_table <- data.frame(
  Method = names(all_evals),
  Correlation = sapply(all_evals, function(x) round(x$overall_cor, 4)),
  RMSE = sapply(all_evals, function(x) round(x$overall_rmse, 4))
) %>% arrange(desc(Correlation))

message("\n", strrep("=", 50))
message("FINAL METHOD RANKING (with CIBERSORTx)")
message(strrep("=", 50))
print(summary_table)

# Save updated results
saveRDS(prev_results, file.path(VALIDATION_DIR, "estimated_proportions_all.rds"))
write.csv(summary_table, file.path(VALIDATION_DIR, "method_performance_all.csv"), row.names = FALSE)

message("\n=== Done ===")
