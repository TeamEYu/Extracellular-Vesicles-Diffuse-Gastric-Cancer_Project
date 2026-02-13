# ============================================================
# Differential expression analysis (DGC vs HC) using limma
# - Input: TCGA_validate_data.xlsx (sheet 1), first column = feature/protein IDs
# - Steps:
#   1) Read and clean data; convert all expression columns to numeric
#   2) Filter features with excessive zeros
#   3) Build design matrix with HC as reference
#   4) Fit limma model + empirical Bayes moderation
#   5) Export DE results and regulation labels
# ============================================================

# -----------------------------
# 0) Packages
# -----------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(limma)
})

# -----------------------------
# 1) Load expression matrix
# -----------------------------
infile <- "TCGA_validate_data.xlsx"
sheet_id <- 1

data_raw <- read_excel(infile, sheet = sheet_id) |> as.data.frame()

# Assume the first column contains feature/protein identifiers
rownames(data_raw) <- data_raw[[1]]
expr_df <- data_raw[, -1, drop = FALSE]

str(expr_df)

# Preserve original row names (feature IDs)
original_rownames <- rownames(expr_df)

# Convert all columns to numeric while preserving row names
expr_mat <- as.data.frame(lapply(expr_df, function(x) as.numeric(as.character(x))))
rownames(expr_mat) <- original_rownames

# Confirm numeric conversion
print(sapply(expr_mat, class))

# -----------------------------
# 2) Filter features by zero inflation
#    - Keep features with number of zeros <= zero_threshold
#    - Note: adjust zero_threshold based on sample size and sparsity
# -----------------------------
zero_threshold <- 52
zero_counts <- rowSums(expr_mat == 0, na.rm = TRUE)

expr_mat_filt <- expr_mat[zero_counts <= zero_threshold, , drop = FALSE]
cat("Features before filtering:", nrow(expr_mat), "\n")
cat("Features after filtering :", nrow(expr_mat_filt), "\n")

# -----------------------------
# 3) Define group labels and subset (DGC vs HC)
#    IMPORTANT: Ensure the order matches the columns of expr_mat_filt
# -----------------------------
group <- c(rep("DGC", 68),
           rep("HC",  36))

if (length(group) != ncol(expr_mat_filt)) {
  stop("Group label length (", length(group), 
       ") does not match number of samples (", ncol(expr_mat_filt), ").")
}

# Select only DGC and HC samples (kept for extensibility if more groups exist)
selected_cols <- group %in% c("DGC", "HC")
cat("Selected samples (DGC/HC):", sum(selected_cols), "\n")

expr_matrix <- expr_mat_filt[, selected_cols, drop = FALSE]
group_sel <- factor(group[selected_cols], levels = c("HC", "DGC"))  # HC as reference

# -----------------------------
# 4) Design matrix (reference-level coding)
#    Model: Expression ~ Group (DGC vs HC)
# -----------------------------
design <- model.matrix(~ group_sel)
colnames(design) <- c("Intercept", "DGCvsHC")

# -----------------------------
# 5) limma differential expression
# -----------------------------
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)

# Extract results for DGC vs HC
results <- topTable(
  fit,
  coef = "DGCvsHC",
  number = Inf,
  adjust.method = "fdr"
)

# -----------------------------
# 6) Add regulation labels
#    logFC threshold: |log2FC| > 0.585 (~1.5-fold change)
#    significance: FDR < 0.05
# -----------------------------
logfc_cutoff <- 0.585
fdr_cutoff <- 0.05

results$expression_level <- "NoSig"
results$expression_level[results$logFC >  logfc_cutoff & results$adj.P.Val < fdr_cutoff] <- "Upregulated"
results$expression_level[results$logFC < -logfc_cutoff & results$adj.P.Val < fdr_cutoff] <- "Downregulated"

cat("Differential expression summary:\n")
print(table(results$expression_level))

# -----------------------------
# 7) Save outputs
# -----------------------------
out_csv <- "DGC_vs_HC_TCGA_results.csv"
write.csv(results, file = out_csv, row.names = TRUE)

cat("Saved results to:", out_csv, "\n")