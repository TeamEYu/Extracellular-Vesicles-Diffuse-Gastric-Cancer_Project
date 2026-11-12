# ============================================================
# 005_run_plasma_or_tissue.R — End-to-end run (example template)
# ============================================================

source("001_setup.R")
source("002_preprocess_matrix.R")
source("003_qc_plots.R")

# -----------------------------
# A) Load matrix
# -----------------------------
# Example (change to your file):
mat <- read_protein_matrix("tissuesample_nofiltered.xlsx", sheet = 1)

# -----------------------------
# B) Preprocess
# -----------------------------
res <- preprocess_pipeline(mat, threshold = 10, max_missing = 0.5, pseudocount = 1)

# Export <10 summary
lt_df <- data.frame(Sample = names(res$lt_counts), Counts_LT10 = as.integer(res$lt_counts))
write.csv(lt_df, "per_sample_counts_below10.csv", row.names = FALSE)

# Export filtered log2 matrix
write.csv(res$filtered, "log2_medianNormalized_filtered50.csv", row.names = TRUE, quote = TRUE)

# -----------------------------
# C) QC plots
# -----------------------------
plot_boxplots(res$masked, res$filtered, "QC_boxplot_before_after.pdf")
plot_density_before_after(res$masked, res$filtered, "QC_density_before_after.pdf")

# Missingness map with group annotations (customize rules)
# Example for tissue:
group_rules <- list(
  DNT = "^DNT",
  GCT = "^GCT",
  PT  = "^PT"
)
plot_missingness_map(res$masked, "QC_missingness_map.pdf", group_rules = group_rules)

# -----------------------------
# D) PCA (customize group labels; MUST match ncol(res$filtered))
# -----------------------------
# Example tissue counts: DNT=20, GCT=20, PT=18
group_vec <- c(rep("DNT", 20), rep("GCT", 20), rep("PT", 18))
run_pca(res$filtered, group_vec, "QC_PCA.pdf", title = "PCA after preprocessing")
