# ============================================================
# Protein expression heatmap (row-wise z-score) with pheatmap
# ============================================================

suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
})

# -----------------------------
# User inputs
# -----------------------------
input_file <- "sheet6.csv"          # e.g., "TCGA_42_tissue.csv" or "工作簿6.csv"
output_pdf <- "protein_heatmap.pdf"  # output figure

# -----------------------------
# 1) Load data (rows = proteins, cols = samples)
# -----------------------------
expr <- read.csv(input_file, row.names = 1, check.names = FALSE)
expr <- as.matrix(expr)

# Ensure numeric matrix
expr <- apply(expr, 2, function(x) as.numeric(as.character(x)))
expr <- as.matrix(expr)
rownames(expr) <- rownames(read.csv(input_file, row.names = 1, check.names = FALSE))
colnames(expr) <- colnames(read.csv(input_file, row.names = 1, check.names = FALSE))

# Basic checks
if (anyNA(expr)) {
  message("Warning: NA values detected in the input matrix. Z-scoring may introduce additional NA for zero-variance rows.")
}

# Remove rows with all NA
keep_nonempty <- rowSums(!is.na(expr)) > 0
expr <- expr[keep_nonempty, , drop = FALSE]

# -----------------------------
# 2) Row-wise z-score scaling (protein-wise)
#    z = (x - mean) / sd
# -----------------------------
row_zscore <- function(m) {
  t(apply(m, 1, function(v) {
    if (all(is.na(v))) return(rep(NA_real_, length(v)))
    s <- sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) {
      # zero variance -> return zeros (neutral), avoids all-NA z-scores
      return(rep(0, length(v)))
    }
    (v - mean(v, na.rm = TRUE)) / s
  }))
}

expr_z <- row_zscore(expr)
rownames(expr_z) <- rownames(expr)
colnames(expr_z) <- colnames(expr)

# -----------------------------
# 3) Sample annotations (top bar)
# -----------------------------
sample_names <- colnames(expr_z)

sample_groups <- factor(
  ifelse(grepl("^DNT", sample_names), "DNT",
         ifelse(grepl("^PT", sample_names) | grepl("PT", sample_names), "PT", "GCT")),
  levels = c("DNT", "PT", "GCT")
)

annotation_col <- data.frame(Group = sample_groups)
rownames(annotation_col) <- sample_names

annotation_colors <- list(
  Group = c(DNT = "#8FAADC", PT = "#F4B183", GCT = "#FFD966")
)

# -----------------------------
# 4) Heatmap color palette
# -----------------------------
heatmap_colors <- colorRampPalette(c("#003366", "white", "#B2182B"))(100)

# -----------------------------
# 5) Plot and save (single call)
# -----------------------------
pheatmap(
  mat = expr_z,
  color = heatmap_colors,
  scale = "none",
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  treeheight_row = 0,
  treeheight_col = 0,
  border_color = NA,
  fontsize_row = 10,
  main = "Protein Expression Heatmap (Row-wise Z-score)",
  filename = output_pdf,
  width = 8,
  height = 10
)

message("Done. Heatmap saved to: ", output_pdf)