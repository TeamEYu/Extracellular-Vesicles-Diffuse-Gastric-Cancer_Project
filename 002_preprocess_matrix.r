# ============================================================
# 002_preprocess_matrix.R — Preprocessing utilities and pipeline
# ============================================================

replace_string_NaN <- function(df) {
  df[] <- lapply(df, function(x) {
    if (is.character(x)) x[x == "NaN"] <- NA
    x
  })
  df
}

to_numeric_preserve_rownames <- function(df) {
  rn <- rownames(df)
  out <- as.data.frame(lapply(df, function(x) as.numeric(as.character(x))))
  rownames(out) <- rn
  out
}

read_protein_matrix <- function(path, sheet = 1) {
  df <- readxl::read_excel(path, sheet = sheet) |> as.data.frame()
  rownames(df) <- df[[1]]
  df <- df[, -1, drop = FALSE]
  df <- replace_string_NaN(df)
  df <- to_numeric_preserve_rownames(df)
  df
}

count_below_threshold <- function(mat, threshold = 10) {
  colSums(mat < threshold, na.rm = TRUE)
}

set_below_threshold_to_NA <- function(mat, threshold = 10) {
  mat[mat < threshold] <- NA
  mat
}

median_normalize <- function(mat) {
  sample_medians <- apply(mat, 2, median, na.rm = TRUE)
  target_median  <- median(sample_medians, na.rm = TRUE)
  norm_factors   <- target_median / sample_medians
  sweep(mat, 2, norm_factors, "*")
}

log2_transform <- function(mat, pseudocount = 1) {
  log2(mat + pseudocount)
}

missing_rate_by_row <- function(mat) {
  rowMeans(is.na(mat))
}

filter_by_missingness <- function(mat, max_missing = 0.5) {
  mr <- missing_rate_by_row(mat)
  keep <- mr <= max_missing
  list(
    data = mat[keep, , drop = FALSE],
    removed = names(mr)[!keep],
    missing_rate = mr
  )
}

preprocess_pipeline <- function(mat, threshold = 10, max_missing = 0.5, pseudocount = 1) {

  message("Summarizing values below threshold (", threshold, ") prior to masking...")
  lt_counts <- count_below_threshold(mat, threshold = threshold)

  message("Masking values < ", threshold, " as NA...")
  mat_masked <- set_below_threshold_to_NA(mat, threshold = threshold)

  message("Performing median normalization (column-wise)...")
  mat_norm <- median_normalize(mat_masked)

  message("Applying log2 transformation with pseudocount = ", pseudocount, " ...")
  mat_log2 <- log2_transform(mat_norm, pseudocount = pseudocount)

  message("Filtering proteins with missingness rate > ", max_missing, " ...")
  filt <- filter_by_missingness(mat_log2, max_missing = max_missing)

  list(
    raw = mat,
    lt_counts = lt_counts,
    masked = mat_masked,
    normalized = mat_norm,
    log2 = mat_log2,
    filtered = filt$data,
    removed_proteins = filt$removed,
    missing_rate = filt$missing_rate
  )
}
