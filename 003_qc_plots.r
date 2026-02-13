# ============================================================
# 02_qc_plots.R — QC visualizations
# ============================================================

plot_boxplots <- function(before_mat, after_mat, out_pdf) {
  close_all_devices()
  pdf(out_pdf, width = 12, height = 4)
  par(mfrow = c(1, 2))
  boxplot(before_mat, main = "Before normalization", ylab = "Abundance", outline = FALSE)
  boxplot(after_mat,  main = "After median normalization + log2", ylab = "log2(Abundance + 1)", outline = FALSE)
  dev.off()
}

plot_density_before_after <- function(before_mat, after_mat, out_pdf) {
  before_vals <- as.numeric(before_mat[!is.na(before_mat)])
  after_vals  <- as.numeric(after_mat[!is.na(after_mat)])

  df <- data.frame(
    value = c(before_vals, after_vals),
    stage = rep(c("Before", "After"), times = c(length(before_vals), length(after_vals)))
  )

  p <- ggplot(df, aes(x = value, color = stage)) +
    geom_density(linewidth = 1.0) +
    labs(title = "Distribution before vs after processing",
         x = "Value", y = "Density") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(out_pdf, p, width = 7, height = 5, dpi = 300)
}

plot_missingness_map <- function(original_mat, out_pdf, group_rules = NULL) {
  # 1 = missing in original, 0 = observed
  miss01 <- matrix(ifelse(is.na(original_mat), 1, 0),
                   nrow = nrow(original_mat), ncol = ncol(original_mat),
                   dimnames = dimnames(original_mat))

  annotation_col <- NULL
  if (!is.null(group_rules)) {
    grp <- rep("Other", ncol(original_mat))
    for (nm in names(group_rules)) {
      idx <- grep(group_rules[[nm]], colnames(original_mat))
      grp[idx] <- nm
    }
    annotation_col <- data.frame(Group = factor(grp))
    rownames(annotation_col) <- colnames(original_mat)
  }

  close_all_devices()
  pdf(out_pdf, width = 12, height = 8)
  pheatmap(
    mat = miss01,
    color = c("#DCEFFF", "#FF4C4C"),
    breaks = c(-0.1, 0.5, 1.1),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    show_rownames = FALSE,
    show_colnames = TRUE,
    fontsize_col = 6,
    main = "Missingness Map (1 = Missing in original)"
  )
  dev.off()
}

run_pca <- function(mat, group_vec, out_pdf, title = "PCA") {
  stopifnot(length(group_vec) == ncol(mat))
  pca <- PCA(t(mat), graph = FALSE)
  p <- fviz_pca_ind(
    pca,
    geom.ind = "point",
    col.ind = as.factor(group_vec),
    addEllipses = TRUE,
    ellipse.type = "confidence",
    ellipse.level = 0.95,
    repel = TRUE,
    title = title
  ) + theme(plot.title = element_text(hjust = 0.5))

  ggsave(out_pdf, p, width = 8, height = 6, device = "pdf")
}