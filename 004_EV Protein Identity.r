# ============================================================
# Proteomics QC Summary: Number of Detected Proteins per Sample
# 1) Plasma proteomics: HC / Pre / HP / Post
# 2) Tissue proteomics: DNT / HPT(PT) / GCT
# Outputs:
#   - Plasma boxplot:        Protein_Boxplot_plasma.pdf
#   - Tissue boxplot:        Protein_Boxplot_tissue.pdf
#   - Per-group trend plots: HC_Group_Plot.pdf, Pre_Group_Plot.pdf, HP_Group_Plot.pdf, Post_Group_Plot.pdf
# ============================================================

rm(list = ls())

# -----------------------------
# 0) Packages
# -----------------------------
pkgs <- c("readxl", "dplyr", "tidyr", "ggplot2", "gridExtra")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# -----------------------------
# 1) Plasma data: import + cleaning
# -----------------------------
setwd("D:/MAS胃癌项目/DGC_DEP_LW/v8/")

plasma_raw <- read_excel("Plasma_data_origin.xlsx", sheet = 1) %>% as.data.frame()
rownames(plasma_raw) <- plasma_raw[[1]]
plasma_raw <- plasma_raw[, -1, drop = FALSE]

# Convert literal "NaN" strings to NA (only for character columns)
plasma_raw[] <- lapply(plasma_raw, function(x) {
  if (is.character(x)) x[x == "NaN"] <- NA
  x
})

# Convert all columns to numeric (keep rownames)
plasma_mat <- as.data.frame(lapply(plasma_raw, function(x) as.numeric(as.character(x))))
rownames(plasma_mat) <- rownames(plasma_raw)

# Optional: inspect count of values < 10 before filtering
lt10_counts <- colSums(plasma_mat < 10, na.rm = TRUE)
message("Counts of values < 10 per sample (first 75 columns shown):")
print(head(data.frame(Sample = colnames(plasma_mat), Counts_LT10 = lt10_counts), 75))

# Apply intensity threshold: set values < 10 to NA
plasma_mat[plasma_mat < 10] <- NA

# -----------------------------
# 2) Helper: compute detected-protein counts per sample for a group
# -----------------------------
calc_detected_proteins <- function(protein_matrix, sample_cols, group_name) {
  if (length(sample_cols) == 0) stop(paste0("No samples matched for group: ", group_name))

  group_data <- protein_matrix[, sample_cols, drop = FALSE]
  proteins_per_sample <- colSums(!is.na(group_data))

  out <- list(
    group = group_name,
    sample_count = length(sample_cols),
    mean_proteins = mean(proteins_per_sample),
    sd_proteins = sd(proteins_per_sample),
    proteins_per_sample = proteins_per_sample
  )

  cat("\n", group_name, " group:\n", sep = "")
  cat("  n = ", out$sample_count, "\n", sep = "")
  cat("  Detected proteins per sample (mean ± SD): ",
      round(out$mean_proteins, 2), " ± ", round(out$sd_proteins, 2), "\n", sep = "")
  cat("  Per-sample counts: ", paste(out$proteins_per_sample, collapse = ", "), "\n", sep = "")

  return(out)
}

# -----------------------------
# 3) Plasma groups: column selection by prefix
# -----------------------------
hc_cols   <- grep("^HC",   colnames(plasma_mat))
post_cols <- grep("^Post", colnames(plasma_mat))
hp_cols   <- grep("^HP",   colnames(plasma_mat))
pre_cols  <- grep("^Pre",  colnames(plasma_mat))

hc_res   <- calc_detected_proteins(plasma_mat, hc_cols,   "HC")
post_res <- calc_detected_proteins(plasma_mat, post_cols, "Post")
hp_res   <- calc_detected_proteins(plasma_mat, hp_cols,   "HP")
pre_res  <- calc_detected_proteins(plasma_mat, pre_cols,  "Pre")

plasma_summary <- data.frame(
  Group = c("HC", "Post", "HP", "Pre"),
  SampleCount = c(hc_res$sample_count, post_res$sample_count, hp_res$sample_count, pre_res$sample_count),
  MeanProteins = c(hc_res$mean_proteins, post_res$mean_proteins, hp_res$mean_proteins, pre_res$mean_proteins),
  SD = c(hc_res$sd_proteins, post_res$sd_proteins, hp_res$sd_proteins, pre_res$sd_proteins)
)

cat("\n===== Summary: mean detected proteins per sample (plasma) =====\n")
print(plasma_summary)

# -----------------------------
# 4) Plasma boxplot + jitter
# -----------------------------
plasma_long <- data.frame(
  Group = factor(
    rep(c("HC", "Post", "HP", "Pre"),
        times = c(length(hc_res$proteins_per_sample),
                  length(post_res$proteins_per_sample),
                  length(hp_res$proteins_per_sample),
                  length(pre_res$proteins_per_sample))),
    levels = c("HC", "Post", "HP", "Pre")
  ),
  Proteins = c(hc_res$proteins_per_sample,
               post_res$proteins_per_sample,
               hp_res$proteins_per_sample,
               pre_res$proteins_per_sample)
)

plasma_colors <- c("HC" = "#355C4C", "Post" = "#794191", "HP" = "#E4862B", "Pre" = "#2E8CC0")

p_plasma_box <- ggplot(plasma_long, aes(x = Group, y = Proteins)) +
  geom_boxplot(
    color = "black", fill = "white", outlier.shape = NA, linewidth = 0.8
  ) +
  geom_jitter(
    aes(color = Group), width = 0.2, alpha = 0.8, size = 6
  ) +
  scale_color_manual(values = plasma_colors) +
  coord_cartesian(ylim = c(2000, 6000)) +
  labs(x = NULL, y = "Number of detected proteins") +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    legend.position = "none",
    plot.title = element_blank()
  )

print(p_plasma_box)

ggsave(
  filename = "Protein_Boxplot_plasma.pdf",
  plot = p_plasma_box,
  width = 6, height = 5, units = "in", dpi = 300
)

# -----------------------------
# 5) Tissue data: import + group QC
# -----------------------------
setwd("D:/MAS胃癌项目/公司给的结果/组织样本的蛋白组搜库结果/")

tissue_mat <- read.csv("tissuedata_origin.csv", row.names = 1, check.names = FALSE)
tissue_mat <- as.data.frame(tissue_mat)

dnt_cols <- grep("^DNT", colnames(tissue_mat))
gct_cols <- grep("^GCT", colnames(tissue_mat))
pt_cols  <- grep("^PT",  colnames(tissue_mat))  # HPT corresponds to PT prefix in filenames

dnt_res <- calc_detected_proteins(tissue_mat, dnt_cols, "DNT")
gct_res <- calc_detected_proteins(tissue_mat, gct_cols, "GCT")
hpt_res <- calc_detected_proteins(tissue_mat, pt_cols,  "HPT")

tissue_summary <- data.frame(
  Group = factor(c("DNT", "HPT", "GCT"), levels = c("DNT", "HPT", "GCT")),
  SampleCount = c(dnt_res$sample_count, hpt_res$sample_count, gct_res$sample_count),
  MeanProteins = c(dnt_res$mean_proteins, hpt_res$mean_proteins, gct_res$mean_proteins),
  SD = c(dnt_res$sd_proteins, hpt_res$sd_proteins, gct_res$sd_proteins)
)

cat("\n===== Summary: mean detected proteins per sample (tissue) =====\n")
print(tissue_summary)

# -----------------------------
# 6) Tissue boxplot + jitter (colored boxes with transparency)
# -----------------------------
tissue_long <- data.frame(
  Group = factor(
    rep(c("DNT", "HPT", "GCT"),
        times = c(length(dnt_res$proteins_per_sample),
                  length(hpt_res$proteins_per_sample),
                  length(gct_res$proteins_per_sample))),
    levels = c("DNT", "HPT", "GCT")
  ),
  Proteins = c(dnt_res$proteins_per_sample,
               hpt_res$proteins_per_sample,
               gct_res$proteins_per_sample)
)

tissue_colors <- c("DNT" = "#355C4C", "HPT" = "#794191", "GCT" = "#2E8CC0")

p_tissue_box <- ggplot(tissue_long, aes(x = Group, y = Proteins, fill = Group)) +
  geom_boxplot(
    color = "black", outlier.shape = NA, linewidth = 0.8, alpha = 0.3
  ) +
  geom_jitter(
    aes(color = Group), width = 0.2, alpha = 0.8, size = 6
  ) +
  scale_fill_manual(values = tissue_colors) +
  scale_color_manual(values = tissue_colors) +
  coord_cartesian(
    ylim = c(min(tissue_long$Proteins, na.rm = TRUE) - 200,
             max(tissue_long$Proteins, na.rm = TRUE) + 200)
  ) +
  labs(x = NULL, y = "Number of detected proteins") +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 1),
    axis.ticks = element_line(color = "black", linewidth = 1),
    axis.ticks.length = unit(0.2, "cm"),
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    legend.position = "none"
  )

print(p_tissue_box)

ggsave(
  filename = "Protein_Boxplot_tissue.pdf",
  plot = p_tissue_box,
  width = 5, height = 4, units = "in", dpi = 300
)

# -----------------------------
# 7) Plasma: per-group trend plots (sample index vs detected proteins)
# -----------------------------
make_group_df <- function(group_name, counts_vec, prefix) {
  data.frame(
    Group = group_name,
    Sample = paste0(prefix, seq_along(counts_vec)),
    Proteins = as.numeric(counts_vec),
    GroupOrder = seq_along(counts_vec)
  )
}

hc_df   <- make_group_df("HC",   hc_res$proteins_per_sample,   "HC")
post_df <- make_group_df("Post", post_res$proteins_per_sample, "Post")
hp_df   <- make_group_df("HP",   hp_res$proteins_per_sample,   "HP")
pre_df  <- make_group_df("Pre",  pre_res$proteins_per_sample,  "Pre")

create_group_plot <- function(group_data, group_name, y_limits = c(2000, 6000)) {
  ggplot(group_data, aes(x = GroupOrder, y = Proteins)) +
    geom_point(size = 2.5, color = "#3366CC", alpha = 0.8) +
    geom_smooth(
      method = "loess", se = TRUE,
      color = "#FF6666", fill = "#FFCCCC",
      span = 0.8, linewidth = 1
    ) +
    geom_hline(
      yintercept = mean(group_data$Proteins, na.rm = TRUE),
      linetype = "dashed", color = "gray50", linewidth = 0.8
    ) +
    coord_cartesian(ylim = y_limits) +
    labs(
      title = paste0(group_name, " Group (n = ", nrow(group_data), ")"),
      x = "Sample index",
      y = "Number of detected proteins"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)
    )
}

p_hc   <- create_group_plot(hc_df,   "HC")
p_post <- create_group_plot(post_df, "Post")
p_hp   <- create_group_plot(hp_df,   "HP")
p_pre  <- create_group_plot(pre_df,  "Pre")

ggsave("HC_Group_Plot.pdf",   p_hc,   width = 5, height = 4, units = "in", dpi = 300)
ggsave("Post_Group_Plot.pdf", p_post, width = 5, height = 4, units = "in", dpi = 300)
ggsave("HP_Group_Plot.pdf",   p_hp,   width = 5, height = 4, units = "in", dpi = 300)
ggsave("Pre_Group_Plot.pdf",  p_pre,  width = 5, height = 4, units = "in", dpi = 300)

# Optional: show a 2x2 panel in RStudio
grid.arrange(p_hc, p_hp, p_pre, p_post, nrow = 2, ncol = 2, padding = unit(1, "lines"))