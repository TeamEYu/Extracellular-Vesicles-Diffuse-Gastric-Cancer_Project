# ============================================================
# GO enrichment visualization (validation cohort)
# - Input: Excel sheet containing GO terms, counts, p.adjust, ontology
# - Output: Two-panel figure (lollipop + -log10(adj P) bars with term labels)
# ============================================================

# -----------------------------
# 0) Working directory & packages
# -----------------------------
setwd("D:/MAS_gastric_cancer_project/DGC_DEP_LW/v8/validation_samples")

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------
# 1) Load GO enrichment results
# -----------------------------
go_file  <- "GO_validate_all167.xlsx"
go_sheet <- 2

go_all_df <- read_excel(go_file, sheet = go_sheet)

# Sanity check for required columns (edit here if your column names differ)
required_cols <- c("ONTOLOGY", "Description", "Count", "p.adjust")
missing_cols <- setdiff(required_cols, colnames(go_all_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# -------------------------------------------------------
# 2) Data preprocessing
#    - group by ONTOLOGY
#    - sort by adjusted p-value within each ontology
#    - compute -log10(p.adjust)
#    - set GO term factor levels for top-to-bottom plotting
# -------------------------------------------------------
go_all_df <- go_all_df %>%
  group_by(ONTOLOGY) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  mutate(
    logP  = -log10(p.adjust),
    GOItem = factor(Description, levels = rev(unique(Description)))
  ) %>%
  ungroup()

print(head(go_all_df))

# -----------------------------
# 3) Color palette (discrete)
# -----------------------------
term_colors <- rev(c(
  "#D65656", "#5FAB5F", "#DDB370", "#E1A99A", "#CFE6A1", "#72C3E3",
  "#D43B63", "#D796C3", "#9B3683", "#294B2E", "#B49267", "#E69F2F",
  "#66A4B8", "#5BC5BA", "#EB2B1E", "#805097", "#6388B4", "#FFAE34"
))

# If the number of GO terms exceeds the palette length, recycle colors safely
n_terms <- nlevels(go_all_df$GOItem)
if (n_terms > length(term_colors)) {
  term_colors <- rep(term_colors, length.out = n_terms)
}
names(term_colors) <- levels(go_all_df$GOItem)

# -----------------------------
# 4) Panel 1: lollipop chart (Gene Count)
# -----------------------------
p1 <- ggplot(go_all_df, aes(x = Count, y = GOItem)) +
  geom_segment(aes(yend = GOItem, x = 0, xend = Count),
               linewidth = 0.8, color = "grey80") +
  geom_point(aes(color = GOItem), size = 6) +
  scale_color_manual(values = term_colors) +
  scale_x_continuous(
    limits = c(0, max(go_all_df$Count, na.rm = TRUE) + 2),
    expand = c(0, 0)
  ) +
  labs(title = "Gene Count", x = NULL, y = NULL) +
  theme_bw(base_size = 18) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x  = element_text(color = "black", size = 18),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "none"
  )

# -----------------------------
# 5) Panel 2: bar chart (-log10 adjusted P value) + term labels
# -----------------------------
p2 <- ggplot(go_all_df, aes(y = GOItem)) +
  geom_bar(
    aes(x = logP, fill = GOItem),
    stat = "identity",
    width = 0.6,
    color = "transparent",
    alpha = 0.7
  ) +
  geom_text(
    aes(x = 0.3, label = Description),
    hjust = 0,
    size = 5.5 * (18 / 16)
  ) +
  labs(title = expression(-log[10]("adjusted P value")), x = NULL, y = NULL) +
  scale_fill_manual(values = term_colors) +
  theme_bw(base_size = 18) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x  = element_text(color = "black", size = 18),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank(),
    legend.position = "none"
  )

# -----------------------------
# 6) Combine panels and export
# -----------------------------
final_plot <- p1 + p2 + plot_layout(widths = c(1, 1.5))
print(final_plot)

out_pdf <- "GO_all167_lollipop_bar_combined_v2.pdf"
ggsave(out_pdf, final_plot, width = 14, height = 9, dpi = 300)