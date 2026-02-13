# ============================================================
# Transwell invasion assay (Table S2): Bar plot with SD and
# significance annotations (NC vs NaCl) across four cell lines
# Output: Transwell_4CellLines_Wider.pdf
# ============================================================

# -----------------------------
# 0) Packages
# -----------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

# -----------------------------
# 1) Data (Table S2 - Transwell invasion)
# -----------------------------
plot_data <- data.frame(
  CellLine = factor(
    rep(c("GES-1", "MKN28", "MKN45", "KATOIII"), each = 2),
    levels = c("GES-1", "MKN28", "MKN45", "KATOIII")
  ),
  Group = factor(rep(c("NC", "NaCl"), times = 4), levels = c("NC", "NaCl")),
  Mean = c(
    84.33, 217.67,  # GES-1
    39.00, 99.33,   # MKN28
    39.33, 89.67,   # MKN45
    57.00, 150.33   # KATOIII
  ),
  SD = c(
    17.95, 37.63,   # GES-1
    7.00, 11.68,    # MKN28
    3.51, 12.58,    # MKN45
    12.53, 16.29    # KATOIII
  )
)

# -----------------------------
# 2) Aesthetics
# -----------------------------
color_palette <- c("NC" = "#1c2d5a", "NaCl" = "#b2182b")

# Significance labels (manual; replace with your exact p-value symbols if needed)
sig_labels <- rep("****", 4)

# Positions for significance bars (manually tuned to match y-range)
sig_y <- c(270, 120, 115, 180)   # bar y-positions for each cell line
sig_y_text <- sig_y + 10         # text y-positions

# -----------------------------
# 3) Plot
# -----------------------------
p <- ggplot(plot_data, aes(x = CellLine, y = Mean, color = Group, group = Group)) +

  # Bars (outlined, white fill)
  geom_bar(
    stat = "identity",
    fill = "white",
    linewidth = 1,
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +

  # Error bars (Mean ± SD)
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.8),
    width = 0.2,
    linewidth = 1
  ) +

  # Significance bars (manual segments across each paired group)
  annotate(
    "segment",
    x = c(0.6, 1.6, 2.6, 3.6),
    xend = c(1.4, 2.4, 3.4, 4.4),
    y = sig_y,
    yend = sig_y,
    linewidth = 1,
    color = "black"
  ) +
  annotate(
    "text",
    x = c(1, 2, 3, 4),
    y = sig_y_text,
    label = sig_labels,
    size = 6,
    color = "black"
  ) +

  # Scales and labels
  scale_color_manual(values = color_palette) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 300)) +
  labs(
    x = NULL,
    y = "Number of invaded cells",
    title = "Transwell Invasion Assay"
  ) +

  # Theme
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),

    axis.line = element_line(linewidth = 1, color = "black"),
    axis.ticks = element_line(linewidth = 1, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),

    axis.text.x = element_text(size = 12, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

print(p)

# -----------------------------
# 4) Save
# -----------------------------
ggsave(
  filename = "Transwell_4CellLines_Wider.pdf",
  plot = p,
  width = 6,
  height = 5
)