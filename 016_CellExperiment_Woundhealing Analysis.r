# ============================================================
# Scratch (Wound-Healing) Migration Assay (72 h)
# Plot: Hollow bars (mean ± SD) for NC vs NaCl across cell lines
# Output: PDF figure
# ============================================================

rm(list = ls())

# -----------------------------
# 0) Packages
# -----------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)
if (!requireNamespace("grid", quietly = TRUE)) install.packages("grid")
library(grid)

# -----------------------------
# 1) Input data (Table S1 - Scratch migration, 72 h)
# -----------------------------
plot_data <- data.frame(
  CellLine = factor(
    rep(c("GES-1", "MKN28", "MKN45"), each = 2),
    levels = c("GES-1", "MKN28", "MKN45")
  ),
  Group = factor(
    rep(c("NC", "NaCl"), 3),
    levels = c("NC", "NaCl")
  ),
  Mean = c(
    68.55, 91.31,  # GES-1
    52.31, 87.67,  # MKN28
    8.25,  39.04   # MKN45
  ),
  SD = c(
    2.22, 2.19,    # GES-1
    2.49, 0.97,    # MKN28
    4.25, 6.83     # MKN45
  )
)

# -----------------------------
# 2) Aesthetics
# -----------------------------
color_palette <- c("NC" = "#1c2d5a", "NaCl" = "#b2182b")

# Significance annotation positions (edit as needed)
sig_x_start <- c(0.65, 1.65, 2.65)
sig_x_end   <- c(1.35, 2.35, 3.35)
sig_y_line  <- c(100, 95, 50)
sig_y_text  <- c(103, 98, 53)
sig_label   <- "****"

# Bar geometry
dodge_width <- 0.7
bar_width   <- 0.6

# -----------------------------
# 3) Plot
# -----------------------------
p <- ggplot(plot_data, aes(x = CellLine, y = Mean, color = Group, group = Group)) +
  # Hollow bars
  geom_bar(
    stat = "identity",
    fill = "white",
    linewidth = 1,
    position = position_dodge(dodge_width),
    width = bar_width
  ) +
  # Error bars (mean ± SD)
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(dodge_width),
    width = 0.2,
    linewidth = 1.5
  ) +
  # Significance brackets + labels (manual)
  annotate(
    "segment",
    x = sig_x_start, xend = sig_x_end,
    y = sig_y_line,  yend = sig_y_line,
    linewidth = 1,
    color = "black"
  ) +
  annotate(
    "text",
    x = c(1, 2, 3),
    y = sig_y_text,
    label = sig_label,
    size = 6,
    color = "black"
  ) +
  # Scales and labels
  scale_color_manual(values = color_palette) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 110)) +
  labs(
    x = NULL,
    y = "Migration rate (%)",
    title = "Scratch (Wound-Healing) Migration Assay (72 h)"
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

    axis.text.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),

    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

print(p)

# -----------------------------
# 4) Save figure
# -----------------------------
ggsave(
  filename = "Fig2_ScratchMigration_72h.pdf",
  plot = p,
  width = 5,
  height = 4,
  device = "pdf"
)