# ============================================================
# 010_TNM stage–stratified expression analysis (ACRG / GSE62254)
# - Input: TMNstage_data.xlsx
#   * One column named "Stage" (I/II/III/IV)
#   * Remaining columns are proteins (expression values)
# - Output: One PDF per protein with jitter + boxplot and
#           adjacent-stage Wilcoxon rank-sum comparisons
# ============================================================

setwd("******************************************************4")

# -----------------------------
# 0) Packages
# -----------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggsignif)
  library(scales)   # axis tick formatting
})

# -----------------------------
# 1) Read data
# -----------------------------
input_file <- "TMNstage_data.xlsx"
df <- read_excel(input_file)

# Basic checks
stopifnot("Stage" %in% colnames(df))
cat("Stage distribution:\n")
print(table(df$Stage, useNA = "ifany"))

# -----------------------------
# 2) Reshape data: wide -> long
# -----------------------------
df_long <- df %>%
  pivot_longer(
    cols = -Stage,
    names_to = "Protein",
    values_to = "Expression"
  ) %>%
  mutate(
    Expression = as.numeric(Expression),
    Stage = factor(Stage, levels = c("I", "II", "III", "IV"))
  )

# -----------------------------
# 3) Color palette (ordered by stage progression)
# -----------------------------
stage_colors <- c(
  "I"   = "#91C7E6",
  "II"  = "#60ACEA",
  "III" = "#C86565",
  "IV"  = "#A13E45"
)

# -----------------------------
# 4) Define comparisons (adjacent stages)
# -----------------------------
comparisons <- list(
  c("I", "II"),
  c("II", "III"),
  c("III", "IV")
)

proteins <- unique(df_long$Protein)

# -----------------------------
# 5) Loop over proteins and generate plots
# -----------------------------
for (protein in proteins) {

  current_data <- df_long %>%
    filter(Protein == protein) %>%
    filter(!is.na(Expression))

  # Keep only comparisons where both stages exist in the data
  valid_comparisons <- comparisons[sapply(comparisons, function(pair) {
    all(pair %in% as.character(current_data$Stage))
  })]

  if (length(valid_comparisons) == 0) {
    message("Warning: protein '", protein, "' has insufficient stage data for comparisons. Skipped.")
    next
  }

  # Wilcoxon rank-sum test for each valid comparison
  p_values <- sapply(valid_comparisons, function(pair) {
    tmp <- current_data %>% filter(Stage %in% pair)
    wilcox.test(Expression ~ Stage, data = tmp)$p.value
  })

  # Format p-values for annotation
  p_annotations <- sapply(p_values, function(p) {
    if (p < 0.001) {
      sprintf("P = %.1e", p)
    } else if (p < 0.01) {
      sprintf("P = %.3f", p)
    } else {
      sprintf("P = %.2f", p)
    }
  })

  # Compute y-axis limits and annotation positions
  y_min_val <- min(current_data$Expression, na.rm = TRUE)
  y_max_val <- max(current_data$Expression, na.rm = TRUE)
  y_range <- y_max_val - y_min_val
  if (y_range == 0) y_range <- 0.1

  y_min <- y_min_val - 0.1 * y_range
  y_max <- y_max_val + 0.4 * y_range

  max_expr <- y_max_val
  offset <- 0.1 * y_range
  y_positions <- max_expr + offset * seq_len(length(valid_comparisons))

  # Plot: jitter + boxplot + significance bars
  p <- ggplot(current_data, aes(x = Stage, y = Expression, color = Stage, fill = Stage)) +
    geom_jitter(width = 0.1, size = 6, alpha = 0.8) +
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      color = "black",
      fill = NA,
      linewidth = 0.8
    ) +
    scale_color_manual(values = stage_colors) +
    scale_fill_manual(values = stage_colors) +
    scale_y_continuous(labels = label_number(accuracy = 0.1)) +
    geom_signif(
      comparisons = valid_comparisons,
      y_position = y_positions,
      tip_length = 0.02,
      size = 0.5,
      color = "black",
      textsize = 4,
      annotations = p_annotations
    ) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(
      x = "TNM Stage",
      y = paste0(protein, " Expression"),
      title = protein
    ) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 1),
      axis.ticks = element_line(color = "black", linewidth = 1),
      axis.ticks.length = unit(0.2, "cm"),
      legend.position = "none",
      axis.text = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      plot.title = element_text(hjust = 0.5, size = 18)
    )

  print(p)

  output_file <- paste0(protein, "_stage_expression_v2.pdf")
  ggsave(filename = output_file, plot = p, width = 5, height = 4, dpi = 300)
  cat("Saved:", output_file, "\n")
}
