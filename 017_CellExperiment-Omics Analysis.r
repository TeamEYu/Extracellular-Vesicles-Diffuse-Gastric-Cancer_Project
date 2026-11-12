# ------------------------------------------------------------
# 017_Differential Protein Analysis Pipeline
# Robust two-group comparison: High-salt vs Non-high-salt
# - Option A: Wilcoxon rank-sum test (nonparametric, robust)
# - Option B: Welch's t-test (parametric, unequal variance)
# Multiple testing correction: Benjamini–Hochberg (BH/FDR)
# Output:
#   - Wilcoxon results: 110_protein_salt_results_Wilcoxon.xlsx
#   - t-test results:   110_protein_salt_results_ttest.xlsx
#   - Volcano plots:    Volcano_Wilcoxon.pdf, Volcano_ttest.pdf
# ------------------------------------------------------------

rm(list = ls())

# -----------------------------
# 0) Packages
# -----------------------------
pkgs <- c("readxl", "tidyverse", "broom", "ggrepel", "writexl")
to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(readxl)
library(tidyverse)
library(broom)
library(ggrepel)
library(writexl)

# -----------------------------
# 1) Import data (wide format)
# -----------------------------
my_data <- read_excel("DGCHC_saltgrams.xlsx", sheet = 2)

# Define metadata columns (non-protein columns)
metadata_cols <- c("id", "salt_grams", "cancer_status", "age")
protein_cols <- setdiff(colnames(my_data), metadata_cols)

message("Detected ", length(protein_cols), " protein columns.")

# -----------------------------
# 2) Reshape to long format and define salt group
# -----------------------------
cutoff_value <- 5

data_long <- my_data %>%
  mutate(
    salt_group = ifelse(salt_grams > cutoff_value, "High_salt", "Non_high_salt"),
    salt_group = factor(salt_group, levels = c("Non_high_salt", "High_salt"))
  ) %>%
  pivot_longer(
    cols = all_of(protein_cols),
    names_to = "protein_name",
    values_to = "expression"
  ) %>%
  mutate(expression = as.numeric(expression))  # ensure numeric

# Optional QC: check group sizes
message("Group sizes:")
print(table(data_long$salt_group))

# -----------------------------
# 3) Descriptive statistics (mean + median per group)
# -----------------------------
descriptive_stats <- data_long %>%
  group_by(protein_name, salt_group) %>%
  summarise(
    mean_expression = mean(expression, na.rm = TRUE),
    median_expression = median(expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = salt_group,
    values_from = c(mean_expression, median_expression)
  )

# Helper to compute log2 fold-change safely
safe_log2_fc <- function(mean_high, mean_low, pseudocount = 1e-8) {
  log2((mean_high + pseudocount) / (mean_low + pseudocount))
}

# ============================================================
# Option A: Wilcoxon rank-sum test (robust)
# ============================================================
message("Running Wilcoxon rank-sum tests for all proteins...")

results_wilcox <- data_long %>%
  group_by(protein_name) %>%
  summarise(
    test_summary = list(tidy(wilcox.test(expression ~ salt_group))),
    .groups = "drop"
  ) %>%
  unnest(test_summary) %>%
  select(protein_name, p.value, statistic, method)

final_wilcox <- results_wilcox %>%
  left_join(descriptive_stats, by = "protein_name") %>%
  mutate(
    p_adj_BH = p.adjust(p.value, method = "BH"),
    mean_high_salt = mean_expression_High_salt,
    mean_non_high_salt = mean_expression_Non_high_salt,
    median_high_salt = median_expression_High_salt,
    median_non_high_salt = median_expression_Non_high_salt,
    log2_FoldChange = safe_log2_fc(mean_high_salt, mean_non_high_salt)
  ) %>%
  select(
    protein_name,
    p.value, p_adj_BH,
    mean_high_salt, mean_non_high_salt,
    median_high_salt, median_non_high_salt,
    log2_FoldChange,
    statistic, method
  ) %>%
  arrange(p_adj_BH)

message("Top 20 proteins (Wilcoxon, sorted by BH-adjusted P):")
print(final_wilcox, n = 20)

write_xlsx(final_wilcox, "110_protein_salt_results_Wilcoxon.xlsx")
message("Wilcoxon results saved to: 110_protein_salt_results_Wilcoxon.xlsx")

# -----------------------------
# Volcano plot (Wilcoxon)
# -----------------------------
plot_data_wilcox <- final_wilcox %>%
  mutate(
    significance = case_when(
      p_adj_BH < 0.05 & log2_FoldChange > 0.5  ~ "Significantly_up",
      p_adj_BH < 0.05 & log2_FoldChange < -0.5 ~ "Significantly_down",
      TRUE ~ "Not_significant"
    ),
    neg_log10_padj = -log10(p_adj_BH),
    protein_label = ifelse(p_adj_BH < 0.01, protein_name, NA_character_)
  )

volcano_wilcox <- ggplot(
  plot_data_wilcox,
  aes(x = log2_FoldChange, y = neg_log10_padj, color = significance)
) +
  geom_point(alpha = 0.7) +
  geom_text_repel(aes(label = protein_label), max.overlaps = 20) +
  scale_color_manual(values = c(
    "Significantly_up" = "red",
    "Significantly_down" = "blue",
    "Not_significant" = "grey"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  labs(
    title = "Volcano Plot: High-salt vs Non-high-salt (Wilcoxon)",
    x = "Log2(Fold Change)",
    y = "-Log10(BH-adjusted P)",
    color = "Significance"
  ) +
  theme_minimal()

print(volcano_wilcox)
ggsave("Volcano_Wilcoxon.pdf", volcano_wilcox, width = 7, height = 5, dpi = 300)

# ============================================================
# Option B: Welch's t-test (unequal variance)
# ============================================================
message("Running Welch's t-tests for all proteins...")

results_ttest <- data_long %>%
  group_by(protein_name) %>%
  summarise(
    test_summary = list(tidy(t.test(expression ~ salt_group, var.equal = FALSE))),
    .groups = "drop"
  ) %>%
  unnest(test_summary) %>%
  select(protein_name, p.value, statistic, method, conf.low, conf.high)

final_ttest <- results_ttest %>%
  left_join(descriptive_stats, by = "protein_name") %>%
  mutate(
    p_adj_BH = p.adjust(p.value, method = "BH"),
    mean_high_salt = mean_expression_High_salt,
    mean_non_high_salt = mean_expression_Non_high_salt,
    log2_FoldChange = safe_log2_fc(mean_high_salt, mean_non_high_salt)
  ) %>%
  select(
    protein_name,
    p.value, p_adj_BH,
    mean_high_salt, mean_non_high_salt,
    log2_FoldChange,
    statistic, conf.low, conf.high, method
  ) %>%
  arrange(p_adj_BH)

message("Top 20 proteins (t-test, sorted by BH-adjusted P):")
print(final_ttest, n = 20)

write_xlsx(final_ttest, "110_protein_salt_results_ttest.xlsx")
message("t-test results saved to: 110_protein_salt_results_ttest.xlsx")

# -----------------------------
# Volcano plot (t-test)
# -----------------------------
plot_data_ttest <- final_ttest %>%
  mutate(
    significance = case_when(
      p_adj_BH < 0.05 & log2_FoldChange > 0.5  ~ "Significantly_up",
      p_adj_BH < 0.05 & log2_FoldChange < -0.5 ~ "Significantly_down",
      TRUE ~ "Not_significant"
    ),
    neg_log10_padj = -log10(p_adj_BH),
    protein_label = ifelse(p_adj_BH < 0.01, protein_name, NA_character_)
  )

volcano_ttest <- ggplot(
  plot_data_ttest,
  aes(x = log2_FoldChange, y = neg_log10_padj, color = significance)
) +
  geom_point(alpha = 0.7) +
  geom_text_repel(aes(label = protein_label), max.overlaps = 20) +
  scale_color_manual(values = c(
    "Significantly_up" = "red",
    "Significantly_down" = "blue",
    "Not_significant" = "grey"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  labs(
    title = "Volcano Plot: High-salt vs Non-high-salt (Welch's t-test)",
    x = "Log2(Fold Change)",
    y = "-Log10(BH-adjusted P)",
    color = "Significance"
  ) +
  theme_minimal()

print(volcano_ttest)
ggsave("Volcano_ttest.pdf", volcano_ttest, width = 7, height = 5, dpi = 300)
