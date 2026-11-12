# ============================================================
# 007_ACRG (GSE62254) diffuse-type only: candidate genes survival
# Pipeline:
#   1) Load expression matrix (gene x samples)
#   2) Keep diffuse-type samples only (remove columns containing "non-diffuse")
#   3) Subset candidate genes by gene symbol intersection
#   4) Collapse duplicated probes per gene (take max across probes)
#   5) Transpose expression to sample x gene and merge with clinical data
#   6) Univariate survival analysis per gene:
#        - Log-rank test (median split: High vs Low)
#        - Cox PH model (HR, 95% CI, Wald P)
#        - Optional multiple testing correction (BH FDR)
#   7) Export results and generate KM curves for significant genes
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(writexl)
  library(survival)
  library(broom)
  library(ggplot2)
  library(ggsurvfit)
})

# -----------------------------
# 0) Working directory & inputs
# -----------------------------
setwd("D:/MAS_gastric_cancer_project/DGC_DEP_LW/v8/public_datasets/ACRG_GSE62254")

expr_file      <- "GSE62254_series_matrix_expression.xlsx"
candidate_file <- "KM_candidate_gene_list.xlsx"
clinical_file  <- "ACRG_GSE62254_clinicaldata.xlsx"
clinical_sheet <- 7

# -----------------------------
# Outputs (all English names)
# -----------------------------
out_expr_diffuse_only   <- "diffuse_only_expression_matrix.xlsx"
out_expr_candidates     <- "candidate_gene_expression_matrix.xlsx"
out_expr_dedup_max      <- "candidate_gene_expression_dedup_max.xlsx"
out_merged              <- "candidate_expression_merged_with_clinical.xlsx"
out_survival_results    <- "candidate_genes_univariate_survival_results.xlsx"
out_survival_significant<- "candidate_genes_significant.xlsx"

km_dir <- "KM_plots"
if (!dir.exists(km_dir)) dir.create(km_dir)

# -----------------------------
# User options
# -----------------------------
alpha <- 0.05
use_fdr <- TRUE
fdr_method <- "BH"

km_colors <- c("High expression" = "#c9496b", "Low expression" = "#396CB1")

# ------------------------------------------------------------
# 1) Load expression matrix and keep diffuse-type samples only
# ------------------------------------------------------------
expr_raw <- read_excel(expr_file)

# Remove columns whose names contain "non-diffuse"
non_diffuse_idx <- grep("non-diffuse", colnames(expr_raw), ignore.case = TRUE)
cat("Number of columns containing 'non-diffuse': ", length(non_diffuse_idx), "\n", sep = "")

expr_diffuse_only <- expr_raw %>%
  select(-all_of(non_diffuse_idx))

write_xlsx(expr_diffuse_only, out_expr_diffuse_only)

# -----------------------------------------
# 2) Candidate gene intersection (by symbol)
# -----------------------------------------
expr_diffuse_only <- read_excel(out_expr_diffuse_only)
cand <- read_excel(candidate_file)

stopifnot("Gene Symbol" %in% colnames(expr_diffuse_only))
stopifnot("Gene" %in% colnames(cand))

candidate_genes <- unique(cand$Gene) %>% na.omit()

expr_candidates <- expr_diffuse_only %>%
  filter(`Gene Symbol` %in% candidate_genes)

cat("Rows in diffuse-only expression matrix: ", nrow(expr_diffuse_only), "\n", sep = "")
cat("Rows after intersecting candidate genes: ", nrow(expr_candidates), "\n", sep = "")

write_xlsx(expr_candidates, out_expr_candidates)

# ----------------------------------------------------
# 3) Collapse duplicated probes per gene (take the max)
#    Apply max only to numeric sample columns
# ----------------------------------------------------
expr_candidates <- read_excel(out_expr_candidates)

sample_cols <- setdiff(colnames(expr_candidates), "Gene Symbol")

expr_candidates <- expr_candidates %>%
  mutate(across(all_of(sample_cols), as.numeric))

expr_dedup_max <- expr_candidates %>%
  group_by(`Gene Symbol`) %>%
  summarise(
    across(all_of(sample_cols), ~ max(.x, na.rm = TRUE)),
    .groups = "drop"
  )

cat("Rows before deduplication: ", nrow(expr_candidates), "\n", sep = "")
cat("Rows after deduplication (max across probes): ", nrow(expr_dedup_max), "\n", sep = "")

write_xlsx(expr_dedup_max, out_expr_dedup_max)

# ------------------------------------------------------------
# 4) Transpose to sample x gene and merge with clinical dataset
#    Expression: gene x sample (columns are sample IDs)
#    Clinical: contains sample ID in "GEO_ID"
# ------------------------------------------------------------
expr_dedup_max <- read_excel(out_expr_dedup_max)
clinical <- read_excel(clinical_file, sheet = clinical_sheet)

stopifnot("GEO_ID" %in% colnames(clinical))
stopifnot("Death" %in% colnames(clinical))
stopifnot("OS.m" %in% colnames(clinical))

expr_mat <- expr_dedup_max %>% as.data.frame()
rownames(expr_mat) <- expr_mat$`Gene Symbol`
expr_mat$`Gene Symbol` <- NULL

expr_t <- as.data.frame(t(expr_mat)) %>%
  rownames_to_column(var = "GEO_ID")

merged <- inner_join(clinical, expr_t, by = "GEO_ID")

cat("Clinical rows: ", nrow(clinical), "\n", sep = "")
cat("Expression samples (diffuse-only): ", nrow(expr_t), "\n", sep = "")
cat("Rows after merge (common samples): ", nrow(merged), "\n", sep = "")

write_xlsx(merged, out_merged)

# ------------------------------------------------------------
# 5) Univariate survival analysis for each candidate gene
#    - Median split: High vs Low
#    - Log-rank P + Cox HR (High vs Low) + 95% CI + Cox P
#    - Optional FDR correction
# ------------------------------------------------------------
merged <- read_excel(out_merged)

clinical_cols <- c("GEO_ID", "Death", "OS.m")
gene_cols <- setdiff(colnames(merged), clinical_cols)

results <- tibble(
  Gene = character(),
  Median_Expression = numeric(),
  HR = numeric(),
  HR_low = numeric(),
  HR_high = numeric(),
  Log_rank_P = numeric(),
  Cox_P = numeric()
)

for (g in gene_cols) {
  tmp <- merged %>%
    select(all_of(c("OS.m", "Death", g))) %>%
    rename(Expression = all_of(g)) %>%
    mutate(Expression = as.numeric(Expression)) %>%
    filter(!is.na(OS.m), !is.na(Death), !is.na(Expression))

  if (nrow(tmp) < 10) next

  med <- median(tmp$Expression, na.rm = TRUE)

  tmp <- tmp %>%
    mutate(
      Group = factor(
        ifelse(Expression > med, "High", "Low"),
        levels = c("Low", "High")
      )
    )

  lr <- survdiff(Surv(OS.m, Death) ~ Group, data = tmp)
  lr_p <- 1 - pchisq(lr$chisq, df = 1)

  cox <- coxph(Surv(OS.m, Death) ~ Group, data = tmp)
  cox_tidy <- tidy(cox, conf.int = TRUE)

  hr      <- exp(cox_tidy$estimate)
  hr_low  <- exp(cox_tidy$conf.low)
  hr_high <- exp(cox_tidy$conf.high)
  cox_p   <- cox_tidy$p.value

  results <- results %>%
    add_row(
      Gene = g,
      Median_Expression = med,
      HR = hr,
      HR_low = hr_low,
      HR_high = hr_high,
      Log_rank_P = lr_p,
      Cox_P = cox_p
    )

  cat("Finished: ", g,
      " | HR=", round(hr, 3),
      " | 95%CI=[", round(hr_low, 3), ", ", round(hr_high, 3), "]",
      " | Cox P=", signif(cox_p, 3),
      " | Log-rank P=", signif(lr_p, 3), "\n",
      sep = "")
}

if (use_fdr) {
  results <- results %>%
    mutate(
      Cox_FDR = p.adjust(Cox_P, method = fdr_method),
      LogRank_FDR = p.adjust(Log_rank_P, method = fdr_method)
    )
}

results <- results %>%
  arrange(Cox_P) %>%
  mutate(
    Cox_Significance = case_when(
      Cox_P < 0.001 ~ "***",
      Cox_P < 0.01 ~ "**",
      Cox_P < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    LogRank_Significance = case_when(
      Log_rank_P < 0.001 ~ "***",
      Log_rank_P < 0.01 ~ "**",
      Log_rank_P < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

print(results, n = nrow(results))
write_xlsx(results, out_survival_results)

if (use_fdr) {
  significant <- results %>% filter(Cox_FDR < alpha)
  cat("Significant genes (Cox FDR < ", alpha, "): ", nrow(significant), "\n", sep = "")
} else {
  significant <- results %>% filter(Cox_P < alpha)
  cat("Significant genes (Cox P < ", alpha, "): ", nrow(significant), "\n", sep = "")
}

write_xlsx(significant, out_survival_significant)

# ------------------------------------------------------------
# 6) Kaplan–Meier curves for significant genes
# ------------------------------------------------------------
try(dev.off(), silent = TRUE)

for (g in significant$Gene) {

  plot_df <- merged %>%
    select(all_of(c("OS.m", "Death", g))) %>%
    rename(Expression = all_of(g)) %>%
    mutate(Expression = as.numeric(Expression)) %>%
    filter(!is.na(OS.m), !is.na(Death), !is.na(Expression)) %>%
    mutate(
      Group = ifelse(Expression > median(Expression, na.rm = TRUE),
                     "High expression", "Low expression"),
      Group = factor(Group, levels = c("High expression", "Low expression"))
    )

  lr <- survdiff(Surv(OS.m, Death) ~ Group, data = plot_df)
  pval <- 1 - pchisq(lr$chisq, df = 1)
  p_txt <- ifelse(pval < 0.001, format(pval, scientific = TRUE, digits = 2), sprintf("%.3f", pval))

  fit <- survfit(Surv(OS.m, Death) ~ Group, data = plot_df)

  p <- fit %>%
    ggsurvfit(linewidth = 1) +
    annotate("text", x = 2.5, y = 0.25, label = "Log-rank", size = 6, hjust = 0) +
    annotate("text", x = 2.5, y = 0.15, label = paste0("p = ", p_txt), size = 6, hjust = 0) +
    add_risktable(
      risktable_height = 0.25,
      risktable_stats = c("{n.risk}"),
      stats_label = list(n.risk = "Number at risk"),
      size = 4
    ) +
    add_censor_mark(size = 3, shape = 73) +
    add_quantile(y_value = 0.5, linetype = "dashed", linewidth = 0.8) +
    labs(
      x = "Survival time (months)",
      y = "Overall survival (%)",
      color = g
    ) +
    scale_x_continuous(breaks = seq(0, 120, 30), limits = c(0, 121), expand = c(0.02, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25), expand = c(0, 0)) +
    scale_color_manual(values = km_colors) +
    scale_fill_manual(values = km_colors) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 12, color = "black"),
      legend.position = "top",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )

  ggsave(
    filename = file.path(km_dir, paste0(g, "_KM_curve.pdf")),
    plot = p,
    width = 8,
    height = 7,
    dpi = 300
  )
}
