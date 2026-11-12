# ============================================================
# 012_Plasma PTM differential analysis (HC vs Pre)
# - Data cleaning: replace "Filtered" with NA, numeric casting,
#   and censor values < 10 as NA
# - Group-wise imputation: fill NA with half of the group-specific minimum
# - Differential testing: Wilcoxon rank-sum test + FDR correction
# - Outputs: a tidy result table saved as CSV
# ============================================================

# -----------------------------
# 0) Environment setup
# -----------------------------
rm(list = ls())

# Set working directory (edit if needed)
setwd("*****************************************************")

# -----------------------------
# 1) Load required packages
# -----------------------------
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")

library(readxl)
library(dplyr)
library(purrr)

# -----------------------------
# 2) Read data
# -----------------------------
ptm_data <- read_excel("Plasma_Curated_Modifications.xlsx", sheet = 7)

str(ptm_data)
colnames(ptm_data)

# -----------------------------
# 3) Data preprocessing
#    Assumption: columns 4:43 are PTM intensity/abundance columns.
#    - Convert "Filtered" to NA
#    - Cast to numeric
#    - Set values < 10 to NA (censoring threshold)
# -----------------------------
ptm_data_clean <- ptm_data %>%
  mutate(across(4:43, ~ ifelse(.x == "Filtered", NA, .x))) %>%
  mutate(across(4:43, ~ as.numeric(.x))) %>%
  mutate(across(4:43, ~ ifelse(!is.na(.x) & .x < 10, NA, .x)))

# -----------------------------
# 4) Define group columns (HC vs Pre)
#    Option A: infer by column name prefix
# -----------------------------
hc_cols_auto  <- grep("^HC",  colnames(ptm_data_clean), value = TRUE)
pre_cols_auto <- grep("^Pre", colnames(ptm_data_clean), value = TRUE)

# If you prefer a strict manual list (recommended for reproducibility),
# uncomment and use the following.
hc_cols <- c(
  "HC1","HC2","HC3","HC4","HC5","HC6","HC7","HC8","HC9","HC10",
  "HC11","HC12","HC13","HC14","HC15","HC16","HC17","HC18","HC19","HC20"
)
pre_cols <- c(
  "Pre1","Pre2","Pre3","Pre4","Pre5","Pre6","Pre7","Pre8","Pre9","Pre10",
  "Pre11","Pre12","Pre13","Pre14","Pre15","Pre16","Pre17","Pre18","Pre19","Pre20"
)

# Safety checks
stopifnot(all(hc_cols %in% colnames(ptm_data_clean)))
stopifnot(all(pre_cols %in% colnames(ptm_data_clean)))

# -----------------------------
# 5) Group-wise minimum-based imputation
#    Fill missing values within each group using (group minimum)/2
# -----------------------------
hc_min <- ptm_data_clean %>%
  select(all_of(hc_cols)) %>%
  unlist(use.names = FALSE) %>%
  (\(x) x[!is.na(x)])() %>%
  min()

pre_min <- ptm_data_clean %>%
  select(all_of(pre_cols)) %>%
  unlist(use.names = FALSE) %>%
  (\(x) x[!is.na(x)])() %>%
  min()

hc_fill_value  <- hc_min  / 2
pre_fill_value <- pre_min / 2

ptm_data_filled <- ptm_data_clean %>%
  mutate(across(all_of(hc_cols),  ~ ifelse(is.na(.x), hc_fill_value,  .x))) %>%
  mutate(across(all_of(pre_cols), ~ ifelse(is.na(.x), pre_fill_value, .x)))

cat("HC group minimum: ", hc_min,  " | Imputation value (min/2): ", hc_fill_value,  "\n", sep = "")
cat("Pre group minimum:", pre_min, " | Imputation value (min/2): ", pre_fill_value, "\n", sep = "")

# Quick inspection (optional)
print(head(ptm_data_filled %>% select(Gene, all_of(hc_cols[1:2]), all_of(pre_cols[1:2]))))

# -----------------------------
# 6) Differential analysis (Wilcoxon rank-sum test)
#    - Maintain original row order
#    - Compute medians per group
#    - Compute fold change (Pre / HC)
#    - FDR correction (Benjamini–Hochberg)
# -----------------------------
diff_result <- ptm_data_filled %>%
  mutate(
    row_id = row_number(),

    # Extract per-row group vectors (numeric)
    hc_values = map(row_id, ~ as.numeric(unlist(ptm_data_filled[.x, hc_cols]))),
    pre_values = map(row_id, ~ as.numeric(unlist(ptm_data_filled[.x, pre_cols]))),

    # Group medians
    hc_median  = map_dbl(hc_values,  ~ median(.x, na.rm = TRUE)),
    pre_median = map_dbl(pre_values, ~ median(.x, na.rm = TRUE)),

    # Wilcoxon rank-sum test (Pre vs HC)
    p_value = map2_dbl(hc_values, pre_values, ~ {
      # Note: x = Pre, y = HC
      wilcox.test(
        x = .y,
        y = .x,
        exact = FALSE,
        correct = TRUE
      )$p.value
    }),

    # Effect size proxies
    fold_change = pre_median / hc_median,
    fdr = p.adjust(p_value, method = "fdr"),

    # Significance label (based on FDR)
    significance = case_when(
      fdr < 0.001 ~ "***",
      fdr < 0.01  ~ "**",
      fdr < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  # Keep original order and arrange key columns first
  select(-row_id) %>%
  select(PTM.ProteinId, Gene, PTM.CollapseKey, everything())

# -----------------------------
# 7) Optional sanity check for a specific row index
# -----------------------------
test_idx <- 31

manual_hc  <- as.numeric(unlist(ptm_data_filled[test_idx, hc_cols]))
manual_pre <- as.numeric(unlist(ptm_data_filled[test_idx, pre_cols]))

cat("\nSanity check (row ", test_idx, ")\n", sep = "")
cat("Manual HC median:  ", median(manual_hc,  na.rm = TRUE), "\n", sep = "")
cat("Manual Pre median: ", median(manual_pre, na.rm = TRUE), "\n", sep = "")
cat("Code HC median:    ", diff_result$hc_median[test_idx], "\n", sep = "")
cat("Code Pre median:   ", diff_result$pre_median[test_idx], "\n", sep = "")
cat("P value:           ", diff_result$p_value[test_idx], "\n", sep = "")

# -----------------------------
# 8) Export results
#    Remove list-columns before writing to CSV
# -----------------------------
diff_result_export <- diff_result %>%
  select(-hc_values, -pre_values)

write.csv(
  diff_result_export,
  file = "Plasma_PTM_DifferentialAnalysis_HC_vs_Pre_v2.csv",
  row.names = FALSE
)
