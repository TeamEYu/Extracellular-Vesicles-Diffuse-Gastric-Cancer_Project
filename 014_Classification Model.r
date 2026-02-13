# ============================================================
# Protein-based classification: cancer vs control
# Data: 27 cancer (Group=1) + 120 control (Group=0), total n=147
# Features: 4 core proteins (COL18A1, SERPINE2, CCT6A, GPNMB)
# Workflow:
#   1) Read and clean data (xlsx)
#   2) Convert types + impute missing values using control-group medians
#   3) Descriptive stats + Wilcoxon tests
#   4) Stratified train/test split
#   5) Standardize features based on training set
#   6) Single-protein logistic regression: AUROC, sensitivity, specificity, accuracy, best threshold
#   7) ROC plots (train vs test) per protein
#   8) Summarize model coefficients + AUCs; export to Excel
# ============================================================

rm(list = ls())
set.seed(123)

# -----------------------------
# 0) Packages
# -----------------------------
pkgs <- c("readxl", "dplyr", "caret", "pROC", "ggplot2", "openxlsx")
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
library(readxl)
library(dplyr)
library(caret)
library(pROC)
library(ggplot2)
library(openxlsx)

# -----------------------------
# 1) Input
# -----------------------------
xlsx_file <- "model_data.xlsx"
xlsx_sheet <- 2

# Expected columns:
#   Group (0/1) + 4 proteins
df <- read_excel(xlsx_file, sheet = xlsx_sheet) %>% as.data.frame()

# Basic checks
cat("Dimensions (rows, cols): ", paste(dim(df), collapse = " x "), "\n")
stopifnot("Group" %in% colnames(df))
cat("Group counts (raw):\n")
print(table(df$Group))

# Define feature names (must match your column names exactly)
feature_names <- c("COL18A1", "SERPINE2", "CCT6A", "GPNMB")
missing_features <- setdiff(feature_names, colnames(df))
if (length(missing_features) > 0) {
  stop("Missing feature columns in input: ", paste(missing_features, collapse = ", "))
}

# -----------------------------
# 2) Type conversion
# -----------------------------
# Convert Group to factor: control (0) as reference, cancer (1)
df$Group <- factor(df$Group, levels = c(0, 1), labels = c("control", "cancer"))

# Convert feature columns to numeric (robust to characters)
df <- df %>%
  mutate(across(all_of(feature_names), ~ suppressWarnings(as.numeric(as.character(.x)))))

cat("\nColumn classes:\n")
print(sapply(df[, c("Group", feature_names)], class))

# -----------------------------
# 3) Missing-value summary and imputation
# -----------------------------
missing_before <- sum(is.na(df[, feature_names]))
cat("\nTotal missing values in features before imputation: ", missing_before, "\n")

# Impute all missing values using the median of the control group for each protein
for (feat in feature_names) {
  control_vals <- df[df$Group == "control", feat, drop = TRUE]
  control_vals <- control_vals[!is.na(control_vals)]
  if (length(control_vals) == 0) stop("No valid values in control group for: ", feat)

  control_median <- median(control_vals)
  miss_idx <- is.na(df[[feat]])
  df[[feat]][miss_idx] <- control_median

  cat(
    "Imputed ", sum(miss_idx), " missing values for ", feat,
    " using control median = ", round(control_median, 4), "\n",
    sep = ""
  )
}

missing_after <- sum(is.na(df[, feature_names]))
if (missing_after == 0) {
  cat("\nAll missing values have been imputed successfully.\n")
} else {
  warning("Remaining missing values after imputation: ", missing_after)
}

cat("\nGroup counts (factor):\n")
print(table(df$Group))

# -----------------------------
# 4) Descriptive statistics + group comparisons
# -----------------------------
desc_stats <- df %>%
  group_by(Group) %>%
  summarise(
    across(all_of(feature_names),
           list(mean = ~ mean(.x), sd = ~ sd(.x)),
           .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

cat("\nDescriptive statistics (mean ± SD):\n")
print(desc_stats)

cat("\nWilcoxon rank-sum tests (each protein):\n")
wilcox_results <- lapply(feature_names, function(feat) {
  wt <- wilcox.test(df[[feat]] ~ df$Group, exact = FALSE)
  data.frame(
    Protein = feat,
    P_value = wt$p.value,
    stringsAsFactors = FALSE
  )
}) %>% bind_rows() %>%
  mutate(
    FDR = p.adjust(P_value, method = "fdr"),
    Significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

print(wilcox_results)

# -----------------------------
# 5) Stratified train/test split
# -----------------------------
# Use p=0.7 as in your code; change to 0.8 if you truly want 80/20.
trainIndex <- createDataPartition(df$Group, p = 0.7, list = FALSE, times = 1)
train_data <- df[trainIndex, ]
test_data  <- df[-trainIndex, ]

cat("\nTraining set counts:\n"); print(table(train_data$Group))
cat("Training set proportions:\n"); print(prop.table(table(train_data$Group)))

cat("\nTesting set counts:\n"); print(table(test_data$Group))
cat("Testing set proportions:\n"); print(prop.table(table(test_data$Group)))

# -----------------------------
# 6) Feature standardization (fit on training only)
# -----------------------------
pre_obj <- preProcess(train_data[, feature_names], method = c("center", "scale"))
train_features <- predict(pre_obj, train_data[, feature_names])
test_features  <- predict(pre_obj, test_data[, feature_names])

train_final <- data.frame(train_features, Group = train_data$Group)
test_final  <- data.frame(test_features,  Group = test_data$Group)

cat("\nTraining features after standardization (mean ~0, sd ~1):\n")
print(sapply(train_features, function(col) round(c(mean = mean(col), sd = sd(col)), 3)))

# -----------------------------
# 7) Single-protein models: logistic regression + performance
# -----------------------------
# Helper: ensure a complete 2x2 confusion matrix and compute metrics
compute_metrics <- function(actual_factor, pred_class) {
  cm <- table(Actual = actual_factor, Predicted = pred_class)

  # Ensure both levels exist in rows/cols
  all_levels <- c("control", "cancer")
  for (lvl in all_levels) {
    if (!lvl %in% rownames(cm)) cm <- rbind(cm, setNames(c(0, 0)[seq_len(ncol(cm))], colnames(cm)))
  }
  for (lvl in all_levels) {
    if (!lvl %in% colnames(cm)) cm <- cbind(cm, setNames(rep(0, nrow(cm)), lvl))
  }
  cm <- cm[all_levels, all_levels, drop = FALSE]

  TN <- cm["control", "control"]
  FP <- cm["control", "cancer"]
  FN <- cm["cancer",  "control"]
  TP <- cm["cancer",  "cancer"]

  sens <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  spec <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  acc  <- (TP + TN) / (TP + TN + FP + FN)

  list(
    CM = cm,
    TP = TP, TN = TN, FP = FP, FN = FN,
    Sensitivity = sens,
    Specificity = spec,
    Accuracy = acc
  )
}

# Result table (ENGLISH column names)
result_df <- data.frame(
  Protein = character(),
  Dataset = character(),
  AUROC = numeric(),
  Sensitivity = numeric(),
  Specificity = numeric(),
  Accuracy = numeric(),
  Best_Threshold = numeric(),
  TP = integer(),
  TN = integer(),
  FP = integer(),
  FN = integer(),
  stringsAsFactors = FALSE
)

for (feat in feature_names) {
  cat("\n==================================================\n")
  cat("Evaluating protein: ", feat, "\n", sep = "")

  # ---- Fit model on training set (single feature)
  train_single <- train_final[, c(feat, "Group")]
  train_single$Group <- relevel(train_single$Group, ref = "control")

  lr_single <- glm(Group ~ ., data = train_single, family = "binomial")

  # ---- Training predictions
  train_single$pred_prob <- predict(lr_single, newdata = train_single, type = "response")

  roc_train <- roc(
    response = as.numeric(train_single$Group) - 1,  # control=0, cancer=1
    predictor = train_single$pred_prob,
    quiet = TRUE
  )
  train_auc <- as.numeric(roc_train$auc)

  train_best_thr <- as.numeric(coords(roc_train, x = "best", ret = "threshold", transpose = FALSE))
  train_single$pred_class <- factor(
    ifelse(train_single$pred_prob > train_best_thr, "cancer", "control"),
    levels = c("control", "cancer")
  )

  m_train <- compute_metrics(train_single$Group, train_single$pred_class)

  cat("Training: AUC=", round(train_auc, 3),
      " | thr=", round(train_best_thr, 4),
      " | Sens=", round(m_train$Sensitivity, 3),
      " | Spec=", round(m_train$Specificity, 3),
      " | Acc=", round(m_train$Accuracy, 3), "\n", sep = "")
  cat("Training CM: TP=", m_train$TP, " TN=", m_train$TN, " FP=", m_train$FP, " FN=", m_train$FN, "\n", sep = "")

  result_df <- rbind(result_df, data.frame(
    Protein = feat,
    Dataset = "Training",
    AUROC = round(train_auc, 3),
    Sensitivity = round(m_train$Sensitivity, 3),
    Specificity = round(m_train$Specificity, 3),
    Accuracy = round(m_train$Accuracy, 3),
    Best_Threshold = train_best_thr,
    TP = m_train$TP, TN = m_train$TN, FP = m_train$FP, FN = m_train$FN,
    stringsAsFactors = FALSE
  ))

  # ---- Testing predictions (use SAME training-fitted model)
  test_single <- test_final[, c(feat, "Group")]
  test_single$pred_prob <- predict(lr_single, newdata = test_single, type = "response")

  roc_test <- roc(
    response = as.numeric(test_single$Group) - 1,
    predictor = test_single$pred_prob,
    quiet = TRUE
  )
  test_auc <- as.numeric(roc_test$auc)

  # Note: you previously used "best" threshold on the TEST ROC.
  # For more honest evaluation, use the TRAIN threshold on test.
  # Here we compute both; by default, we APPLY the TRAIN threshold to test.
  test_best_thr <- as.numeric(coords(roc_test, x = "best", ret = "threshold", transpose = FALSE))

  test_single$pred_class <- factor(
    ifelse(test_single$pred_prob > train_best_thr, "cancer", "control"),
    levels = c("control", "cancer")
  )
  m_test <- compute_metrics(test_single$Group, test_single$pred_class)

  cat("Testing:  AUC=", round(test_auc, 3),
      " | thr(train)=", round(train_best_thr, 4),
      " | Sens=", round(m_test$Sensitivity, 3),
      " | Spec=", round(m_test$Specificity, 3),
      " | Acc=", round(m_test$Accuracy, 3), "\n", sep = "")
  cat("Testing CM:  TP=", m_test$TP, " TN=", m_test$TN, " FP=", m_test$FP, " FN=", m_test$FN, "\n", sep = "")
  cat("Testing best threshold (not applied by default) = ", round(test_best_thr, 4), "\n", sep = "")

  result_df <- rbind(result_df, data.frame(
    Protein = feat,
    Dataset = "Testing",
    AUROC = round(test_auc, 3),
    Sensitivity = round(m_test$Sensitivity, 3),
    Specificity = round(m_test$Specificity, 3),
    Accuracy = round(m_test$Accuracy, 3),
    Best_Threshold = train_best_thr,  # store the applied threshold (train-derived)
    TP = m_test$TP, TN = m_test$TN, FP = m_test$FP, FN = m_test$FN,
    stringsAsFactors = FALSE
  ))
}

cat("\n\nFull performance table:\n")
print(result_df)

# -----------------------------
# 8) ROC plot function (train vs test) per protein
# -----------------------------
plot_roc_train_test <- function(feat, train_final, test_final, result_df,
                                out_pdf = NULL, width = 6.2, height = 5.7) {
  # Fit on training
  train_single <- train_final[, c(feat, "Group")]
  train_single$Group <- relevel(train_single$Group, ref = "control")
  model <- glm(Group ~ ., data = train_single, family = "binomial")

  train_single$pred_prob <- predict(model, newdata = train_single, type = "response")
  test_single <- test_final[, c(feat, "Group")]
  test_single$pred_prob <- predict(model, newdata = test_single, type = "response")

  roc_train <- roc(as.numeric(train_single$Group) - 1, train_single$pred_prob, quiet = TRUE)
  roc_test  <- roc(as.numeric(test_single$Group) - 1,  test_single$pred_prob,  quiet = TRUE)

  train_row <- result_df[result_df$Protein == feat & result_df$Dataset == "Training", ]
  test_row  <- result_df[result_df$Protein == feat & result_df$Dataset == "Testing", ]

  p <- ggroc(list(Training = roc_train, Testing = roc_test), aes = c("color")) +
    geom_line(linewidth = 1, alpha = 0.7) +
    annotate("text", x = 0.65, y = 0.45,
             label = paste0("Training AUC = ", train_row$AUROC),
             hjust = 0, color = "#218380", size = 6, fontface = "bold") +
    annotate("text", x = 0.65, y = 0.35,
             label = paste0("Training Accuracy = ", train_row$Accuracy),
             hjust = 0, color = "#218380", size = 6) +
    annotate("text", x = 0.65, y = 0.25,
             label = paste0("Testing AUC = ", test_row$AUROC),
             hjust = 0, color = "#B5179E", size = 6, fontface = "bold") +
    annotate("text", x = 0.65, y = 0.15,
             label = paste0("Testing Accuracy = ", test_row$Accuracy),
             hjust = 0, color = "#B5179E", size = 6) +
    scale_color_manual(values = c("#218380", "#B5179E")) +
    # ROC convention in pROC: x-axis is specificity by default for ggroc
    xlab("Specificity") +
    ylab("Sensitivity") +
    ggtitle(paste0(feat, " ROC Curve")) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.text = element_text(size = 16, color = "black"),
      axis.title = element_text(size = 18, color = "black"),
      panel.border = element_rect(linewidth = 1.2, color = "black"),
      legend.position = "none"
    )

  if (!is.null(out_pdf)) {
    ggsave(out_pdf, plot = p, width = width, height = height, device = "pdf")
  }
  p
}

# Generate ROC PDFs for all proteins
for (feat in feature_names) {
  p <- plot_roc_train_test(
    feat = feat,
    train_final = train_final,
    test_final = test_final,
    result_df = result_df,
    out_pdf = paste0(feat, "_ROC.pdf")
  )
  print(p)
}

# -----------------------------
# 9) Summarize model coefficients + AUCs and export
# -----------------------------
model_params <- data.frame(
  Protein = character(),
  Model_Equation = character(),
  Intercept = numeric(),
  Intercept_SE = numeric(),
  Intercept_P = numeric(),
  Coefficient = numeric(),
  Coefficient_SE = numeric(),
  Coefficient_P = numeric(),
  Coef_Significance = character(),
  Training_AUC = numeric(),
  Testing_AUC = numeric(),
  stringsAsFactors = FALSE
)

for (prot in feature_names) {
  train_data_prot <- train_final[, c(prot, "Group")]
  train_data_prot$Group <- relevel(train_data_prot$Group, ref = "control")
  model_prot <- glm(Group ~ ., data = train_data_prot, family = "binomial")

  coef_tab <- summary(model_prot)$coefficients

  intercept_val <- coef_tab["(Intercept)", "Estimate"]
  intercept_se  <- coef_tab["(Intercept)", "Std. Error"]
  intercept_p   <- coef_tab["(Intercept)", "Pr(>|z|)"]

  prot_coef     <- coef_tab[prot, "Estimate"]
  prot_coef_se  <- coef_tab[prot, "Std. Error"]
  prot_coef_p   <- coef_tab[prot, "Pr(>|z|)"]

  coef_sig <- dplyr::case_when(
    prot_coef_p < 0.001 ~ "***",
    prot_coef_p < 0.01  ~ "**",
    prot_coef_p < 0.05  ~ "*",
    TRUE ~ "ns"
  )

  # Note: features are standardized => coefficient is per 1 SD increase
  model_equation <- sprintf(
    "logit(P(cancer)) = %.4f + %.4f × %s (standardized)",
    intercept_val, prot_coef, prot
  )

  train_auc <- result_df[result_df$Protein == prot & result_df$Dataset == "Training", "AUROC"]
  test_auc  <- result_df[result_df$Protein == prot & result_df$Dataset == "Testing",  "AUROC"]

  model_params <- rbind(model_params, data.frame(
    Protein = prot,
    Model_Equation = model_equation,
    Intercept = round(intercept_val, 4),
    Intercept_SE = round(intercept_se, 4),
    Intercept_P = signif(intercept_p, 4),
    Coefficient = round(prot_coef, 4),
    Coefficient_SE = round(prot_coef_se, 4),
    Coefficient_P = signif(prot_coef_p, 4),
    Coef_Significance = coef_sig,
    Training_AUC = as.numeric(train_auc),
    Testing_AUC = as.numeric(test_auc),
    stringsAsFactors = FALSE
  ))
}

cat("\n==================================================\n")
cat("Single-protein logistic model parameter summary\n")
cat("==================================================\n")
print(model_params, row.names = FALSE)

# Export to Excel
write.xlsx(model_params, "SingleProtein_ModelParameter_Summary.xlsx", rowNames = FALSE)
cat("\nSaved: SingleProtein_ModelParameter_Summary.xlsx\n")