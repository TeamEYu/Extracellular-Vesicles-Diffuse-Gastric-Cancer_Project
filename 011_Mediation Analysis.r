# ============================================================
# Mediation analysis (three-step approach) with Firth logistic
# regression and bootstrap inference for the indirect effect
# Dataset: plasma significant PTM mediation data
# ============================================================

# -----------------------------
# 0) Environment setup
# -----------------------------
rm(list = ls())

# Set working directory (edit as needed)
setwd("D:/MAS_GastricCancer_Project/Company_Results/Protein_Modifications")

# -----------------------------
# 1) Load required packages
# -----------------------------
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("logistf", quietly = TRUE)) install.packages("logistf")
if (!requireNamespace("boot", quietly = TRUE)) install.packages("boot")

library(readxl)
library(dplyr)
library(logistf)
library(boot)

# -----------------------------
# 2) Read data and preprocess
# -----------------------------
medidata <- read_excel("Plasma_SignificantModification_MediationData_Supplement.xlsx")

# Standardize variable names
# - Q15828: protein abundance (exposure)
# - Q15828M3: modified peptide / PTM feature (mediator)
medidata <- medidata %>%
  rename(
    Q15828   = CRISPLD2,
    Q15828M3 = CRISPLD2_M80_M2
  )

# Convert selected columns to numeric (adjust the column range if needed)
medidata[, 3:69] <- lapply(medidata[, 3:69], function(x) as.numeric(as.character(x)))

# Log2-transform mediator
# If zeros are possible, use log2(x + 1) instead of log2(x)
medidata$Q15828M3 <- log2(medidata$Q15828M3)

# Confirm outcome y is binary (0/1)
print(table(medidata$y))

# -----------------------------
# 3) Three-step mediation models
# -----------------------------
# Model 1: total effect (c) of exposure on outcome using Firth logistic regression
model1_firth <- logistf(
  formula = y ~ Q15828,
  data = medidata,
  pl = TRUE,
  plcontrol = logistpl.control(maxit = 1000)
)
print(summary(model1_firth))

# Model 2: effect (a) of exposure on mediator using linear regression
model2_lm <- lm(Q15828M3 ~ Q15828, data = medidata)
print(summary(model2_lm))

# Model 3: direct effect (c') of exposure on outcome adjusting for mediator,
# and mediator effect (b), using Firth logistic regression
model3_firth <- logistf(
  formula = y ~ Q15828 + Q15828M3,
  data = medidata,
  pl = TRUE,
  control = logistf.control(maxit = 5000, maxstep = 0.05)
)
print(summary(model3_firth))

# -----------------------------
# 4) Bootstrap inference for indirect effect (a*b)
#    Compatible with Firth logistic regression in Model 3
# -----------------------------
boot_fun <- function(data, indices) {
  d <- data[indices, , drop = FALSE]

  # Path a: exposure -> mediator
  a_model <- lm(Q15828M3 ~ Q15828, data = d)
  a <- coef(a_model)[["Q15828"]]

  # Path b: mediator -> outcome (adjusting for exposure), Firth logistic regression
  b_model <- logistf(
    y ~ Q15828 + Q15828M3,
    data = d,
    pl = FALSE,
    control = logistf.control(maxit = 5000, maxstep = 0.05),
    plcontrol = logistpl.control(maxit = 5000, eps = 1e-6)
  )
  b <- coef(b_model)[["Q15828M3"]]

  return(a * b)
}

set.seed(123)
boot_result <- boot(
  data = medidata,
  statistic = boot_fun,
  R = 1000
)

# Remove failed resamples (NA) while preserving matrix structure
valid_idx <- !is.na(boot_result$t)
boot_result$t <- boot_result$t[valid_idx, , drop = FALSE]

# BCa confidence interval for the indirect effect
boot_ci <- boot.ci(boot_result, type = "bca")

# -----------------------------
# 5) Extract core effects and report
# -----------------------------
# Direct effect (ADE): coefficient of Q15828 in Model 3
ade <- coef(model3_firth)[["Q15828"]]
ade_p <- model3_firth$prob[["Q15828"]]

# Total effect: coefficient of Q15828 in Model 1
total_effect <- coef(model1_firth)[["Q15828"]]
total_p <- model1_firth$prob[["Q15828"]]

# Indirect effect (ACME): mean of bootstrap a*b
acme_est <- mean(boot_result$t)
acme_ci <- boot_ci$bca[4:5]

# Proportion mediated (on the log-odds scale; interpret cautiously)
prop_mediated <- (acme_est / total_effect) * 100

cat("\nMediation analysis results:\n")
cat("1) Indirect effect (a*b): ", acme_est,
    " | 95% CI (BCa): ", acme_ci[1], ", ", acme_ci[2], "\n", sep = "")
cat("2) Direct effect (c'):   ", ade,
    " | P value: ", ade_p, "\n", sep = "")
cat("3) Total effect (c):     ", total_effect,
    " | P value: ", total_p, "\n", sep = "")
cat("4) Proportion mediated:  ", round(prop_mediated, 2), "%\n", sep = "")

# -----------------------------
# Notes
# -----------------------------
# - The proportion mediated computed as (a*b)/c is on the log-odds scale because
#   Models 1 and 3 are logistic; report and interpret it accordingly.
# - If Q15828M3 contains zeros/negatives, consider using log2(Q15828M3 + 1) or
#   another offset to avoid -Inf/NaN.