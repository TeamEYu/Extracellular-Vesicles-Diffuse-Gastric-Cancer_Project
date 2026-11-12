# ============================================================
# 001_setup.R — Global setup (packages, reproducibility)
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(FactoMineR)
  library(factoextra)
  library(limma)
  library(remotes)
})

set.seed(12345)

close_all_devices <- function() {
  while (!is.null(dev.list())) dev.off()
}
