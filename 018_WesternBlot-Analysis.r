# ============================================================
# Full R script: per-protein HC vs DGC barplot + mean±SD
# Welch t-test; stars: *** <0.001, ** <0.01, * <0.05
# Colors: same hue (light vs dark); brackets per protein pair
# ============================================================

# -----------------------------
# 0) Packages
# -----------------------------
# install.packages(c("tidyverse", "ggpubr", "rstatix"))  # run once if needed
library(tidyverse)
library(ggpubr)
library(rstatix)

# -----------------------------
# 1) Input data (from your sheet)
# -----------------------------
df <- tribble(
  ~pos, ~group, ~protein, ~value,
  
 ******
)

# Ensure order: HC left, DGC right; protein order same as your figure
df <- df %>%
  mutate(
    group = factor(group, levels = c("HC", "DGC")),
    protein = factor(protein, levels = c("CCT6A","COL18A1","GPNMB","SERPINE2","SYNGR1","SYPL1"))
  )

# -----------------------------
# 2) Summary: mean ± SD
# -----------------------------
df_summary <- df %>%
  group_by(protein, group) %>%
  summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value),
    .groups = "drop"
  )

# -----------------------------
# 3) Welch t-test per protein + star labels (your thresholds)
# -----------------------------
stat_tbl <- df %>%
  group_by(protein) %>%
  t_test(value ~ group, var.equal = FALSE) %>%   # Welch t-test
  mutate(
    p.signif = case_when(
      p < 1e-3 ~ "***",
      p < 1e-2 ~ "**",
      p < 5e-2 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  add_xy_position(x = "protein", dodge = 0.70)    # puts brackets above each pair

# -----------------------------
# 4) Plot: same-hue light vs dark (clean, publication-safe)
# -----------------------------
# Blue-grey: light vs dark (recommended)
cols <- c("HC" = "#EAD6D1", "DGC" = "#9C6B63")

p <- ggplot(df_summary, aes(x = protein, y = mean, fill = group)) +
  geom_col(
    position = position_dodge(width = 0.70),
    width = 0.62,
    color = "black",
    linewidth = 0.6
  ) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.70),
    width = 0.18,
    linewidth = 0.6
  ) +
  stat_pvalue_manual(
    stat_tbl,
    label = "p.signif",
    tip.length = 0.01,
    bracket.size = 0.6,
    size = 4.2,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = cols) +
  labs(x = NULL, y = "ImageJ area (a.u.)", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 1.4, colour = "black"),
    axis.ticks = element_line(linewidth = 1.2, colour = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    legend.position = "right"
  )

print(p)

# -----------------------------
# 5) Save figure (high-res)
# -----------------------------
ggsave("ImageJ_area_HC_vs_DGC_stars.png",
       p,
       width = 6.6,
       height = 5,
       dpi = 600)

# -----------------------------
# 6) Optional: export a summary + p table
# -----------------------------
out_tbl <- df_summary %>%
  mutate(mean_sd = sprintf("%.2f ± %.2f", mean, sd)) %>%
  select(protein, group, mean_sd) %>%
  pivot_wider(names_from = group, values_from = mean_sd) %>%
  left_join(stat_tbl %>% select(protein, p, p.signif), by = "protein")

print(out_tbl)