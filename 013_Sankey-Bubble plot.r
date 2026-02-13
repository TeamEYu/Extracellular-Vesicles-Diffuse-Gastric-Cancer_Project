# ============================================================
# Integrated Sankey + Dotplot visualization (Metascape/KEGG-like)
# - Input: one Excel file containing pathway table (sheet 5) and gene logFC table (sheet 6)
# - Output: a custom-layout PDF combining:
#     (1) Sankey diagram: Pathway -> Gene, gene nodes colored by logFC
#     (2) Dotplot: GeneRatio vs pathway position, colored/sized by -log10(pvalue)
# ============================================================

rm(list = ls())

# -----------------------------
# 0) Environment setup
# -----------------------------
setwd("D:/MAS_GastricCancer_Project/DGC_DEP_LW/v8/Public_Datasets/ACRG_Cohort_GSE62254")

# -----------------------------
# 1) Load required packages
# -----------------------------
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("ggsankey", quietly = TRUE)) install.packages("ggsankey")
if (!requireNamespace("grid", quietly = TRUE)) install.packages("grid")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsankey)
library(grid)
library(gridExtra)

# -----------------------------
# 2) Read input tables
# -----------------------------
xlsx_file <- "metascape_tissue_specific_93.txji85f6l.xlsx"

pathway_tbl <- read_excel(xlsx_file, sheet = 5)
gene_logfc  <- read_excel(xlsx_file, sheet = 6)

# Quick checks (optional)
head(pathway_tbl)
head(gene_logfc)
str(pathway_tbl)
str(gene_logfc)

# -----------------------------
# 3) Preprocess: expand Gene_Names into one gene per row and merge logFC
#    Assumption:
#      pathway_tbl columns: Description, Gene_Names (comma-separated)
#      gene_logfc columns: Gene, logFC (or a numeric column to color gene nodes)
# -----------------------------
kegg_data <- pathway_tbl %>%
  select(1, 2) %>%                         # Description + Gene_Names
  rename(
    Description = 1,
    Gene_Names  = 2
  ) %>%
  separate_rows(Gene_Names, sep = ",") %>%
  mutate(Gene_Names = trimws(Gene_Names)) %>%
  left_join(gene_logfc, by = c("Gene_Names" = "Gene"))

# -----------------------------
# 4) Build Sankey data
#    - Rename to: path (pathway), element (gene), num (logFC or score)
# -----------------------------
sankey <- kegg_data %>%
  transmute(
    path    = Description,
    element = Gene_Names,
    num     = suppressWarnings(as.numeric(as.character(logFC)))
  ) %>%
  arrange(path) %>%
  mutate(path = factor(path, levels = unique(path)))

# -----------------------------
# 5) Define node colors
#    - Pathway nodes: distinct categorical palette (extended if needed)
#    - Gene nodes: gradient by num (logFC)
# -----------------------------
# Gene node gradient: blue -> light gray -> orange
gene_palette_fun <- colorRampPalette(c("#008bd0", "#eeeeee", "#ffa61d"))

# Pathway node palette (provide enough distinct colors)
path_palette <- c(
  "#00a08f", "#dbc251", "#f8b251", "#e88e8f", "#9a99e1",
  "#008bd0", "#55b89e", "#f09c96", "#d4b97f", "#8ecae6",
  "#6a4c93", "#ff6f61", "#2a9d8f", "#e9c46a", "#264653"
)

pathway_colors <- tibble(
  path  = levels(sankey$path),
  color = path_palette[seq_len(length(levels(sankey$path)))]
)

# Map gene colors by binning num into 50 buckets
# NOTE: adjust breaks to reflect your actual logFC range
gene_color_map <- sankey %>%
  mutate(
    num = as.numeric(num),
    color = gene_palette_fun(50)[as.integer(cut(
      num,
      breaks = seq(-2, 2, length.out = 51),  # recommended for typical logFC
      include.lowest = TRUE
    ))]
  ) %>%
  select(element, color) %>%
  distinct()

# -----------------------------
# 6) Convert to ggsankey long format and attach colors
# -----------------------------
sankey_long <- sankey %>%
  make_long(path, element) %>%
  mutate(
    color = ifelse(
      node %in% sankey$element,
      gene_color_map$color[match(node, gene_color_map$element)],
      pathway_colors$color[match(node, pathway_colors$path)]
    )
  )

# -----------------------------
# 7) Draw the Sankey diagram
# -----------------------------
p_sankey <- ggplot(
  sankey_long,
  aes(
    x = x,
    next_x = next_x,
    node = node,
    next_node = next_node,
    fill = color,
    label = node
  )
) +
  scale_fill_identity() +
  geom_sankey(
    flow.alpha = 0.2,
    smooth = 8,
    width = 0.08
  ) +
  geom_sankey_text(size = 4, color = "black", hjust = 0) +
  theme_void() +
  theme(legend.position = "none")

# -----------------------------
# 8) Extract pathway node Y-positions from ggplot build
#    We keep x == 1 (left-side pathway nodes) to align dotplot.
# -----------------------------
tmp_build <- ggplot_build(p_sankey)
node_layer <- NULL

for (i in seq_along(tmp_build$data)) {
  if (all(c("x", "ymin", "ymax", "node") %in% names(tmp_build$data[[i]]))) {
    node_layer <- tmp_build$data[[i]]
    break
  }
}

stopifnot(!is.null(node_layer))

sankey_nodes <- tibble(
  node = node_layer$node,
  x = node_layer$x,
  ymin = node_layer$ymin,
  ymax = node_layer$ymax,
  y_center = (node_layer$ymin + node_layer$ymax) / 2
) %>%
  filter(x == 1)

# Attach y_center to the original pathway table for dotplot alignment
# Assumption: pathway_tbl has a column named "Description"
pathway_tbl_aligned <- pathway_tbl %>%
  left_join(sankey_nodes %>% select(node, y_center), by = c("Description" = "node"))

# -----------------------------
# 9) Ensure GeneRatio is numeric
#    Convert "a/b" to a numeric ratio (a / b)
# -----------------------------
pathway_tbl_aligned <- pathway_tbl_aligned %>%
  separate(
    col = GeneRatio,
    into = c("num", "den"),
    sep = "/",
    remove = FALSE,
    fill = "right"
  ) %>%
  mutate(
    num = suppressWarnings(as.numeric(num)),
    den = suppressWarnings(as.numeric(den)),
    GeneRatio = num / den
  ) %>%
  select(-num, -den)

# Optional checks
class(pathway_tbl_aligned$GeneRatio)
table(is.na(pathway_tbl_aligned$GeneRatio))

# -----------------------------
# 10) Dotplot (bubble plot)
#     - x: y_center (aligned to Sankey pathway nodes)
#     - y: GeneRatio
#     - size/color: -log10(pvalue)
# -----------------------------
dot_palette <- colorRampPalette(c("#393781", "#eeeeee", "#f22942"))(100)

p_dot <- ggplot(pathway_tbl_aligned, aes(x = y_center, y = GeneRatio)) +
  coord_flip() +
  geom_point(
    aes(size = -log10(pvalue), color = -log10(pvalue)),
    shape = 16
  ) +
  scale_size_continuous(range = c(2, 8)) +
  scale_color_gradientn(colors = dot_palette, name = "-log10(pvalue)") +
  labs(x = NULL, y = "GeneRatio") +
  guides(size = "none") +
  scale_y_continuous(
    limits = c(
      min(pathway_tbl_aligned$GeneRatio, na.rm = TRUE) - 0.01,
      max(pathway_tbl_aligned$GeneRatio, na.rm = TRUE) + 0.01
    )
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text = element_text(color = "black", size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.position = "left"
  )

# -----------------------------
# 11) Custom layout export (Sankey + Dotplot) to PDF
# -----------------------------
pdf("custom_layout_v2.pdf", width = 8.5, height = 5.5)

grid.newpage()
pushViewport(viewport(layout = grid.layout(
  nrow = 70,
  ncol = 100,
  heights = unit(rep(1, 70), "null"),
  widths  = unit(rep(1, 100), "null")
)))

# Sankey main panel
pushViewport(viewport(
  layout.pos.row = 1:69,
  layout.pos.col = 35:99,
  name = "sankey_viewport"
))
print(p_sankey, newpage = FALSE)
popViewport()

# Dotplot inset panel
pushViewport(viewport(
  layout.pos.row = 18:57,
  layout.pos.col = 15:50,
  name = "dotplot_viewport"
))
print(p_dot, newpage = FALSE)
popViewport()

dev.off()