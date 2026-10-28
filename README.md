# Extracellular-Vesicles-Diffuse-Gastric-Cancer_Project

## Study Overview

We designed a multi-tier pipeline to identify **DGC-specific, tumor-traceable EV biomarkers** with improved specificity by:

1. Discovery in plasma EVs (DGC vs HC; DGC vs HP; intersection)
2. Benign-stability check (HP vs HC)
3. Subtype-specific filtering (remove EVDEPs shared with IGC)
4. Tissue traceability (tumor-enriched in DGCT vs DNT/HPT; stable baseline DNT vs HPT)
5. Plasma–tissue concordance and directional consistency
6. Network-level prioritization (STRING + Cytoscape/cytoHubba)
7. External validation (TCGA-STAD; GEO-ACRG survival)
8. Dynamic monitoring (paired pre-/post-operative plasma in DGC)
9. Mechanistic linkage (cell high-salt EV Analysis)

## Repository Structure
├── 001_setup.r
├── 002_preprocess_matrix.r
├── 003_qc_plots.r
├── 004_EV Protein Identity.r
├── 005_run_plasma_or_tissue.r
├── 006_heatmap.r
├── 007_GEO Survival Analysis.r
├── 008_GO Analysis.r
├── 009_TCGA Analysis.r
├── 010_TMN-box plot Analysis.r
├── 011_Mediation Analysis.r
├── 012_Post-translational-modificati...
├── 013_Sankey-Bubble plot.r
├── 014_Classification Model.r
├── 015_CellExperiment_Transwell An...
├── 016_CellExperiment_Woundhealin...
├── 017_CellExperiment-Omics Analys...
├── 018_WesternBlot-Analysis.r
└── README.md

## Requirements

### R (primary analysis)
- R >= 4.3.1
- Recommended packages:
  - `tidyverse`, `data.table`
  - `limma`, `edgeR` (if needed)
  - `DreamAI` (imputation)
  - `pROC`
  - `survival`, `survminer`
  - `WGCNA`
  - `clusterProfiler`
  - `rtracklayer`
  - `ComplexHeatmap` or `pheatmap`

### External software / web tools
- Spectronaut (directDIA / PTM workflows)
- Skyline (PRM curation)
- Cytoscape + cytoHubba
- STRING (v12.0 used for PPI export)

---Done
