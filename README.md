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

│   ├── 00_config.R
│   ├── 01_qc_plasma.R
│   ├── 02_imputation_dreamai.R
│   ├── 03_limma_plasma_DE.R
│   ├── 04_subtype_filter_IGC.R
│   ├── 05_tissue_traceability.R
│   ├── 06_concordance_plasma_tissue.R
│   ├── 07_prm_validation.R
│   ├── 08_tcga_processing.R
│   ├── 09_geo_acrg_survival.R
│   ├── 10_prepost_surgery_limma.R
│   ├── 11_network_string_cytoscape.R
│   ├── 12_wgcna_diet_modules.R
│   └── 13_ptm_mediation.R

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
