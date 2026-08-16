# Pan-Cancer Proteomics Pipeline

[![R](https://img.shields.io/badge/R-4.3%2B-blue.svg)](https://www.r-project.org/)
[![Python](https://img.shields.io/badge/Python-3.10%2B-green.svg)](https://www.python.org/)
[![Status](https://img.shields.io/badge/Status-Completed-success.svg)]()

An end-to-end bioinformatics pipeline for cross-tumor proteomic processing, differential expression analysis, multi-class machine learning classification, and functional pathway enrichment across Breast (**BRCA**), Lung (**LUAD**), and Colorectal (**COAD**) carcinomas using CPTAC datasets.

---

## Key Findings

- **Cross-Cohort Proteomic Alignment:** Intersected and normalized **6,633 common proteins** across 3 tumor types after handling pipe-delimited accession headers, k-NN imputation, and quantile normalization.
- **Machine Learning Classification:** Stratified Random Forest classifier achieved **97% overall accuracy** (Macro F1 = 0.97) in distinguishing BRCA, LUAD, and COAD cohorts.
- **Key Biomarker Drivers:** Identified diagnostic proteomic features including `EHD3`, `MYL6`, `CAV1`, `CHGB`, `ERBB2` (HER2), `MYL9`, and `TAGLN`.
- **Enriched Functional Pathways:** Biomarkers mapped significantly to **Sema4D/Semaphorin signaling** ($p = 5.44 \times 10^{-6}$), **Actomyosin dynamics** ($p = 6.54 \times 10^{-6}$), **RHO GTPase activation**, and **Focal Adhesion**.

---

## 📁 Repository Structure

├── data/
│   ├── raw/                      # Raw CPTAC quantification files (.tsv)
│   └── processed/                # Normalized pan-cancer matrix
├── scripts/
│   ├── 01_preprocessing.R        # Identifier cleaning & cross-cohort matrix alignment
│   ├── 02_analysis.R             # PCA & limma differential expression
│   ├── 03_visualization.R        # Volcano plots & hierarchical heatmaps
│   ├── 04_ml_classifier.py       # Random Forest multi-class model & feature ranking
│   └── 05_pathway_enrichment.R   # g:Profiler functional enrichment (GO/KEGG/Reactome)
└── results/                      # Generated figures (.png) and tabular outputs (.csv)

---

## 🛠️ Pipeline Execution

1. **Preprocessing & Quantile Normalization:**
   ```bash
   Rscript scripts/01_preprocessing.R
