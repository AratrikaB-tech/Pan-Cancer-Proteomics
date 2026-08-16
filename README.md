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

```text
Pan-Cancer-Proteomics/
├── README.md
├── data/
│   └── processed/
│       └── normalized_proteomics.csv
├── scripts/
│   ├── 01_preprocessing.R
│   ├── 02_analysis.R
│   ├── 03_visualization.R
│   ├── 04_ml_classifier.py
│   └── 05_pathway_enrichment.R
└── results/
    ├── pca_plot.png
    ├── volcano_brca.png
    ├── heatmap_top50_proteins.png
    ├── confusion_matrix.png
    ├── top20_biomarkers.png
    ├── pathway_enrichment_plot.png
    ├── de_brca_vs_others.csv
    ├── de_luad_vs_others.csv
    ├── de_coad_vs_others.csv
    ├── top20_biomarker_features.csv
    └── pathway_enrichment_results.csv

```

## 🛠️ Pipeline Execution

1. **Preprocessing & Quantile Normalization:**
   ```bash
   Rscript scripts/01_preprocessing.R

2. **Differential Expression & PCA:**
   ```bash
   Rscript scripts/02_analysis.R

3. **Generating Visualizations (Heatmaps & Volcano Plots):**
   ```bash
   Rscript scripts/03_visualization.R

4. **Multi-Class Machine Learning Classifier:**
   ```bash
   python3 scripts/04_ml_classifier.py

5. **Functional Pathway Enrichment Analysis:**
   ```bash
   Rscript scripts/05_pathway_enrichment.R
