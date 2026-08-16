# scripts/01_preprocessing.R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("impute", quietly = TRUE)) BiocManager::install("impute")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")

library(tidyverse)
library(impute)
library(limma)

cat("Loading raw CPTAC datasets from data/raw/...\n")

# Helper function to parse pipe-delimited CPTAC tables
load_cptac <- function(filepath, prefix) {
  df <- read.delim(filepath, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

  # Extracting gene symbol prior to the '|' character
  raw_ids <- df[, 1]
  gene_symbols <- sub("\\|.*$", "", raw_ids)

  # Extracting expression data
  sample_cols <- colnames(df)[-1]
  mat <- as.matrix(df[, sample_cols])
  mat <- apply(mat, 2, as.numeric)

  # Combining and aggregating duplicate gene symbols by taking the mean
  df_clean <- as.data.frame(mat)
  df_clean$Gene <- gene_symbols

  df_aggregated <- df_clean %>%
    filter(!is.na(Gene) & Gene != "" & Gene != "Composite.Element.REF") %>%
    group_by(Gene) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  res_mat <- as.matrix(df_aggregated[, -1])
  rownames(res_mat) <- df_aggregated$Gene
  colnames(res_mat) <- paste0(prefix, "_", colnames(res_mat))

  return(res_mat)
}

# Loading cohorts
brca <- load_cptac("data/raw/brca_cptac.tsv", "BRCA")
luad <- load_cptac("data/raw/luad_cptac.tsv", "LUAD")
coad <- load_cptac("data/raw/coad_cptac.tsv", "COAD")

# Intersecting common gene symbols across all 3 cohorts
common_genes <- reduce(list(rownames(brca), rownames(luad), rownames(coad)), intersect)
cat(paste("Found", length(common_genes), "common gene symbols across all 3 datasets.\n"))

# Merging expression matrices
merged_matrix <- cbind(
  brca[common_genes, ],
  luad[common_genes, ],
  coad[common_genes, ]
)

# 1. Filtering out proteins missing in > 30% of total samples
n_samples <- ncol(merged_matrix)
keep_proteins <- rowSums(is.na(merged_matrix)) < (0.3 * n_samples)
filtered_matrix <- merged_matrix[keep_proteins, ]

# 2. k-NN Imputation for missing values
imputed_data <- impute.knn(as.matrix(filtered_matrix))$data

# 3. Quantile Normalization
norm_matrix <- normalizeBetweenArrays(imputed_data)

# 4. Saving Processed Matrix
write.csv(norm_matrix, "data/processed/normalized_proteomics.csv")
cat("Preprocessing complete! Saved to data/processed/normalized_proteomics.csv\n")
