# scripts/02_analysis.R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")

library(tidyverse)
library(limma)

cat("Loading normalized proteomics dataset...\n")
data <- read.csv("data/processed/normalized_proteomics.csv", row.names = 1, check.names = FALSE)

# 1. PCA Analysis
cat("Performing Principal Component Analysis (PCA)...\n")
pca_res <- prcomp(t(data), scale. = TRUE)

# Metadata parsing
sample_names <- colnames(data)
cancer_types <- sub("_.*$", "", sample_names)
pca_df <- data.frame(
  Sample = sample_names,
  CancerType = cancer_types,
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2]
)

# Plotting PCA
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = CancerType)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "Pan-Cancer Proteomic PCA Profile", x = "PC1", y = "PC2")

ggsave("results/pca_plot.png", plot = p_pca, width = 7, height = 5)
cat("PCA plot saved to results/pca_plot.png\n")

# 2. Differential Expression Analysis (limma)
cat("Running Differential Expression Analysis via limma...\n")
design <- model.matrix(~ 0 + factor(cancer_types))
colnames(design) <- levels(factor(cancer_types))

fit <- lmFit(data, design)
contrast_matrix <- makeContrasts(
  BRCA_vs_Others = BRCA - (LUAD + COAD)/2,
  LUAD_vs_Others = LUAD - (BRCA + COAD)/2,
  COAD_vs_Others = COAD - (BRCA + LUAD)/2,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- dBayer <- eBayes(fit2)

# Extracting top DE proteins for BRCA
de_brca <- topTable(fit2, coef = "BRCA_vs_Others", number = Inf, adjust.method = "BH")
write.csv(de_brca, "results/de_brca_vs_others.csv")

cat("Differential expression complete. Results saved to results/de_brca_vs_others.csv\n")
