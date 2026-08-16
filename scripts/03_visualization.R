# scripts/03_visualization.R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")

library(tidyverse)
library(limma)
library(pheatmap)

cat("Loading normalized expression data and running contrasts...\n")
data <- read.csv("data/processed/normalized_proteomics.csv", row.names = 1, check.names = FALSE)

sample_names <- colnames(data)
cancer_types <- sub("_.*$", "", sample_names)

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
fit2 <- eBayes(fit2)

# 1. Exporting DE Tables for LUAD and COAD
de_luad <- topTable(fit2, coef = "LUAD_vs_Others", number = Inf, adjust.method = "BH")
de_coad <- topTable(fit2, coef = "COAD_vs_Others", number = Inf, adjust.method = "BH")

write.csv(de_luad, "results/de_luad_vs_others.csv")
write.csv(de_coad, "results/de_coad_vs_others.csv")
cat("Exported DE tables for LUAD and COAD.\n")

# 2. Volcano Plot for BRCA vs Others
de_brca <- topTable(fit2, coef = "BRCA_vs_Others", number = Inf, adjust.method = "BH")
de_brca$Gene <- rownames(de_brca)
de_brca$Significance <- case_when(
  de_brca$logFC > 1 & de_brca$adj.P.Val < 0.05 ~ "Up-regulated",
  de_brca$logFC < -1 & de_brca$adj.P.Val < 0.05 ~ "Down-regulated",
  TRUE ~ "Not Significant"
)

p_volcano <- ggplot(de_brca, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up-regulated" = "#d95f02", "Down-regulated" = "#7570b3", "Not Significant" = "grey70")) +
  theme_minimal() +
  labs(title = "Volcano Plot: BRCA vs. Pan-Cancer Cohort", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value")

ggsave("results/volcano_brca.png", plot = p_volcano, width = 7, height = 5)
cat("Volcano plot saved to results/volcano_brca.png\n")

# 3. Expression Heatmap for Top 50 Variable Proteins
top_proteins <- head(rownames(de_brca[order(de_brca$adj.P.Val), ]), 50)
heatmap_mat <- data[top_proteins, ]

# Annotation dataframe
annotation_col <- data.frame(CancerType = cancer_types)
rownames(annotation_col) <- sample_names

png("results/heatmap_top50_proteins.png", width = 900, height = 800, res = 120)
pheatmap(heatmap_mat,
         scale = "row",
         show_colnames = FALSE,
         annotation_col = annotation_col,
         main = "Top 50 Differentially Expressed Proteomic Markers")
dev.off()
cat("Heatmap saved to results/heatmap_top50_proteins.png\n")
