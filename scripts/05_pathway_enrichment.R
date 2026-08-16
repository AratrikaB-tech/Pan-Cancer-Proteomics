# scripts/05_pathway_enrichment.R
if (!requireNamespace("gprofiler2", quietly = TRUE)) install.packages("gprofiler2", repos = "https://cloud.r-project.org")

library(tidyverse)
library(gprofiler2)

cat("Loading top biomarker proteins from machine learning model...\n")
biomarkers <- read.csv("results/top20_biomarker_features.csv", header = TRUE)
gene_list <- biomarkers[[1]]

cat(paste("Running pathway enrichment for", length(gene_list), "top biomarkers...\n"))

# Perform functional profiling against Homo sapiens databases
gost_res <- gost(
  query = gene_list,
  organism = "hsapiens",
  sources = c("GO:BP", "KEGG", "REAC"),
  user_threshold = 0.05,
  correction_method = "g_SCS"
)

if (!is.null(gost_res$result)) {
  enrichment_df <- gost_res$result %>%
    select(source, term_id, term_name, p_value, term_size, intersection_size) %>%
    arrange(p_value)
  
  write.csv(enrichment_df, "results/pathway_enrichment_results.csv", row.names = FALSE)
  cat("Pathways successfully mapped! Saved to results/pathway_enrichment_results.csv\n")
  
  # Bar plot of top 10 enriched terms
  top_terms <- head(enrichment_df, 10)
  p_terms <- ggplot(top_terms, aes(x = reorder(term_name, -p_value), y = -log10(p_value), fill = source)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top Enriched Biological Pathways", x = "Pathway Term", y = "-Log10 p-value")
  
  ggsave("results/pathway_enrichment_plot.png", plot = p_terms, width = 8, height = 5)
  cat("Pathway plot saved to results/pathway_enrichment_plot.png\n")
} else {
  cat("No statistically significant pathways detected at p < 0.05 threshold.\n")
}
