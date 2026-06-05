#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

dir.create("figures/osmoadaptation", recursive = TRUE, showWarnings = FALSE)
dir.create("results/osmoadaptation", recursive = TRUE, showWarnings = FALSE)

df <- read_tsv(
  "results/tables/Table3_Osmoadaptation_Candidates.tsv",
  show_col_types = FALSE
)

df <- df %>%
  distinct(ve_protein, pa_protein, .keep_all = TRUE) %>%
  mutate(
    combined_score = abs(ve_log2FC) + abs(pa_log2FC),
    gene_label = case_when(
      !is.na(ve_symbol) & ve_symbol != "" ~ paste0(ve_symbol, " | ", ve_gene),
      !is.na(pa_symbol) & pa_symbol != "" ~ paste0(pa_symbol, " | ", ve_gene),
      TRUE ~ ve_gene
    )
  )

top_df <- df %>%
  arrange(desc(combined_score)) %>%
  slice_head(n = 30)

write_tsv(
  top_df,
  "results/osmoadaptation/Osmoadaptation_Top30_for_Heatmap.tsv"
)

mat <- top_df %>%
  select(gene_label, ve_log2FC, pa_log2FC) %>%
  column_to_rownames("gene_label") %>%
  as.matrix()

colnames(mat) <- c("B. velezensis", "B. paralicheniformis")

annotation_row <- top_df %>%
  select(gene_label, response_class, functional_category) %>%
  column_to_rownames("gene_label") %>%
  as.data.frame()

png(
  "figures/osmoadaptation/Osmoadaptation_Heatmap_Top30.png",
  width = 2500,
  height = 3200,
  res = 300
)

pheatmap(
  mat,
  scale = "none",
  annotation_row = annotation_row,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 12,
  main = "Conserved osmoadaptation-associated orthologs",
  border_color = NA
)

dev.off()

pdf(
  "figures/osmoadaptation/Osmoadaptation_Heatmap_Top30.pdf",
  width = 9,
  height = 11
)

pheatmap(
  mat,
  scale = "none",
  annotation_row = annotation_row,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 8,
  fontsize_col = 12,
  main = "Conserved osmoadaptation-associated orthologs",
  border_color = NA
)

dev.off()

cat("\nGenerated:\n")
cat("figures/osmoadaptation/Osmoadaptation_Heatmap_Top30.png\n")
cat("figures/osmoadaptation/Osmoadaptation_Heatmap_Top30.pdf\n")
cat("results/osmoadaptation/Osmoadaptation_Top30_for_Heatmap.tsv\n")
