#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

dir.create("figures/heatmap", recursive = TRUE, showWarnings = FALSE)
dir.create("results/heatmap", recursive = TRUE, showWarnings = FALSE)

df <- read_tsv(
  "results/tables/master_conserved_ortholog_deg_table.tsv",
  show_col_types = FALSE
)

# Usar una fila por ortólogo; si hay múltiples KO por proteína, conservar la primera
df_unique <- df %>%
  distinct(ve_protein, .keep_all = TRUE) %>%
  mutate(
    combined_abs_log2FC = abs(ve_log2FC) + abs(pa_log2FC),
    label = case_when(
      !is.na(ve_symbol) & ve_symbol != "" ~ paste0(ve_symbol, " | ", ve_gene),
      !is.na(pa_symbol) & pa_symbol != "" ~ paste0(pa_symbol, " | ", ve_gene),
      TRUE ~ ve_gene
    )
  )

top50 <- df_unique %>%
  filter(response_class %in% c("conserved_up", "conserved_down")) %>%
  arrange(desc(combined_abs_log2FC)) %>%
  slice_head(n = 50)

mat <- top50 %>%
  select(label, ve_log2FC, pa_log2FC) %>%
  column_to_rownames("label") %>%
  as.matrix()

colnames(mat) <- c("B. velezensis", "B. paralicheniformis")

ann_row <- top50 %>%
  select(label, response_class, functional_category) %>%
  column_to_rownames("label")

write_tsv(
  top50,
  "results/heatmap/top50_conserved_orthologs_selected.tsv"
)

write_tsv(
  as.data.frame(mat) %>% rownames_to_column("gene_label"),
  "results/heatmap/top50_conserved_orthologs_heatmap_matrix.tsv"
)

png(
  "figures/heatmap/top50_conserved_orthologs_heatmap.png",
  width = 2600,
  height = 3200,
  res = 300
)

pheatmap(
  mat,
  annotation_row = ann_row,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  fontsize_row = 7,
  fontsize_col = 12,
  main = "Top 50 conserved salt-responsive orthologs",
  border_color = NA
)

dev.off()

pdf(
  "figures/heatmap/top50_conserved_orthologs_heatmap.pdf",
  width = 8,
  height = 10
)

pheatmap(
  mat,
  annotation_row = ann_row,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  fontsize_row = 7,
  fontsize_col = 12,
  main = "Top 50 conserved salt-responsive orthologs",
  border_color = NA
)

dev.off()

cat("\nHeatmap generated:\n")
cat("figures/heatmap/top50_conserved_orthologs_heatmap.png\n")
cat("figures/heatmap/top50_conserved_orthologs_heatmap.pdf\n")
