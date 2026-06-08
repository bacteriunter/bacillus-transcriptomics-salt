#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

dir.create("figures/osmoadaptation", recursive = TRUE, showWarnings = FALSE)
dir.create("results/osmoadaptation", recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript/figures", recursive = TRUE, showWarnings = FALSE)

infile <- "results/tables/Table3_Osmoadaptation_Candidates.tsv"

df <- read_tsv(infile, show_col_types = FALSE)

df <- df %>%
  mutate(
    raw_label = case_when(
      !is.na(ve_symbol) & ve_symbol != "" ~ ve_symbol,
      !is.na(pa_symbol) & pa_symbol != "" ~ pa_symbol,
      TRUE ~ ve_gene
    ),
    label = gsub("_", " ", raw_label),
    label = make.unique(label),
    combined_abs_log2FC = abs(as.numeric(ve_log2FC)) + abs(as.numeric(pa_log2FC))
  ) %>%
  arrange(desc(combined_abs_log2FC)) %>%
  slice_head(n = 30)

mat <- df %>%
  select(label, ve_log2FC, pa_log2FC) %>%
  mutate(
    ve_log2FC = as.numeric(ve_log2FC),
    pa_log2FC = as.numeric(pa_log2FC)
  ) %>%
  column_to_rownames("label") %>%
  as.matrix()

colnames(mat) <- c("Bvelezensis", "Bparalicheniformis")

ann_row <- df %>%
  select(label, response_class, functional_category) %>%
  mutate(
    `Response class` = case_when(
      response_class == "conserved_up" ~ "Conserved up",
      response_class == "conserved_down" ~ "Conserved down",
      response_class == "opposite" ~ "Opposite",
      TRUE ~ response_class
    ),
    `Functional category` = case_when(
      str_detect(functional_category, "glycine_betaine_choline") ~ "Glycine betaine / choline",
      str_detect(functional_category, "sodium_homeostasis") ~ "Sodium homeostasis",
      str_detect(functional_category, "proline") ~ "Proline",
      str_detect(functional_category, "oxidative_stress") ~ "Oxidative stress",
      str_detect(functional_category, "poly_gamma_glutamate") ~ "Poly-γ-glutamate",
      str_detect(functional_category, "biofilm_cell_envelope") ~ "Biofilm / cell envelope",
      str_detect(functional_category, "regulation") ~ "Regulation",
      str_detect(functional_category, "transport") ~ "Transport",
      TRUE ~ "Unclassified"
    )
  ) %>%
  select(label, `Functional category`, `Response class`) %>%
  column_to_rownames("label")

ann_colors <- list(
  `Response class` = c(
    "Conserved up" = "#4DAF4A",
    "Conserved down" = "#377EB8",
    "Opposite" = "#E41A1C"
  ),
  `Functional category` = c(
    "Glycine betaine / choline" = "#CC79A7",
    "Sodium homeostasis" = "#009E73",
    "Proline" = "#F781BF",
    "Oxidative stress" = "#A6AD00",
    "Poly-γ-glutamate" = "#CAB2D6",
    "Biofilm / cell envelope" = "#FB9A99",
    "Regulation" = "#56B4E9",
    "Transport" = "#E69F00",
    "Unclassified" = "#BDBDBD"
  )
)

write_tsv(
  df,
  "results/osmoadaptation/osmoadaptation_top30_selected.tsv"
)

write_tsv(
  as.data.frame(mat) %>% rownames_to_column("gene_label"),
  "results/osmoadaptation/osmoadaptation_top30_heatmap_matrix.tsv"
)

col_labels <- as.expression(c(
  bquote(italic("B. velezensis")),
  bquote(italic("B. paralicheniformis"))
))

draw_heatmap <- function(filename, type = c("png", "pdf")) {
  type <- match.arg(type)

  if (type == "png") {
    png(filename, width = 3000, height = 3800, res = 300)
  } else {
    pdf(filename, width = 8.8, height = 10.5)
  }

  pheatmap(
    mat,
    annotation_row = ann_row,
    annotation_colors = ann_colors,
    labels_col = col_labels,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "none",
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 13,
    main = "Conserved osmoadaptation-associated orthologs",
    border_color = NA,
    angle_col = "0",
    annotation_names_row = FALSE,
    annotation_names_col = FALSE
  )

  dev.off()
}

draw_heatmap(
  "figures/osmoadaptation/Osmoadaptation_Heatmap_Top30.png",
  "png"
)

draw_heatmap(
  "figures/osmoadaptation/Osmoadaptation_Heatmap_Top30.pdf",
  "pdf"
)

file.copy(
  "figures/osmoadaptation/Osmoadaptation_Heatmap_Top30.png",
  "manuscript/figures/Figure5_OsmoadaptationHeatmapTop30.png",
  overwrite = TRUE
)

file.copy(
  "figures/osmoadaptation/Osmoadaptation_Heatmap_Top30.pdf",
  "manuscript/figures/Figure5_OsmoadaptationHeatmapTop30.pdf",
  overwrite = TRUE
)

cat("\nHeatmap generated:\n")
cat("figures/osmoadaptation/Osmoadaptation_Heatmap_Top30.png\n")
cat("figures/osmoadaptation/Osmoadaptation_Heatmap_Top30.pdf\n")
cat("manuscript/figures/Figure5_OsmoadaptationHeatmapTop30.png\n")
cat("manuscript/figures/Figure5_OsmoadaptationHeatmapTop30.pdf\n")
