#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

comparisons <- c(
  "velezensis_control_vs_salt",
  "paralicheniformis_control_vs_salt"
)

dir.create("figures/volcano", recursive = TRUE, showWarnings = FALSE)
dir.create("results/volcano", recursive = TRUE, showWarnings = FALSE)

label_genes <- c(
  "proB", "proC", "putP",
  "betB", "gbsB", "opuAA", "opuAB",
  "pgsB", "pgsC", "nhaC"
)

for (comp in comparisons) {

  message("Volcano for: ", comp)

  res <- read_tsv(
    file.path("results/deseq2", comp, "deseq2_results.tsv"),
    show_col_types = FALSE
  )

  if (grepl("velezensis", comp)) {
    annot <- read_tsv("results/annotation/velezensis_annotation.tsv", show_col_types = FALSE)
  } else {
    annot <- read_tsv("results/annotation/paralicheniformis_annotation.tsv", show_col_types = FALSE)
  }

  plot_df <- res %>%
    left_join(annot, by = "gene_id") %>%
    mutate(
      neg_log10_padj = -log10(padj),
      deg_class = case_when(
        padj < 0.05 & log2FoldChange >= 1 ~ "Up",
        padj < 0.05 & log2FoldChange <= -1 ~ "Down",
        TRUE ~ "Not significant"
      ),
      label = if_else(gene_symbol %in% label_genes, gene_symbol, NA_character_)
    )

  write_tsv(
    plot_df,
    file.path("results/volcano", paste0(comp, "_volcano_table.tsv"))
  )

  p <- ggplot(plot_df, aes(log2FoldChange, neg_log10_padj, color = deg_class)) +
    geom_point(alpha = 0.65, size = 1.4, na.rm = TRUE) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(
      aes(label = label),
      max.overlaps = 50,
      size = 3,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    labs(
      title = comp,
      x = "log2 fold change (salt / control)",
      y = "-log10 adjusted p-value",
      color = "DEG class"
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path("figures/volcano", paste0(comp, "_volcano.pdf")),
    p,
    width = 6,
    height = 5
  )

  ggsave(
    file.path("figures/volcano", paste0(comp, "_volcano.png")),
    p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

message("Volcano plots completed.")
