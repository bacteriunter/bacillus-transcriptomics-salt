#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

comparisons <- c(
  "velezensis_control_vs_salt",
  "paralicheniformis_control_vs_salt"
)

for (comp in comparisons) {

  message("Filtering DEGs for: ", comp)

  infile <- file.path("results/deseq2", comp, "deseq2_results.tsv")
  outdir <- file.path("results/deseq2", comp, "filtered")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  res <- read_tsv(infile, show_col_types = FALSE)

  sig <- res %>%
    filter(!is.na(padj)) %>%
    filter(padj < 0.05, abs(log2FoldChange) >= 1)

  up <- sig %>%
    filter(log2FoldChange >= 1)

  down <- sig %>%
    filter(log2FoldChange <= -1)

  write_tsv(sig, file.path(outdir, "DEGs_padj0.05_log2FC1.tsv"))
  write_tsv(up, file.path(outdir, "upregulated_padj0.05_log2FC1.tsv"))
  write_tsv(down, file.path(outdir, "downregulated_padj0.05_log2FC1.tsv"))

  summary <- tibble(
    comparison_group = comp,
    significant_abs_log2FC_1 = nrow(sig),
    upregulated = nrow(up),
    downregulated = nrow(down)
  )

  write_tsv(summary, file.path(outdir, "filtered_deg_summary.tsv"))
}

message("Filtering completed.")
