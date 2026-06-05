#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(tidyverse)
})

samples <- read_tsv("metadata/sample_annotation/samples.tsv", show_col_types = FALSE)

core_comparisons <- c(
  "velezensis_control_vs_salt",
  "paralicheniformis_control_vs_salt"
)

dir.create("results/deseq2", recursive = TRUE, showWarnings = FALSE)
dir.create("results/counts", recursive = TRUE, showWarnings = FALSE)

for (comp in core_comparisons) {

  message("Running DESeq2 for: ", comp)

  meta <- samples %>%
    filter(comparison_group == comp) %>%
    mutate(
      condition = factor(condition, levels = c("control", "salt")),
      quant_file = file.path("results/salmon_quant", run_accession, "quant.sf")
    )

  stopifnot(nrow(meta) == 6)
  stopifnot(all(table(meta$condition) == 3))
  stopifnot(all(file.exists(meta$quant_file)))

  files <- meta$quant_file
  names(files) <- meta$sample_id

  txi <- tximport(
    files,
    type = "salmon",
    txOut = TRUE,
    countsFromAbundance = "no"
  )

  dds <- DESeqDataSetFromTximport(
    txi,
    colData = as.data.frame(meta),
    design = ~ condition
  )

  keep <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep, ]

  dds <- DESeq(dds)

  res <- results(dds, contrast = c("condition", "salt", "control")) %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    arrange(padj)

  outdir <- file.path("results/deseq2", comp)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  write_tsv(res, file.path(outdir, "deseq2_results.tsv"))

  norm_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("gene_id")

  raw_counts <- counts(dds, normalized = FALSE) %>%
    as.data.frame() %>%
    rownames_to_column("gene_id")

  write_tsv(norm_counts, file.path(outdir, "normalized_counts.tsv"))
  write_tsv(raw_counts, file.path(outdir, "raw_counts.tsv"))

  summary_tbl <- tibble(
    comparison_group = comp,
    n_genes_tested = nrow(res),
    significant_padj_0_05 = sum(res$padj < 0.05, na.rm = TRUE),
    up_padj_0_05_log2FC_1 = sum(res$padj < 0.05 & res$log2FoldChange >= 1, na.rm = TRUE),
    down_padj_0_05_log2FC_minus1 = sum(res$padj < 0.05 & res$log2FoldChange <= -1, na.rm = TRUE)
  )

  write_tsv(summary_tbl, file.path(outdir, "deseq2_summary.tsv"))
}

message("DESeq2 core analyses completed.")
