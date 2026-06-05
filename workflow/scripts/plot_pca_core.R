#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(tidyverse)
})

samples <- read_tsv("metadata/sample_annotation/samples.tsv", show_col_types = FALSE)

comparisons <- c(
  "velezensis_control_vs_salt",
  "paralicheniformis_control_vs_salt"
)

dir.create("figures/pca", recursive = TRUE, showWarnings = FALSE)
dir.create("results/pca", recursive = TRUE, showWarnings = FALSE)

for (comp in comparisons) {

  message("PCA for: ", comp)

  meta <- samples %>%
    filter(comparison_group == comp) %>%
    mutate(
      condition = factor(condition, levels = c("control", "salt")),
      quant_file = file.path("results/salmon_quant", run_accession, "quant.sf")
    )

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
  vsd <- vst(dds, blind = FALSE)

  pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
  percentVar <- round(100 * attr(pca, "percentVar"), 1)

  pca <- pca %>%
    rownames_to_column("sample_id")

  write_tsv(
    pca,
    file.path("results/pca", paste0(comp, "_pca_coordinates.tsv"))
  )

  p <- ggplot(pca, aes(PC1, PC2, color = condition)) +
    geom_point(size = 4) +
    geom_text(
      aes(label = sample_id),
      vjust = -1,
      size = 3,
      show.legend = FALSE
    ) +
    labs(
      title = comp,
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      color = "Condition"
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path("figures/pca", paste0(comp, "_PCA.pdf")),
    p,
    width = 5.5,
    height = 4.5
  )

  ggsave(
    file.path("figures/pca", paste0(comp, "_PCA.png")),
    p,
    width = 5.5,
    height = 4.5,
    dpi = 300
  )
}

message("PCA plots completed.")
