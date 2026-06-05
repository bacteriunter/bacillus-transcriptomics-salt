#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

dir.create("figures/kegg_enrichment", recursive = TRUE, showWarnings = FALSE)

files <- c(
  conserved_up = "results/enrichment/kegg_ko/conserved_up_kegg_ko_enrichment.tsv",
  conserved_down = "results/enrichment/kegg_ko/conserved_down_kegg_ko_enrichment.tsv",
  opposite = "results/enrichment/kegg_ko/opposite_kegg_ko_enrichment.tsv"
)

df <- purrr::imap_dfr(
  files,
  ~ read_tsv(.x, show_col_types = FALSE) %>%
    mutate(response_set = .y)
)

plot_df <- df %>%
  group_by(response_set) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(
    label = paste0(KO, " | ", stringr::str_trunc(definition, 45)),
    neg_log10_pvalue = -log10(pvalue),
    significant_fdr = if_else(padj < 0.05, "FDR < 0.05", "FDR ≥ 0.05")
  )

write_tsv(
  plot_df,
  "results/enrichment/kegg_ko/top10_nominal_kegg_ko_enrichment.tsv"
)

p <- ggplot(
  plot_df,
  aes(
    x = reorder(label, neg_log10_pvalue),
    y = neg_log10_pvalue,
    fill = significant_fdr
  )
) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ response_set, scales = "free_y") +
  labs(
    title = "Top nominal KEGG KO enrichment results",
    x = "KO",
    y = "-log10 nominal p-value",
    fill = "Adjusted significance"
  ) +
  theme_bw(base_size = 12)

ggsave(
  "figures/kegg_enrichment/top10_nominal_kegg_ko_enrichment.png",
  p,
  width = 12,
  height = 8,
  dpi = 600
)

ggsave(
  "figures/kegg_enrichment/top10_nominal_kegg_ko_enrichment.pdf",
  p,
  width = 12,
  height = 8
)

cat("\nGenerated KEGG KO enrichment figure.\n")
