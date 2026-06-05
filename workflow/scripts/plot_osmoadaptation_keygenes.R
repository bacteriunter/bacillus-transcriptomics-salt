#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

dir.create("figures/osmoadaptation", recursive = TRUE, showWarnings = FALSE)
dir.create("results/osmoadaptation", recursive = TRUE, showWarnings = FALSE)

df <- read_tsv(
  "results/tables/Table3_Osmoadaptation_Candidates.tsv",
  show_col_types = FALSE
)

key_genes <- c(
  "proB", "proC", "putP",
  "gbsB", "betB", "opuAA", "opuAB",
  "nhaC", "pgsB", "pgsC"
)

df_key <- df %>%
  distinct(ve_protein, pa_protein, .keep_all = TRUE) %>%
  filter(
    ve_symbol %in% key_genes |
    pa_symbol %in% key_genes
  ) %>%
  mutate(
    gene_label = case_when(
      !is.na(ve_symbol) & ve_symbol != "" ~ ve_symbol,
      !is.na(pa_symbol) & pa_symbol != "" ~ pa_symbol,
      TRUE ~ ve_gene
    ),
    module = case_when(
      gene_label %in% c("proB", "proC", "putP") ~ "Proline",
      gene_label %in% c("gbsB", "betB", "opuAA", "opuAB") ~ "Glycine betaine / choline",
      gene_label %in% c("nhaC") ~ "Sodium homeostasis",
      gene_label %in% c("pgsB", "pgsC") ~ "Poly-gamma-glutamate",
      TRUE ~ "Other"
    )
  )

plot_df <- df_key %>%
  select(gene_label, module, ve_log2FC, pa_log2FC) %>%
  pivot_longer(
    cols = c(ve_log2FC, pa_log2FC),
    names_to = "species",
    values_to = "log2FC"
  ) %>%
  mutate(
    species = recode(
      species,
      ve_log2FC = "B. velezensis",
      pa_log2FC = "B. paralicheniformis"
    )
  )

write_tsv(
  plot_df,
  "results/osmoadaptation/Osmoadaptation_KeyGenes_Barplot_Data.tsv"
)

p <- ggplot(
  plot_df,
  aes(
    x = gene_label,
    y = log2FC,
    fill = species
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  facet_wrap(
    ~ module,
    scales = "free_x"
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  labs(
    title = "Key conserved osmoadaptation genes under salt stress",
    x = "Gene",
    y = "log2 fold change (salt / control)",
    fill = "Species"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

ggsave(
  "figures/osmoadaptation/Osmoadaptation_KeyGenes_Barplot.png",
  p,
  width = 11,
  height = 6,
  dpi = 600
)

ggsave(
  "figures/osmoadaptation/Osmoadaptation_KeyGenes_Barplot.pdf",
  p,
  width = 11,
  height = 6
)

cat("\nGenerated:\n")
cat("figures/osmoadaptation/Osmoadaptation_KeyGenes_Barplot.png\n")
cat("figures/osmoadaptation/Osmoadaptation_KeyGenes_Barplot.pdf\n")
cat("results/osmoadaptation/Osmoadaptation_KeyGenes_Barplot_Data.tsv\n")
