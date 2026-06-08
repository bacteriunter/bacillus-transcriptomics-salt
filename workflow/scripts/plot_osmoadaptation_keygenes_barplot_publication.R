#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

dir.create("figures/osmoadaptation", recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript/figures", recursive = TRUE, showWarnings = FALSE)

infile <- "results/tables/Table2_Core_Osmoadaptation_Genes.tsv"

df <- read_tsv(infile, show_col_types = FALSE)

df <- df %>%
  transmute(
    gene_label = Gene,
    functional_category_clean = case_when(
      str_detect(Category, "glycine|betaine|choline") ~ "Glycine betaine / choline",
      str_detect(Category, "poly|gamma|glutamate") ~ "Poly-γ-glutamate",
      str_detect(Category, "proline") ~ "Proline",
      str_detect(Category, "sodium|homeostasis") ~ "Sodium homeostasis",
      str_detect(Category, "transport") ~ "Transport",
      str_detect(Category, "oxidative") ~ "Oxidative stress",
      str_detect(Category, "biofilm|cell|envelope") ~ "Biofilm / cell envelope",
      str_detect(Category, "regulation") ~ "Regulation",
      TRUE ~ Category
    ),
    `B. velezensis` = as.numeric(log2FC_B_velezensis),
    `B. paralicheniformis` = as.numeric(log2FC_B_paralicheniformis)
  ) %>%
  mutate(
    gene_label = gsub("_", " ", gene_label)
  )

plot_df <- df %>%
  pivot_longer(
    cols = c(`B. velezensis`, `B. paralicheniformis`),
    names_to = "species",
    values_to = "log2FC"
  ) %>%
  filter(!is.na(log2FC)) %>%
  mutate(
    species = factor(
      species,
      levels = c("B. paralicheniformis", "B. velezensis")
    ),
    functional_category_clean = factor(
      functional_category_clean,
      levels = c(
        "Glycine betaine / choline",
        "Poly-γ-glutamate",
        "Proline",
        "Sodium homeostasis",
        "Transport",
        "Oxidative stress",
        "Biofilm / cell envelope",
        "Regulation"
      )
    )
  )

p <- ggplot(
  plot_df,
  aes(x = gene_label, y = log2FC, fill = species)
) +
  geom_hline(yintercept = 0, linewidth = 0.45, color = "black") +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  facet_wrap(
    ~ functional_category_clean,
    scales = "free_x",
    ncol = 2
  ) +
  scale_fill_manual(
    values = c(
      "B. paralicheniformis" = "#F8766D",
      "B. velezensis" = "#00BFC4"
    ),
    labels = c(
      expression(italic("B. paralicheniformis")),
      expression(italic("B. velezensis"))
    )
  ) +
  labs(
    title = "Key conserved osmoadaptation genes under salt stress",
    x = "Gene",
    y = expression(log[2]~"FC (Salt / Control)"),
    fill = "Species"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = "grey85", color = "grey30"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 12, 10, 12)
  )

ggsave("figures/osmoadaptation/Osmoadaptation_KeyGenes_Barplot.png", p, width = 11, height = 7, dpi = 600)
ggsave("figures/osmoadaptation/Osmoadaptation_KeyGenes_Barplot.pdf", p, width = 11, height = 7)

ggsave("manuscript/figures/Figure6_KeyOsmoadaptationGenes.png", p, width = 11, height = 7, dpi = 600)
ggsave("manuscript/figures/Figure6_KeyOsmoadaptationGenes.pdf", p, width = 11, height = 7)

cat("\nFigure 6 generated:\n")
cat("figures/osmoadaptation/Osmoadaptation_KeyGenes_Barplot.png\n")
cat("figures/osmoadaptation/Osmoadaptation_KeyGenes_Barplot.pdf\n")
cat("manuscript/figures/Figure6_KeyOsmoadaptationGenes.png\n")
cat("manuscript/figures/Figure6_KeyOsmoadaptationGenes.pdf\n")
