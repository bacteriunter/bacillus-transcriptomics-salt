#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

dir.create("results/correlation", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/correlation", recursive = TRUE, showWarnings = FALSE)

infile <- "results/tables/master_conserved_ortholog_deg_table.tsv"

df <- read_tsv(infile, show_col_types = FALSE) %>%
  mutate(
    ve_log2FC = as.numeric(ve_log2FC),
    pa_log2FC = as.numeric(pa_log2FC),
    response_class = factor(
      response_class,
      levels = c("conserved_up", "conserved_down", "opposite")
    ),
    abs_sum = abs(ve_log2FC) + abs(pa_log2FC),
    gene_label = case_when(
      !is.na(ve_symbol) & ve_symbol != "" ~ ve_symbol,
      !is.na(pa_symbol) & pa_symbol != "" ~ pa_symbol,
      TRUE ~ ve_protein
    )
  ) %>%
  filter(
    !is.na(ve_log2FC),
    !is.na(pa_log2FC),
    !is.na(response_class)
  ) %>%
  arrange(ve_protein, pa_protein, desc(abs_sum)) %>%
  distinct(ve_protein, pa_protein, .keep_all = TRUE)

pearson <- cor.test(df$ve_log2FC, df$pa_log2FC, method = "pearson")
spearman <- cor.test(df$ve_log2FC, df$pa_log2FC, method = "spearman", exact = FALSE)

summary_df <- tibble(
  n_orthologs = nrow(df),
  pearson_r = unname(pearson$estimate),
  pearson_pvalue = pearson$p.value,
  spearman_rho = unname(spearman$estimate),
  spearman_pvalue = spearman$p.value,
  conserved_up_n = sum(df$response_class == "conserved_up"),
  conserved_down_n = sum(df$response_class == "conserved_down"),
  opposite_n = sum(df$response_class == "opposite")
)

write_tsv(summary_df, "results/correlation/ortholog_log2FC_correlation_summary.tsv")
write_tsv(df, "results/correlation/ortholog_log2FC_scatter_data.tsv")

label_df <- df %>%
  arrange(desc(abs_sum)) %>%
  slice_head(n = 15)

p <- ggplot(
  df,
  aes(
    x = ve_log2FC,
    y = pa_log2FC,
    color = response_class
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    linetype = "dashed",
    color = "grey45"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    linetype = "dashed",
    color = "grey45"
  ) +
  geom_point(alpha = 0.75, size = 2.1) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.6,
    color = "black"
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = gene_label),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.linewidth = 0.25,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = min(df$ve_log2FC, na.rm = TRUE),
    y = max(df$pa_log2FC, na.rm = TRUE),
    hjust = 0,
    vjust = 1,
    size = 4,
    label = paste0(
      "n = ", nrow(df),
      "\nPearson r = ", round(unname(pearson$estimate), 3),
      ", p = ", formatC(pearson$p.value, format = "e", digits = 2),
      "\nSpearman rho = ", round(unname(spearman$estimate), 3),
      ", p = ", formatC(spearman$p.value, format = "e", digits = 2)
    )
  ) +
  scale_color_manual(
    values = c(
      conserved_up = "#D55E00",
      conserved_down = "#0072B2",
      opposite = "#009E73"
    ),
    labels = c(
      conserved_up = "Conserved up",
      conserved_down = "Conserved down",
      opposite = "Opposite"
    )
  ) +
  labs(
    x = expression(log[2]~FC~italic("B. velezensis")),
    y = expression(log[2]~FC~italic("B. paralicheniformis")),
    color = "Response class",
    title = "Ortholog-level transcriptional conservation under salt stress"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/correlation/Figure_Ortholog_log2FC_Correlation.png",
  p,
  width = 8.5,
  height = 6.5,
  dpi = 600
)

ggsave(
  "figures/correlation/Figure_Ortholog_log2FC_Correlation.pdf",
  p,
  width = 8.5,
  height = 6.5
)

cat("\nGenerated corrected deduplicated correlation analysis:\n")
print(summary_df)
