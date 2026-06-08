#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

dir.create("figures/response_summary", recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript/figures", recursive = TRUE, showWarnings = FALSE)

infile <- "results/orthology/response_classes/response_class_summary.tsv"

df <- read_tsv(infile, show_col_types = FALSE)

df <- df %>%
  transmute(
    response_class_raw = response_class,
    n_orthologs = as.integer(n_orthologs)
  ) %>%
  mutate(
    response_class = case_when(
      response_class_raw == "conserved_up" ~ "Conserved up",
      response_class_raw == "conserved_down" ~ "Conserved down",
      response_class_raw == "opposite" ~ "Opposite",
      TRUE ~ response_class_raw
    ),
    response_class = factor(
      response_class,
      levels = c("Conserved up", "Conserved down", "Opposite")
    )
  ) %>%
  arrange(response_class) %>%
  mutate(
    percent = 100 * n_orthologs / sum(n_orthologs),
    label = paste0(n_orthologs, " (", sprintf("%.1f", percent), "%)")
  )

p <- ggplot(
  df,
  aes(x = response_class, y = n_orthologs, fill = response_class)
) +
  geom_col(width = 0.68, show.legend = FALSE) +
  geom_text(
    aes(label = label),
    vjust = -0.45,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "Conserved up" = "#4DAF4A",
      "Conserved down" = "#377EB8",
      "Opposite" = "#E41A1C"
    )
  ) +
  scale_y_continuous(
    limits = c(0, max(df$n_orthologs) * 1.18),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Shared orthologous DEG response classes",
    x = "Ortholog response class",
    y = "Number of orthologs"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 12, 10, 12)
  )

ggsave(
  "figures/response_summary/response_class_summary_barplot_publication.png",
  p,
  width = 8,
  height = 5.5,
  dpi = 600
)

ggsave(
  "figures/response_summary/response_class_summary_barplot_publication.pdf",
  p,
  width = 8,
  height = 5.5
)

ggsave(
  "manuscript/figures/Figure2_ResponseClassSummary.png",
  p,
  width = 8,
  height = 5.5,
  dpi = 600
)

ggsave(
  "manuscript/figures/Figure2_ResponseClassSummary.pdf",
  p,
  width = 8,
  height = 5.5
)

cat("Saved:\n")
cat("figures/response_summary/response_class_summary_barplot_publication.png\n")
cat("figures/response_summary/response_class_summary_barplot_publication.pdf\n")
cat("manuscript/figures/Figure2_ResponseClassSummary.png\n")
cat("manuscript/figures/Figure2_ResponseClassSummary.pdf\n")
