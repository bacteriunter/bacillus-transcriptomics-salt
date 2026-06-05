#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

dir.create("figures/response_summary", recursive = TRUE, showWarnings = FALSE)
dir.create("results/response_summary", recursive = TRUE, showWarnings = FALSE)

df <- read_tsv(
  "results/orthology/conserved_deg_orthologs.tsv",
  show_col_types = FALSE
)

summary <- df %>%
  count(response_class, name = "n_orthologs") %>%
  mutate(
    response_class = factor(
      response_class,
      levels = c("conserved_up", "conserved_down", "opposite")
    ),
    percent = 100 * n_orthologs / sum(n_orthologs),
    label = paste0(n_orthologs, " (", round(percent, 1), "%)")
  ) %>%
  arrange(response_class)

write_tsv(
  summary,
  "results/response_summary/response_class_summary.tsv"
)

p_bar <- ggplot(
  summary,
  aes(x = response_class, y = n_orthologs, fill = response_class)
) +
  geom_col(width = 0.7) +
  geom_text(
    aes(label = label),
    vjust = -0.4,
    size = 5,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      conserved_up = "#4DAF4A",
      conserved_down = "#377EB8",
      opposite = "#E41A1C"
    )
  ) +
  labs(
    title = "Shared orthologous DEG response classes",
    x = "Response class",
    y = "Number of orthologs"
  ) +
  theme_bw(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  ylim(0, max(summary$n_orthologs) * 1.15)

ggsave(
  "figures/response_summary/response_class_summary_barplot.png",
  p_bar,
  width = 7,
  height = 5,
  dpi = 600
)

ggsave(
  "figures/response_summary/response_class_summary_barplot.pdf",
  p_bar,
  width = 7,
  height = 5
)

p_stack <- ggplot(
  summary,
  aes(x = "Shared orthologous DEGs", y = n_orthologs, fill = response_class)
) +
  geom_col(width = 0.5) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 5,
    color = "white",
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      conserved_up = "#4DAF4A",
      conserved_down = "#377EB8",
      opposite = "#E41A1C"
    )
  ) +
  labs(
    title = "Composition of shared orthologous DEGs",
    x = NULL,
    y = "Number of orthologs",
    fill = "Response class"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  "figures/response_summary/response_class_summary_stacked.png",
  p_stack,
  width = 6,
  height = 6,
  dpi = 600
)

ggsave(
  "figures/response_summary/response_class_summary_stacked.pdf",
  p_stack,
  width = 6,
  height = 6
)

cat("\nGenerated response class summary figures:\n")
cat("figures/response_summary/response_class_summary_barplot.png\n")
cat("figures/response_summary/response_class_summary_stacked.png\n")
