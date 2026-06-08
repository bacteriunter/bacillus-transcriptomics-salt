#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("manuscript/figures", recursive = TRUE, showWarnings = FALSE)

input_file <- "results/tables/functional_category_response_matrix.tsv"

df <- read_tsv(input_file, show_col_types = FALSE)

df <- df %>%
  mutate(
    response_class = case_when(
      response_class == "conserved_up" ~ "Conserved up",
      response_class == "conserved_down" ~ "Conserved down",
      response_class == "opposite" ~ "Opposite",
      TRUE ~ response_class
    ),
    response_class = factor(
      response_class,
      levels = c("Conserved up", "Conserved down", "Opposite")
    ),
    functional_category = case_when(
      functional_category == "transport" ~ "Transport",
      functional_category == "regulation" ~ "Regulation",
      functional_category == "glycine_betaine_choline" ~ "Glycine betaine / choline",
      functional_category == "sodium_homeostasis" ~ "Sodium homeostasis",
      functional_category == "proline" ~ "Proline",
      functional_category == "oxidative_stress" ~ "Oxidative stress",
      functional_category == "biofilm_cell_envelope" ~ "Biofilm / cell envelope",
      functional_category == "poly_gamma_glutamate" ~ "Poly-γ-glutamate",
      functional_category == "unclassified" ~ "Unclassified",
      TRUE ~ str_replace_all(functional_category, "_", " ")
    )
  )

if ("n_orthologs" %in% names(df)) {
  count_col <- "n_orthologs"
} else if ("count" %in% names(df)) {
  count_col <- "count"
} else if ("n" %in% names(df)) {
  count_col <- "n"
} else {
  stop("No count column found. Expected n_orthologs, count, or n.")
}

df <- df %>%
  rename(n_orthologs = all_of(count_col)) %>%
  group_by(functional_category) %>%
  mutate(total = sum(n_orthologs)) %>%
  ungroup() %>%
  mutate(
    functional_category = fct_reorder(functional_category, total)
  )

p <- ggplot(
  df,
  aes(
    x = n_orthologs,
    y = functional_category,
    fill = response_class
  )
) +
  geom_col(width = 0.72) +
  scale_fill_manual(
    values = c(
      "Conserved up" = "#4DAF4A",
      "Conserved down" = "#377EB8",
      "Opposite" = "#E41A1C"
    )
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Functional categories of orthologous response classes",
    x = "Number of orthologs",
    y = "Functional category",
    fill = "Response class"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 12, 10, 12)
  )

ggsave(
  "figures/functional_categories_conserved_response_publication.png",
  p,
  width = 8.5,
  height = 6,
  dpi = 600
)

ggsave(
  "figures/functional_categories_conserved_response_publication.pdf",
  p,
  width = 8.5,
  height = 6
)

ggsave(
  "manuscript/figures/Figure3_FunctionalCategories.png",
  p,
  width = 8.5,
  height = 6,
  dpi = 600
)

ggsave(
  "manuscript/figures/Figure3_FunctionalCategories.pdf",
  p,
  width = 8.5,
  height = 6
)

cat("Saved:\n")
cat("figures/functional_categories_conserved_response_publication.png\n")
cat("figures/functional_categories_conserved_response_publication.pdf\n")
cat("manuscript/figures/Figure3_FunctionalCategories.png\n")
cat("manuscript/figures/Figure3_FunctionalCategories.pdf\n")
