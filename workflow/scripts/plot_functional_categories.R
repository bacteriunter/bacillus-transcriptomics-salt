#!/usr/bin/env Rscript

library(tidyverse)

input_file <- "results/tables/functional_category_response_matrix.tsv"

output_png <- "figures/functional_categories_conserved_response.png"
output_pdf <- "figures/functional_categories_conserved_response.pdf"

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

df <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

df$response_class <- factor(
  df$response_class,
  levels = c(
    "conserved_up",
    "conserved_down",
    "opposite"
  )
)

df$functional_category <- factor(
  df$functional_category,
  levels = rev(
    c(
      "transport",
      "regulation",
      "glycine_betaine_choline",
      "sodium_homeostasis",
      "proline",
      "oxidative_stress",
      "biofilm_cell_envelope",
      "poly_gamma_glutamate"
    )
  )
)

p <- ggplot(
  df,
  aes(
    x = functional_category,
    y = n,
    fill = response_class
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      conserved_up = "#4DAF4A",
      conserved_down = "#377EB8",
      opposite = "#E41A1C"
    )
  ) +
  labs(
    title = "Conserved ortholog responses across functional categories",
    x = "Functional category",
    y = "Number of orthologs",
    fill = "Response class"
  ) +
  theme_bw(base_size = 18) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.position = "right"
  )

ggsave(
  filename = output_png,
  plot = p,
  width = 10,
  height = 6,
  dpi = 600
)

ggsave(
  filename = output_pdf,
  plot = p,
  width = 10,
  height = 6
)

cat("\nFigure generated:\n")
cat(output_png, "\n")
cat(output_pdf, "\n")
