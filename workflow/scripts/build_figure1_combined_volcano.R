#!/usr/bin/env Rscript

library(grid)

dir.create(
  "figures/manuscript",
  recursive = TRUE,
  showWarnings = FALSE
)

img1 <- rasterImage <- png::readPNG
