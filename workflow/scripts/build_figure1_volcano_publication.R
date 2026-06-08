#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(grid)
  library(gtable)
})

dir.create("figures/manuscript", recursive = TRUE, showWarnings = FALSE)

threshold_padj <- 0.05
threshold_lfc  <- 1

x_limits <- c(-15, 15)
y_limits <- c(0, 310)

files <- list(
  velezensis = "results/volcano/velezensis_control_vs_salt_volcano_table.tsv",
  paralicheniformis = "results/volcano/paralicheniformis_control_vs_salt_volcano_table.tsv"
)

panel_titles <- list(
  velezensis = expression(italic("B. velezensis")~"(Control vs Salt)"),
  paralicheniformis = expression(italic("B. paralicheniformis")~"(Control vs Salt)")
)

read_volcano <- function(path, species_key) {
  df <- read_tsv(path, show_col_types = FALSE)

  gene_col <- intersect(c("gene", "gene_id", "id", "Gene", "GeneID"), names(df))[1]
  lfc_col  <- intersect(c("log2FoldChange", "log2FC", "log2_fold_change", "logFC"), names(df))[1]
  padj_col <- intersect(c("padj", "adjusted_pvalue", "adjusted_p_value", "FDR"), names(df))[1]
  label_col <- intersect(c("symbol", "gene_symbol", "label", "GeneSymbol"), names(df))[1]

  if (is.na(gene_col) || is.na(lfc_col) || is.na(padj_col)) {
    stop(
      "Required columns not found in: ", path, "\n",
      "Available columns: ", paste(names(df), collapse = ", ")
    )
  }

  df <- df %>%
    mutate(
      gene_id = .data[[gene_col]],
      log2FC = as.numeric(.data[[lfc_col]]),
      padj = as.numeric(.data[[padj_col]]),
      neg_log10_padj = -log10(pmax(padj, .Machine$double.xmin)),
      neg_log10_padj_plot = pmin(neg_log10_padj, y_limits[2]),
      expression_class = case_when(
        !is.na(padj) & !is.na(log2FC) & padj < threshold_padj & log2FC >= threshold_lfc ~ "Upregulated",
        !is.na(padj) & !is.na(log2FC) & padj < threshold_padj & log2FC <= -threshold_lfc ~ "Downregulated",
        TRUE ~ "Not significant"
      ),
      species = species_key
    )

  if (!is.na(label_col)) {
    df <- df %>% mutate(label = .data[[label_col]])
  } else {
    df <- df %>% mutate(label = gene_id)
  }

  df
}

df <- bind_rows(
  read_volcano(files$velezensis, "velezensis"),
  read_volcano(files$paralicheniformis, "paralicheniformis")
)

target_genes <- c(
  "betB", "gbsB", "opuAA", "opuAB", "opuA",
  "proB", "proC", "putP", "nhaC", "pgsB", "pgsC"
)

label_df <- df %>%
  filter(label %in% target_genes) %>%
  group_by(species, label) %>%
  slice_max(order_by = abs(log2FC), n = 1, with_ties = FALSE) %>%
  ungroup()

make_panel <- function(species_key, keep_legend = FALSE) {
  d <- df %>% filter(species == species_key)
  lab <- label_df %>% filter(species == species_key)

  ggplot(d, aes(x = log2FC, y = neg_log10_padj_plot, color = expression_class)) +
    geom_point(alpha = 0.65, size = 1.75) +
    geom_vline(
      xintercept = c(-threshold_lfc, threshold_lfc),
      linetype = "dashed",
      linewidth = 0.45
    ) +
    geom_hline(
      yintercept = -log10(threshold_padj),
      linetype = "dashed",
      linewidth = 0.45
    ) +
    geom_text_repel(
      data = lab,
      aes(label = label),
      size = 3.2,
      max.overlaps = Inf,
      box.padding = 0.55,
      point.padding = 0.35,
      force = 2,
      force_pull = 0.4,
      min.segment.length = 0,
      segment.size = 0.25,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Downregulated" = "#F8766D",
        "Not significant" = "#00BA38",
        "Upregulated" = "#619CFF"
      ),
      breaks = c("Downregulated", "Not significant", "Upregulated")
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = seq(-15, 15, 5),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    coord_cartesian(
      ylim = y_limits,
      clip = "off"
    ) +
    labs(
      title = panel_titles[[species_key]],
      x = expression(log[2]~"FC (Salt / Control)"),
      y = expression(-log[10]~"(adjusted p-value)"),
      color = "Expression response"
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      panel.grid.minor = element_line(linewidth = 0.2),
      panel.grid.major = element_line(linewidth = 0.35),
      legend.position = if (keep_legend) "right" else "none",
      plot.margin = margin(10, 10, 10, 10)
    )
}

get_legend <- function(plot) {
  g <- ggplotGrob(plot)
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  if (length(idx) == 0) return(nullGrob())
  g$grobs[[idx]]
}

p_left <- make_panel("velezensis", keep_legend = FALSE)
p_right <- make_panel("paralicheniformis", keep_legend = FALSE)
p_for_legend <- make_panel("paralicheniformis", keep_legend = TRUE)

legend <- get_legend(p_for_legend)

g_left <- ggplotGrob(p_left)
g_right <- ggplotGrob(p_right)

max_widths <- unit.pmax(g_left$widths, g_right$widths)
g_left$widths <- max_widths
g_right$widths <- max_widths

combined <- gtable(
  widths = unit.c(unit(1, "null"), unit(1, "null"), unit(2.3, "in")),
  heights = unit(1, "null")
)

combined <- gtable_add_grob(
  combined,
  g_left,
  t = 1, l = 1
)

combined <- gtable_add_grob(
  combined,
  g_right,
  t = 1, l = 2
)

combined <- gtable_add_grob(
  combined,
  legend,
  t = 1, l = 3
)

png(
  "figures/manuscript/Figure1_VolcanoPlots_Combined_publication.png",
  width = 16,
  height = 7,
  units = "in",
  res = 600
)
grid.newpage()
grid.draw(combined)
dev.off()

pdf(
  "figures/manuscript/Figure1_VolcanoPlots_Combined_publication.pdf",
  width = 16,
  height = 7
)
grid.newpage()
grid.draw(combined)
dev.off()

cat("Saved:\n")
cat("figures/manuscript/Figure1_VolcanoPlots_Combined_publication.png\n")
cat("figures/manuscript/Figure1_VolcanoPlots_Combined_publication.pdf\n")
