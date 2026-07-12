suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

config <- list(
  table_dir = file.path(".", "tables"),
  figure_dir = file.path(".", "figures")
)

theme_qc_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size * 0.9),
      legend.key.size = grid::unit(0.38, "cm"),
      legend.position = "right",
      plot.margin = margin(6, 8, 6, 6)
    )
}

save_plot_pair <- function(plot, filename, width, height) {
  pdf_path <- file.path(config$figure_dir, paste0(filename, ".pdf"))
  png_path <- file.path(config$figure_dir, paste0(filename, ".png"))
  if (file.exists(pdf_path)) file.remove(pdf_path)
  if (file.exists(png_path)) file.remove(png_path)
  ggsave(pdf_path, plot, width = width, height = height, units = "in", device = cairo_pdf)
  ggsave(png_path, plot, width = width, height = height, units = "in", dpi = 300)
}

plot_stacked_cells <- function(df, x, y, fill, palette, y_label = "Cells") {
  df[[fill]] <- factor(df[[fill]], levels = names(palette))
  ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    geom_col(width = 0.76, color = "white", linewidth = 0.16) +
    scale_fill_manual(values = palette, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(x = NULL, y = y_label) +
    theme_qc_pub()
}

plot_stacked_percent <- function(df, x, y, fill, palette, y_label = "Cell proportion after QC (%)") {
  df[[fill]] <- factor(df[[fill]], levels = names(palette))
  ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    geom_col(width = 0.76, color = "white", linewidth = 0.12) +
    scale_fill_manual(values = palette, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = y_label) +
    theme_qc_pub(base_size = 9) +
    theme(legend.text = element_text(size = 8.2))
}

qc_status_palette <- c(
  Pass_pre_doublet = "#2F6F73",
  Fail_MAD = "#C69214",
  Fail_MT20 = "#B24745",
  Fail_other = "#9A9A9A"
)

doublet_palette <- c(
  Singlet = "#2F6F73",
  Doublet = "#B24745",
  Skipped_low_cells = "#8E8E8E"
)

final_status_palette <- c(
  Final_pass = "#2F6F73",
  Filtered = "#B8A27A"
)

celltype_palette <- c(
  "T cells" = "#3E6C8E",
  "Macrophages" = "#B45D4C",
  "Radial glial" = "#6D8F4E",
  "Microglial" = "#8C6BB1",
  "Mural cells" = "#C58A2F",
  "OPCs" = "#4C9A9A",
  "Oligodendrocytes" = "#7A6A55",
  "Ambiguous" = "#9A9A9A",
  "Astrocytes" = "#5A7FBB",
  "NK cells" = "#D08A7A",
  "Monocytes" = "#7C9D68",
  "Endothelial" = "#A66D9B",
  "Neurons" = "#C9A227",
  "cDCs" = "#6EA6B8",
  "B cells" = "#8A7AA8",
  "Ependymal cells" = "#B18462",
  "pDCs" = "#5F6F80"
)

qc_summary <- read.csv(file.path(config$table_dir, "QC_filtering_summary_by_sample.csv"), check.names = FALSE)
doublet_summary <- read.csv(file.path(config$table_dir, "DoubletFinder_summary_by_sample.csv"), check.names = FALSE)
final_summary <- read.csv(file.path(config$table_dir, "final_filtering_summary_by_sample.csv"), check.names = FALSE)
composition <- read.csv(file.path(config$table_dir, "cell_type_composition_after_QC.csv"), check.names = FALSE)

sample_levels <- mixedsort_unique <- function(x) {
  x[order(as.integer(sub("^Pt", "", x)))]
}
sample_levels <- mixedsort_unique(unique(qc_summary$sample))

qc_summary$sample <- factor(qc_summary$sample, levels = sample_levels)
doublet_summary$sample <- factor(doublet_summary$sample, levels = sample_levels)
final_summary$sample <- factor(final_summary$sample, levels = sample_levels)
composition$sample <- factor(composition$sample, levels = sample_levels)

save_plot_pair(
  plot_stacked_cells(qc_summary, "sample", "n_cells", "qc_status", qc_status_palette),
  "QC_filtering_summary_by_sample",
  width = 8.4,
  height = 4.6
)

save_plot_pair(
  plot_stacked_cells(doublet_summary, "sample", "n_cells", "doubletfinder_class", doublet_palette),
  "DoubletFinder_summary_by_sample",
  width = 8.4,
  height = 4.6
)

save_plot_pair(
  plot_stacked_cells(final_summary, "sample", "n_cells", "final_status", final_status_palette),
  "final_filtering_summary_by_sample",
  width = 8.4,
  height = 4.6
)

save_plot_pair(
  plot_stacked_percent(composition, "sample", "percent", "annotation", celltype_palette),
  "cell_type_composition_after_QC_percent",
  width = 10.2,
  height = 5.2
)

message("Publication-style QC barplots redrawn.")
