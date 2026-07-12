suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
})

config <- list(
  input_qs2 = normalizePath(file.path("outputs", "GBM.RNA.qc_doubletfinder.infercnv_immune_reference_calls.qs2"), winslash = "\\", mustWork = FALSE),
  out_dir = normalizePath(".", winslash = "\\", mustWork = TRUE),
  annotation_col = "anno_ident",
  min_cells_per_group = 20,
  copy_to_figure_table = TRUE,
  figure_table_dir = normalizePath(file.path("..", "..", "图片表格"), winslash = "\\", mustWork = FALSE)
)

dir.create(file.path(config$out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

theme_marker_pub <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = base_size - 1),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.key.width = grid::unit(0.45, "cm")
    )
}

make_dotplot_data <- function(obj, features, group_by) {
  features <- intersect(features, rownames(obj))
  if (length(features) == 0) {
    stop("None of the requested marker genes are present in the object.", call. = FALSE)
  }
  DotPlot(obj, features = features, group.by = group_by, assay = "RNA")$data
}

make_gene_group_df <- function(gene_groups, gene_order) {
  bind_rows(lapply(names(gene_groups), function(group_name) {
    data.frame(
      features.plot = gene_groups[[group_name]],
      gene_group = group_name,
      stringsAsFactors = FALSE
    )
  })) |>
    filter(features.plot %in% gene_order) |>
    distinct(features.plot, .keep_all = TRUE) |>
    mutate(
      features.plot = factor(features.plot, levels = gene_order),
      gene_group = factor(gene_group, levels = names(gene_groups))
    )
}

plot_grouped_bubble <- function(plot_data, gene_group_df, group_levels, file, width, height) {
  plot_data <- plot_data |>
    left_join(gene_group_df, by = "features.plot") |>
    mutate(
      features.plot = factor(features.plot, levels = levels(gene_group_df$features.plot)),
      gene_group = factor(gene_group, levels = levels(gene_group_df$gene_group)),
      id = factor(id, levels = rev(group_levels))
    )

  p <- ggplot(plot_data, aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp)) +
    geom_point(alpha = 0.95) +
    facet_grid(. ~ gene_group, scales = "free_x", space = "free_x") +
    scale_color_gradient2(
      low = "#3B4CC0",
      mid = "#F5F5F5",
      high = "#B40426",
      midpoint = 0,
      name = "Avg expr\n(z-scored)"
    ) +
    scale_size(range = c(0, 5.8), limits = c(0, 100), name = "% expressed") +
    labs(
      x = NULL,
      y = NULL,
      caption = "Color values are z-scored within this panel's grouping set."
    ) +
    theme_marker_pub(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8.5),
      plot.caption = element_text(size = 7, color = "grey35", hjust = 0),
      legend.position = "right"
    )

  ggsave(file, p, width = width, height = height, units = "in", dpi = 300, device = grDevices::cairo_pdf)
  plot_data
}

if (!file.exists(config$input_qs2)) {
  stop("Immune-reference inferCNV object not found. Run 02_summarize_infercnv_immune_reference_calls.R first: ", config$input_qs2, call. = FALSE)
}

msg("Loading immune-reference inferCNV object:", config$input_qs2)
obj <- qs2::qs_read(config$input_qs2)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}
if (!"data" %in% Layers(obj[["RNA"]])) {
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
}

stopifnot(config$annotation_col %in% colnames(obj@meta.data))
stopifnot("infercnv_immune_call" %in% colnames(obj@meta.data))
stopifnot("is_malignant_for_downstream_immune" %in% colnames(obj@meta.data))

## Plot 1: original annotation marker bubble plot.
annotation_gene_groups <- list(
  "Malignant/CNV\ncandidate" = c("EGFR", "PDGFRA", "SOX2", "OLIG2", "GFAP", "NES"),
  "GBM state /\ninvasion" = c("CD44", "CHI3L1", "VIM", "COL1A1", "MMP9", "PLAU", "PLAUR"),
  "Proliferation" = c("MKI67", "TOP2A"),
  "Astro / OPC /\nOligodendrocyte" = c("AQP4", "ALDH1L1", "CSPG4", "MBP", "MAG", "MOBP"),
  "Neuronal /\nradial glial" = c("RBFOX3", "SNAP25", "DCX", "HOPX", "VIM"),
  "Immune /\nmyeloid" = c("PTPRC", "CD3D", "NKG7", "CD68", "CX3CR1", "AIF1", "LST1"),
  "Vascular /\nmural" = c("PECAM1", "VWF", "PDGFRB", "RGS5")
)

annotation_gene_order <- unique(unname(unlist(annotation_gene_groups)))
annotation_gene_order <- intersect(annotation_gene_order, rownames(obj))
annotation_gene_group_df <- make_gene_group_df(annotation_gene_groups, annotation_gene_order)

annotation_counts <- table(obj@meta.data[[config$annotation_col]])
annotation_levels <- names(sort(annotation_counts, decreasing = TRUE))
obj@meta.data[[config$annotation_col]] <- factor(obj@meta.data[[config$annotation_col]], levels = annotation_levels)

annotation_plot_data <- make_dotplot_data(obj, annotation_gene_order, config$annotation_col)
annotation_plot_data <- plot_grouped_bubble(
  annotation_plot_data,
  annotation_gene_group_df,
  annotation_levels,
  file.path(config$out_dir, "figures", "marker_bubble_original_annotation_lineage_malignancy_reference.pdf"),
  width = 14,
  height = 6.8
)
write.csv(annotation_plot_data, file.path(config$out_dir, "tables", "marker_bubble_original_annotation_lineage_malignancy_reference_data.csv"), row.names = FALSE)

## Plot 2: malignant high-confidence vs non-malignant/CNV-low validation.
obj$malignant_validation_group <- dplyr::case_when(
  obj$is_malignant_for_downstream_immune %in% TRUE ~ "Malignant high-confidence",
  obj$infercnv_immune_call %in% c("non_malignant_reference", "non_malignant_like_CNV_low") ~ "Non-malignant / CNV-low",
  obj$infercnv_immune_call == "malignant_like_CNV_burden_only" ~ "Ambiguous burden-only",
  obj$infercnv_immune_call == "undetermined_low_reference" ~ "Undetermined",
  TRUE ~ NA_character_
)

main_cells <- colnames(obj)[obj$malignant_validation_group %in% c("Malignant high-confidence", "Non-malignant / CNV-low")]
obj_main <- subset(obj, cells = main_cells)
obj_main$malignant_validation_group <- factor(
  obj_main$malignant_validation_group,
  levels = c("Non-malignant / CNV-low", "Malignant high-confidence")
)

validation_gene_groups <- list(
  "Tumor lineage\n(non-specific)" = c("EGFR", "PDGFRA", "SOX2", "OLIG2", "GFAP", "NES"),
  "GBM state /\ninvasion" = c("CD44", "CHI3L1", "VIM", "COL1A1", "MMP9", "PLAU", "PLAUR"),
  "Proliferation" = c("MKI67", "TOP2A"),
  "GBM CNV\nregions" = c("PTEN", "CDKN2A"),
  "Immune /\nvascular ref" = c("PTPRC", "CD3D", "NKG7", "CD68", "CX3CR1", "PECAM1", "VWF"),
  "Stress\nsanity" = c("HSPA1A", "HSPA1B", "JUN", "FOS")
)

validation_gene_order <- unique(unname(unlist(validation_gene_groups)))
validation_gene_order <- intersect(validation_gene_order, rownames(obj_main))
validation_gene_group_df <- make_gene_group_df(validation_gene_groups, validation_gene_order)

validation_plot_data <- make_dotplot_data(obj_main, validation_gene_order, "malignant_validation_group")
validation_plot_data <- plot_grouped_bubble(
  validation_plot_data,
  validation_gene_group_df,
  levels(obj_main$malignant_validation_group),
  file.path(config$out_dir, "figures", "marker_bubble_malignant_high_confidence_vs_non_malignant.pdf"),
  width = 12,
  height = 3.1
)
write.csv(validation_plot_data, file.path(config$out_dir, "tables", "marker_bubble_malignant_high_confidence_vs_non_malignant_data.csv"), row.names = FALSE)

heat_data <- validation_plot_data |>
  select(features.plot, id, avg.exp.scaled, pct.exp, gene_group) |>
  mutate(
    features.plot = factor(features.plot, levels = validation_gene_order),
    id = factor(id, levels = levels(obj_main$malignant_validation_group))
  )

p_heat <- ggplot(heat_data, aes(x = features.plot, y = id, fill = avg.exp.scaled)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.0f%%", pct.exp)), size = 2.4, color = "black") +
  facet_grid(. ~ gene_group, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "#F5F5F5",
    high = "#B40426",
    midpoint = 0,
    name = "Avg expr\n(z-scored)"
  ) +
  labs(
    x = NULL,
    y = NULL,
    caption = "Color values are z-scored within this panel's grouping set; tile text shows % expressed."
  ) +
  theme_marker_pub(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    plot.caption = element_text(size = 7, color = "grey35", hjust = 0),
    legend.position = "right"
  )

ggsave(
  file.path(config$out_dir, "figures", "marker_heatmap_malignant_high_confidence_vs_non_malignant.pdf"),
  p_heat,
  width = 12,
  height = 2.7,
  units = "in",
  dpi = 300,
  device = grDevices::cairo_pdf
)

group_summary <- obj@meta.data |>
  count(malignant_validation_group, name = "n_cells") |>
  mutate(percent = n_cells / sum(n_cells) * 100)
write.csv(group_summary, file.path(config$out_dir, "tables", "malignant_validation_group_summary.csv"), row.names = FALSE)

if (config$copy_to_figure_table) {
  msg("Figure/table mirror dir:", config$figure_table_dir, "exists:", dir.exists(config$figure_table_dir))
}

if (config$copy_to_figure_table && dir.exists(config$figure_table_dir)) {
  fig_main <- file.path(config$figure_table_dir, "Figures")
  fig_sup <- file.path(config$figure_table_dir, "Supplementary_Figures")
  tab_sup <- file.path(config$figure_table_dir, "Supplementary_Tables")
  dir.create(fig_main, showWarnings = FALSE, recursive = TRUE)
  dir.create(fig_sup, showWarnings = FALSE, recursive = TRUE)
  dir.create(tab_sup, showWarnings = FALSE, recursive = TRUE)

  copy_map <- c(
    marker_bubble_original_annotation_lineage_malignancy_reference = "Supp_marker_bubble_original_annotation_lineage_malignancy_reference.pdf",
    marker_bubble_malignant_high_confidence_vs_non_malignant = "Fig_marker_bubble_malignant_high_confidence_vs_non_malignant.pdf",
    marker_heatmap_malignant_high_confidence_vs_non_malignant = "Fig_marker_heatmap_malignant_high_confidence_vs_non_malignant.pdf"
  )
  for (nm in names(copy_map)) {
    src <- file.path(config$out_dir, "figures", paste0(nm, ".pdf"))
    dst_dir <- if (startsWith(copy_map[[nm]], "Supp_")) fig_sup else fig_main
    if (file.exists(src)) {
      file.copy(src, file.path(dst_dir, copy_map[[nm]]), overwrite = TRUE)
    }
  }

  table_files <- c(
    "marker_bubble_original_annotation_lineage_malignancy_reference_data.csv",
    "marker_bubble_malignant_high_confidence_vs_non_malignant_data.csv",
    "malignant_validation_group_summary.csv"
  )
  for (table_file in table_files) {
    src <- file.path(config$out_dir, "tables", table_file)
    if (file.exists(src)) {
      file.copy(src, file.path(tab_sup, table_file), overwrite = TRUE)
    }
  }
}

msg("Done.")
