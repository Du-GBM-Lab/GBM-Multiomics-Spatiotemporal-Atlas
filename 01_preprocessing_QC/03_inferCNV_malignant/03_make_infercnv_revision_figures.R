suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(png)
  library(grid)
})

config <- list(
  infercnv_dir = normalizePath(file.path(".", "outputs"), winslash = "\\", mustWork = TRUE),
  score_csv = normalizePath(file.path("tables", "infercnv_cell_cnv_burden_and_calls.csv"), winslash = "\\", mustWork = TRUE),
  object_qs2 = normalizePath(file.path("outputs", "GBM.RNA.qc_doubletfinder.infercnv_calls.qs2"), winslash = "\\", mustWork = TRUE),
  gene_order_file = normalizePath(file.path("reference", "hg38_gencode_v27.txt"), winslash = "\\", mustWork = TRUE),
  out_dir = normalizePath(".", winslash = "\\", mustWork = TRUE),
  annotation_col = "anno_ident",
  selected_scatter_n = 6,
  max_cells_per_sample_heatmap = 500,
  min_malignant_cells_for_heatmap = 10,
  random_seed = 20260507
)

dir.create(file.path(config$out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

read_object <- function(path) {
  if (grepl("\\.qs2$", path, ignore.case = TRUE)) {
    return(qs2::qs_read(path))
  }
  readRDS(path)
}

find_infercnv_rds <- function(sample_dir) {
  candidates <- c(
    list.files(sample_dir, pattern = "_infercnv_obj\\.qs2$", full.names = TRUE),
    list.files(sample_dir, pattern = "_infercnv_obj\\.rds$", full.names = TRUE),
    list.files(sample_dir, pattern = "^run\\.final\\.infercnv_obj$", full.names = TRUE)
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) NA_character_ else candidates[1]
}

call_levels <- c(
  "non_malignant_reference",
  "non_malignant_like_CNV_low",
  "malignant_like_CNV_burden_only",
  "malignant_like_CNV_high_confidence",
  "undetermined_low_reference"
)

call_palette <- c(
  non_malignant_reference = "#587A8C",
  non_malignant_like_CNV_low = "#8E9A5B",
  malignant_like_CNV_burden_only = "#D8A03D",
  malignant_like_CNV_high_confidence = "#B24745",
  undetermined_low_reference = "#A7A7A7"
)

theme_cnv_pub <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold"),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.45, "cm")
    )
}

score_df <- read.csv(config$score_csv, stringsAsFactors = FALSE)
score_df$infercnv_call <- factor(score_df$infercnv_call, levels = call_levels)

## Main-ready two-axis scatter: choose representative samples by malignant fraction.
call_summary <- score_df |>
  count(sample, infercnv_call, name = "n_cells") |>
  group_by(sample) |>
  mutate(percent = n_cells / sum(n_cells) * 100) |>
  ungroup()

mal_pct <- call_summary |>
  filter(infercnv_call == "malignant_like_CNV_high_confidence") |>
  transmute(sample, malignant_percent = percent)

pick_closest <- function(target) {
  mal_pct$sample[which.min(abs(mal_pct$malignant_percent - target))]
}

selected_samples <- unique(c(
  mal_pct |> arrange(desc(malignant_percent)) |> slice_head(n = 2) |> pull(sample),
  mal_pct |> arrange(malignant_percent) |> slice_head(n = 2) |> pull(sample),
  pick_closest(stats::median(mal_pct$malignant_percent, na.rm = TRUE)),
  pick_closest(30)
))
selected_samples <- selected_samples[seq_len(min(length(selected_samples), config$selected_scatter_n))]

write.csv(
  mal_pct |> mutate(selected_for_main_scatter = sample %in% selected_samples),
  file.path(config$out_dir, "tables", "infercnv_main_scatter_selected_samples.csv"),
  row.names = FALSE
)

p_scatter_main <- score_df |>
  filter(sample %in% selected_samples, !is.na(cnv_burden_z), !is.na(cnv_correlation_ref_z)) |>
  ggplot(aes(x = cnv_burden_z, y = cnv_correlation_ref_z, color = infercnv_call)) +
  geom_point(size = 0.35, alpha = 0.55) +
  geom_vline(xintercept = 3, linetype = "dashed", linewidth = 0.25) +
  geom_hline(yintercept = 3, linetype = "dashed", linewidth = 0.25) +
  facet_wrap(~ sample, scales = "free", ncol = 3) +
  scale_color_manual(values = call_palette, drop = FALSE) +
  labs(x = "CNV burden z-score", y = "CNV correlation z-score", color = NULL) +
  theme_cnv_pub(base_size = 9) +
  theme(legend.position = "bottom")

ggsave(
  file.path(config$out_dir, "figures", "infercnv_two_axis_selected_samples_main.pdf"),
  p_scatter_main,
  width = 7.5,
  height = 5.2,
  units = "in",
  dpi = 300,
  device = grDevices::cairo_pdf
)

## DotPlot validation v2: grouped genes and publication palette.
msg("Redrawing marker validation DotPlot.")
obj <- qs2::qs_read(config$object_qs2)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}
if (!"data" %in% Layers(obj[["RNA"]])) {
  obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
}

gene_groups <- list(
  "Tumor lineage\n(non-specific)" = c("EGFR", "SOX2", "OLIG2", "GFAP", "PDGFRA", "CD44", "CHI3L1", "VIM"),
  "Proliferation" = c("MKI67"),
  "GBM CNV regions" = c("PTEN", "CDKN2A"),
  "Immune\n(reference)" = c("PTPRC", "CD3D", "NKG7", "CD68", "CX3CR1"),
  "Vascular\n(reference)" = c("PECAM1"),
  "Stress\n(sanity)" = c("HSPA1A", "HSPA1B", "JUN", "FOS")
)

gene_order <- unname(unlist(gene_groups))
gene_order <- intersect(gene_order, rownames(obj))
gene_group_df <- data.frame(
  features.plot = unname(unlist(gene_groups)),
  gene_group = rep(names(gene_groups), lengths(gene_groups)),
  stringsAsFactors = FALSE
) |>
  filter(features.plot %in% gene_order) |>
  mutate(
    features.plot = factor(features.plot, levels = gene_order),
    gene_group = factor(gene_group, levels = names(gene_groups))
  )

obj$infercnv_call <- factor(obj$infercnv_call, levels = call_levels)
dot_base <- DotPlot(
  obj,
  features = gene_order,
  group.by = "infercnv_call",
  assay = "RNA"
)

plot_data <- dot_base$data |>
  left_join(gene_group_df, by = "features.plot") |>
  mutate(
    features.plot = factor(features.plot, levels = gene_order),
    id = factor(id, levels = rev(call_levels))
  )

write.csv(plot_data, file.path(config$out_dir, "tables", "marker_validation_by_infercnv_call_dotplot_data.csv"), row.names = FALSE)

p_dot_v2 <- ggplot(plot_data, aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp)) +
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
  labs(x = NULL, y = NULL) +
  theme_cnv_pub(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(config$out_dir, "figures", "marker_validation_by_infercnv_call_v2.pdf"),
  p_dot_v2,
  width = 11,
  height = 3.4,
  units = "in",
  dpi = 300,
  device = grDevices::cairo_pdf
)

## Downgrade the too-dense split DotPlot to a supplementary table.
if ("infercnv_call_annotation" %in% colnames(obj@meta.data)) {
  obj$infercnv_call_annotation <- obj$infercnv_call_annotation
} else {
  obj$infercnv_call_annotation <- ifelse(
    is.na(obj$infercnv_call) | is.na(obj@meta.data[[config$annotation_col]]),
    NA_character_,
    paste(obj$infercnv_call, obj@meta.data[[config$annotation_col]], sep = " | ")
  )
}

group_sizes <- sort(table(obj$infercnv_call_annotation), decreasing = TRUE)
keep_groups <- names(group_sizes)[group_sizes >= 20]
if (length(keep_groups) > 0) {
  split_dot_data <- DotPlot(
    subset(obj, cells = rownames(obj@meta.data)[obj$infercnv_call_annotation %in% keep_groups]),
    features = gene_order,
    group.by = "infercnv_call_annotation",
    assay = "RNA"
  )$data
  write.csv(split_dot_data, file.path(config$out_dir, "tables", "marker_validation_by_infercnv_call_and_annotation_dotplot_data.csv"), row.names = FALSE)
}

split_pdf <- file.path(config$out_dir, "figures", "marker_validation_by_infercnv_call_and_annotation.pdf")
if (file.exists(split_pdf)) {
  unlink(split_pdf)
}

## Combined inferCNV chromosome heatmap for high-confidence malignant cells.
msg("Building combined inferCNV chromosome heatmap.")
set.seed(config$random_seed)
samples <- sort(unique(score_df$sample))
mal_expr_list <- list()
sample_cell_counts <- list()

for (sample_id in samples) {
  sample_dir <- file.path(config$infercnv_dir, sample_id)
  infercnv_rds <- find_infercnv_rds(sample_dir)
  if (is.na(infercnv_rds)) next

  mal_cells <- score_df$cell[
    score_df$sample == sample_id &
      score_df$infercnv_call == "malignant_like_CNV_high_confidence"
  ]
  if (length(mal_cells) < config$min_malignant_cells_for_heatmap) next

  infercnv_obj <- read_object(infercnv_rds)
  expr <- infercnv_obj@expr.data
  mal_cells <- intersect(mal_cells, colnames(expr))
  if (length(mal_cells) < config$min_malignant_cells_for_heatmap) next

  if (length(mal_cells) > config$max_cells_per_sample_heatmap) {
    mal_cells <- sample(mal_cells, config$max_cells_per_sample_heatmap)
  }

  expr_centered <- sweep(expr[, mal_cells, drop = FALSE], 1, rowMeans(expr), "-")
  colnames(expr_centered) <- paste(sample_id, colnames(expr_centered), sep = "__")
  mal_expr_list[[sample_id]] <- expr_centered
  sample_cell_counts[[sample_id]] <- data.frame(
    sample = sample_id,
    n_high_confidence_malignant = sum(score_df$sample == sample_id & score_df$infercnv_call == "malignant_like_CNV_high_confidence"),
    n_cells_used_for_heatmap = ncol(expr_centered),
    stringsAsFactors = FALSE
  )
}

if (length(mal_expr_list) < 2) {
  stop("Too few samples with high-confidence malignant cells for combined heatmap.", call. = FALSE)
}

common_genes <- Reduce(intersect, lapply(mal_expr_list, rownames))
gene_order_tbl <- read.table(
  config$gene_order_file,
  sep = "\t",
  header = FALSE,
  col.names = c("gene", "chr", "start", "end"),
  stringsAsFactors = FALSE
) |>
  filter(gene %in% common_genes, chr %in% paste0("chr", c(1:22, "X", "Y"))) |>
  distinct(gene, .keep_all = TRUE) |>
  mutate(chr = factor(chr, levels = paste0("chr", c(1:22, "X", "Y")))) |>
  arrange(chr, start)

common_genes <- gene_order_tbl$gene
mal_expr_combined <- do.call(cbind, lapply(mal_expr_list, function(x) x[common_genes, , drop = FALSE]))
mal_expr_combined <- mal_expr_combined[gene_order_tbl$gene, , drop = FALSE]
mal_expr_combined <- pmax(pmin(mal_expr_combined, 0.35), -0.35)

sample_anno <- sub("__.*$", "", colnames(mal_expr_combined))
sample_levels <- unique(sample_anno)
sample_cols <- setNames(colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(sample_levels)), sample_levels)

write.csv(
  bind_rows(sample_cell_counts),
  file.path(config$out_dir, "tables", "infercnv_combined_heatmap_cells_used_by_sample.csv"),
  row.names = FALSE
)

pdf(
  file.path(config$out_dir, "figures", "infercnv_combined_chromosome_heatmap_high_confidence_malignant.pdf"),
  width = 12,
  height = 7.5,
  useDingbats = FALSE
)
draw(
  Heatmap(
    mal_expr_combined,
    name = "CNV\nsignal",
    col = circlize::colorRamp2(c(-0.3, 0, 0.3), c("#2C4C9A", "#F7F7F7", "#B2182B")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_split = gene_order_tbl$chr,
    column_split = factor(sample_anno, levels = sample_levels),
    top_annotation = HeatmapAnnotation(
      sample = factor(sample_anno, levels = sample_levels),
      col = list(sample = sample_cols),
      show_annotation_name = FALSE
    ),
    use_raster = TRUE,
    raster_quality = 2,
    row_title_gp = grid::gpar(fontsize = 6),
    column_title_gp = grid::gpar(fontsize = 7),
    heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 8), labels_gp = grid::gpar(fontsize = 7))
  ),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

## Supplementary per-sample native inferCNV heatmap manifest.
native_png_manifest <- lapply(samples, function(sample_id) {
  sample_dir <- file.path(config$infercnv_dir, sample_id)
  data.frame(
    sample = sample_id,
    infercnv_png = file.path(sample_dir, "infercnv.png"),
    infercnv_subclusters_png = file.path(sample_dir, "infercnv_subclusters.png"),
    infercnv_png_exists = file.exists(file.path(sample_dir, "infercnv.png")),
    infercnv_subclusters_png_exists = file.exists(file.path(sample_dir, "infercnv_subclusters.png")),
    stringsAsFactors = FALSE
  )
}) |> bind_rows()

write.csv(native_png_manifest, file.path(config$out_dir, "tables", "infercnv_native_heatmap_png_manifest.csv"), row.names = FALSE)

native_png_manifest <- native_png_manifest |> filter(infercnv_png_exists)
if (nrow(native_png_manifest) > 0) {
  pdf(
    file.path(config$out_dir, "figures", "infercnv_native_chromosome_heatmaps_per_sample_supplementary.pdf"),
    width = 11,
    height = 8.5,
    useDingbats = FALSE
  )
  for (i in seq_len(nrow(native_png_manifest))) {
    img <- png::readPNG(native_png_manifest$infercnv_png[i])
    grid::grid.newpage()
    grid::grid.text(
      native_png_manifest$sample[i],
      x = 0.02,
      y = 0.98,
      just = c("left", "top"),
      gp = grid::gpar(fontsize = 12, fontface = "bold")
    )
    grid::grid.raster(img, x = 0.5, y = 0.48, width = 0.98, height = 0.9, interpolate = TRUE)
  }
  dev.off()
}

msg("Done.")
