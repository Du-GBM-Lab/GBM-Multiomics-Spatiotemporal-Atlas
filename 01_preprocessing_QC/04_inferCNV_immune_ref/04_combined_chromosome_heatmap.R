## Build a publication-quality combined chromosome heatmap from sample-wise
## immune-reference inferCNV outputs.

suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(RColorBrewer)
})

config <- list(
  input_qs2 = normalizePath(file.path("outputs", "GBM.RNA.qc_doubletfinder.infercnv_immune_reference_calls.qs2"), winslash = "\\", mustWork = FALSE),
  infercnv_dir = normalizePath(file.path("outputs"), winslash = "\\", mustWork = TRUE),
  out_dir = normalizePath(".", winslash = "\\", mustWork = TRUE),
  sample_col = "Pt_number",
  call_col = "infercnv_immune_call",
  reference_call = "non_malignant_reference",
  malignant_call = "malignant_like_CNV_high_confidence",
  max_cells_per_sample_malignant = 120,
  max_cells_per_sample_reference = 60,
  gene_min_sample_fraction = 0.8,
  heatmap_color_breaks = c(-0.28, -0.14, 0, 0.14, 0.28),
  heatmap_palette = c("#174A8B", "#78ADD2", "white", "#E89572", "#992C24"),
  raster_device = if (requireNamespace("Cairo", quietly = TRUE)) "CairoPNG" else "png",
  raster_quality = 3,
  fig_width_in = 11.2,
  fig_height_in = 5.6,
  random_seed = 20260507
)

dir.create(file.path(config$out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

set.seed(config$random_seed)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

read_object <- function(path) {
  if (grepl("\\.qs2$", path, ignore.case = TRUE)) qs2::qs_read(path) else readRDS(path)
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

standardize_gene_order <- function(gene_order) {
  gene_order <- as.data.frame(gene_order)
  if (!"gene" %in% colnames(gene_order)) {
    gene_order$gene <- rownames(gene_order)
  }
  if (!all(c("chr", "start") %in% colnames(gene_order))) {
    if (ncol(gene_order) >= 3) {
      colnames(gene_order)[seq_len(min(3, ncol(gene_order)))] <- c("chr", "start", "stop")[seq_len(min(3, ncol(gene_order)))]
    }
  }
  if (!all(c("gene", "chr", "start") %in% colnames(gene_order))) {
    stop("inferCNV gene_order must contain gene/chr/start information.", call. = FALSE)
  }
  gene_order
}

align_block <- function(mat, genes) {
  missing <- setdiff(genes, rownames(mat))
  if (length(missing) > 0) {
    pad <- matrix(0, nrow = length(missing), ncol = ncol(mat), dimnames = list(missing, colnames(mat)))
    mat <- rbind(mat, pad)
  }
  mat[genes, , drop = FALSE]
}

if (!file.exists(config$input_qs2)) {
  stop("Immune-reference call object not found: ", config$input_qs2, call. = FALSE)
}

msg("Loading immune-reference call object:", config$input_qs2)
obj <- read_object(config$input_qs2)
md <- obj@meta.data
stopifnot(config$sample_col %in% colnames(md))
stopifnot(config$call_col %in% colnames(md))

md$cell <- rownames(md)
md$sample <- as.character(md[[config$sample_col]])
md$call <- as.character(md[[config$call_col]])
samples <- sort(unique(md$sample))

msg("Collecting and centering per-sample inferCNV matrices.")
expr_blocks_malignant <- list()
expr_blocks_reference <- list()
sample_summary <- list()
gene_order_master <- NULL

for (sample_id in samples) {
  sample_dir <- file.path(config$infercnv_dir, sample_id)
  infercnv_path <- find_infercnv_rds(sample_dir)
  if (is.na(infercnv_path)) {
    msg("[SKIP] No inferCNV object for", sample_id)
    next
  }

  infercnv_obj <- read_object(infercnv_path)
  expr <- infercnv_obj@expr.data
  if (!is.matrix(expr)) expr <- as.matrix(expr)

  if (is.null(gene_order_master) && !is.null(infercnv_obj@gene_order) && nrow(infercnv_obj@gene_order) > 0) {
    gene_order_master <- standardize_gene_order(infercnv_obj@gene_order)
  }

  cells_in_expr <- colnames(expr)
  sample_md <- md[md$sample == sample_id & md$cell %in% cells_in_expr, , drop = FALSE]
  ref_cells <- sample_md$cell[sample_md$call == config$reference_call]
  mal_cells <- sample_md$cell[sample_md$call == config$malignant_call]

  if (length(ref_cells) < 10) {
    msg("[SKIP]", sample_id, "has too few reference cells (", length(ref_cells), ")")
    next
  }

  ref_mean <- rowMeans(expr[, ref_cells, drop = FALSE])
  expr_centered <- sweep(expr, 1, ref_mean, "-")

  if (length(ref_cells) > config$max_cells_per_sample_reference) {
    ref_cells <- sample(ref_cells, config$max_cells_per_sample_reference)
  }
  if (length(mal_cells) > config$max_cells_per_sample_malignant) {
    mal_cells <- sample(mal_cells, config$max_cells_per_sample_malignant)
  }

  if (length(ref_cells) > 0) {
    ref_block <- expr_centered[, ref_cells, drop = FALSE]
    colnames(ref_block) <- paste(sample_id, colnames(ref_block), sep = "__")
    expr_blocks_reference[[sample_id]] <- ref_block
  }
  if (length(mal_cells) > 0) {
    mal_block <- expr_centered[, mal_cells, drop = FALSE]
    colnames(mal_block) <- paste(sample_id, colnames(mal_block), sep = "__")
    expr_blocks_malignant[[sample_id]] <- mal_block
  }

  sample_summary[[sample_id]] <- data.frame(
    sample = sample_id,
    n_genes_in_expr = nrow(expr),
    n_reference_used = length(ref_cells),
    n_malignant_used = length(mal_cells),
    stringsAsFactors = FALSE
  )
}

if (length(expr_blocks_malignant) == 0) {
  stop("No malignant cells available for the combined heatmap. Check inferCNV outputs.", call. = FALSE)
}
if (is.null(gene_order_master)) {
  stop("Could not retrieve @gene_order from any inferCNV object.", call. = FALSE)
}

sample_summary_df <- bind_rows(sample_summary)
write.csv(sample_summary_df, file.path(config$out_dir, "tables", "combined_chromosome_heatmap_sample_summary.csv"), row.names = FALSE)

all_blocks <- c(expr_blocks_reference, expr_blocks_malignant)
gene_lists <- lapply(all_blocks, rownames)
gene_universe <- unique(unlist(gene_lists))
gene_coverage <- vapply(
  gene_universe,
  function(g) mean(vapply(gene_lists, function(gl) g %in% gl, logical(1))),
  numeric(1)
)
common_genes <- names(gene_coverage)[gene_coverage >= config$gene_min_sample_fraction]
msg(sprintf(
  "Gene coverage filter: %d / %d genes pass (>= %.0f%% of blocks).",
  length(common_genes), length(gene_universe), 100 * config$gene_min_sample_fraction
))

gene_order_use <- gene_order_master |>
  filter(gene %in% common_genes, chr %in% paste0("chr", c(1:22, "X", "Y"))) |>
  distinct(gene, .keep_all = TRUE) |>
  mutate(chr = factor(chr, levels = paste0("chr", c(1:22, "X", "Y")))) |>
  arrange(chr, start)

common_genes <- gene_order_use$gene
msg(sprintf("Final gene set: %d genes across %d chromosomes.", length(common_genes), length(unique(gene_order_use$chr))))

ref_mat <- if (length(expr_blocks_reference) > 0) {
  do.call(cbind, lapply(expr_blocks_reference, align_block, genes = common_genes))
} else {
  matrix(numeric(0), nrow = length(common_genes), ncol = 0, dimnames = list(common_genes, NULL))
}
mal_mat <- do.call(cbind, lapply(expr_blocks_malignant, align_block, genes = common_genes))

ref_sample_anno <- rep(names(expr_blocks_reference), vapply(expr_blocks_reference, ncol, integer(1)))
mal_sample_anno <- rep(names(expr_blocks_malignant), vapply(expr_blocks_malignant, ncol, integer(1)))

expr_combined <- cbind(ref_mat, mal_mat)
cell_sample <- c(ref_sample_anno, mal_sample_anno)
cell_type <- c(rep("Reference", ncol(ref_mat)), rep("Malignant", ncol(mal_mat)))
expr_combined <- pmax(pmin(expr_combined, max(config$heatmap_color_breaks)), min(config$heatmap_color_breaks))

stopifnot(ncol(expr_combined) == length(cell_sample))
stopifnot(nrow(expr_combined) == length(common_genes))

msg(sprintf(
  "Combined matrix: %d genes x %d cells (%d ref, %d malignant).",
  nrow(expr_combined), ncol(expr_combined), ncol(ref_mat), ncol(mal_mat)
))

chr_factor <- factor(gene_order_use$chr, levels = paste0("chr", c(1:22, "X", "Y")))
sample_levels <- samples
sample_palette <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(length(sample_levels)), sample_levels)
type_palette <- c(Reference = "#4F6F7F", Malignant = "#B24745")

row_anno <- rowAnnotation(
  Type = factor(cell_type, levels = c("Reference", "Malignant")),
  col = list(Type = type_palette),
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 8),
  simple_anno_size = unit(3, "mm"),
  show_legend = c(Type = TRUE)
)

col_fun <- colorRamp2(config$heatmap_color_breaks, config$heatmap_palette)
mat_for_plot <- t(expr_combined)
row_split_factor <- factor(cell_type, levels = c("Reference", "Malignant"))

ht <- Heatmap(
  mat_for_plot,
  name = "Modified\nexpression",
  col = col_fun,
  cluster_columns = FALSE,
  column_split = chr_factor,
  column_title_gp = gpar(fontsize = 5.4, fontface = "bold"),
  column_title_rot = 0,
  column_gap = unit(0.28, "mm"),
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  row_split = row_split_factor,
  row_title_rot = 90,
  row_title_gp = gpar(fontsize = 8.5, fontface = "bold"),
  row_gap = unit(1.8, "mm"),
  show_row_names = FALSE,
  show_row_dend = FALSE,
  left_annotation = row_anno,
  use_raster = TRUE,
  raster_quality = config$raster_quality,
  raster_device = config$raster_device,
  raster_resize_mat = TRUE,
  border = FALSE,
  heatmap_legend_param = list(
    title = "Relative\nexpression",
    direction = "horizontal",
    legend_width = unit(3.6, "cm"),
    at = c(-0.28, -0.14, 0, 0.14, 0.28),
    labels = c("-0.28", "-0.14", "0", "0.14", "0.28"),
    title_gp = gpar(fontsize = 8),
    labels_gp = gpar(fontsize = 7)
  )
)

out_pdf <- file.path(config$out_dir, "figures", "Fig_combined_chromosome_heatmap_24samples.pdf")
out_png <- file.path(config$out_dir, "figures", "Fig_combined_chromosome_heatmap_24samples.png")

cairo_pdf(out_pdf, width = config$fig_width_in, height = config$fig_height_in)
draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)
dev.off()

if (FALSE && requireNamespace("pdftools", quietly = TRUE)) {
  pdftools::pdf_convert(out_pdf, format = "png", dpi = 300, filenames = out_png)
} else {
  msg("pdftools is not installed; skipping PNG conversion and keeping PDF only.")
}

cell_table <- data.frame(
  cell = colnames(expr_combined),
  sample = cell_sample,
  cell_type = cell_type,
  stringsAsFactors = FALSE
)
write.csv(cell_table, file.path(config$out_dir, "tables", "combined_chromosome_heatmap_cell_table.csv"), row.names = FALSE)
write.csv(gene_order_use, file.path(config$out_dir, "tables", "combined_chromosome_heatmap_gene_order.csv"), row.names = FALSE)

msg("Wrote:", out_pdf)
if (file.exists(out_png)) msg("Wrote:", out_png)
msg("Done.")
