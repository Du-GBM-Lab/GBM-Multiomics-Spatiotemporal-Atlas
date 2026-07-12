options(stringsAsFactors = FALSE)
set.seed(42)

if (basename(getwd()) == "06_恶性细胞拟时序") {
  project_root <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(file.path(project_root, "06_恶性细胞拟时序"))

dir.create("figures/final_panels/supplement", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/final_panels/source_data/supplement", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(qs2)
  library(Seurat)
  library(Matrix)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

log_file <- "logs/22_make_S7A_three_lineage_driver_heatmap.log"
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(...))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

lineage_labels <- c(
  "lineage1" = "L1: NPC-P to MES-I",
  "lineage2" = "L2: NPC-P to MES-V",
  "lineage3" = "L3: NPC-P to OPC-M"
)
lineage_colors <- c(
  "L1: NPC-P to MES-I" = "#BC3C29",
  "L2: NPC-P to MES-V" = "#20854E",
  "L3: NPC-P to OPC-M" = "#E18727"
)
obj_path <- normalizePath(file.path("..", "05_恶性细胞分亚群与Neftel对照", "outputs", "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"), winslash = "/", mustWork = TRUE)
pseudotime_path <- file.path("tables", "slingshot_pseudotime_per_cell.csv")

assoc <- read_csv("tables/tradeseq_association_per_lineage.csv", show_col_types = FALSE) |>
  filter(p_adj_BH < 0.05, !is.na(pattern), pattern != "flat_or_subtle") |>
  mutate(
    lineage_label = recode(lineage, !!!lineage_labels)
  )

binned_raw <- read_csv("figures/source_data/08_binned_expression_for_patterns.csv", show_col_types = FALSE) |>
  inner_join(assoc, by = c("gene", "lineage"))

candidate_metrics <- binned_raw |>
  group_by(gene, lineage, lineage_label) |>
  summarise(
    early_mean = mean(expr_mean[bin <= 4], na.rm = TRUE),
    late_mean = mean(expr_mean[bin >= 17], na.rm = TRUE),
    mid_mean = mean(expr_mean[bin >= 9 & bin <= 12], na.rm = TRUE),
    min_expr = min(expr_mean, na.rm = TRUE),
    max_expr = max(expr_mean, na.rm = TRUE),
    dynamic_range = max_expr - min_expr,
    early_late_delta = late_mean - early_mean,
    abs_delta = abs(early_late_delta),
    mean_expr = mean(expr_mean, na.rm = TRUE),
    waldStat = first(waldStat),
    p_adj_BH = first(p_adj_BH),
    pattern = first(pattern),
    .groups = "drop"
  ) |>
  mutate(
    timing = if_else(early_late_delta < 0, "early_high", "late_high"),
    endpoint_balance = abs_delta / pmax(dynamic_range, 1e-6),
    visual_score = dynamic_range * 0.62 + abs_delta * 0.20 + log1p(waldStat) * 0.13 + endpoint_balance * 0.05,
    lineage_short = recode(lineage, lineage1 = "L1", lineage2 = "L2", lineage3 = "L3")
  ) |>
  arrange(lineage_short, desc(visual_score))

write_csv(
  candidate_metrics,
  "figures/final_panels/source_data/supplement/S7A_driver_gene_visual_candidates.csv"
)

candidate_pool <- candidate_metrics |>
  mutate(
    row_id = paste(lineage, gene, sep = "|"),
    display_gene = gene
  )

log_msg("Loading object for high-resolution expression binning:", obj_path)
obj <- qs2::qs_read(obj_path)
pst <- readr::read_csv(pseudotime_path, show_col_types = FALSE)
cells <- intersect(colnames(obj), pst$cellID)
pst <- pst |>
  filter(.data$cellID %in% cells) |>
  arrange(match(.data$cellID, cells))
stopifnot(identical(pst$cellID, cells))

n_real_bins <- 120
probe_genes <- unique(c(candidate_pool$gene, "PLAUR"))
probe_genes <- intersect(probe_genes, rownames(obj))
expr_mat <- SeuratObject::LayerData(obj[["RNA"]], layer = "data")[probe_genes, cells, drop = FALSE]

make_lineage_bins <- function(lin, n_bins = n_real_bins) {
  pst_col <- paste0(lin, "_pst")
  weight_col <- paste0(lin, "_weight")
  pst |>
    filter(.data$lineage_assignment == lin, .data[[weight_col]] > 0.5, is.finite(.data[[pst_col]])) |>
    transmute(
      cellID = .data$cellID,
      lineage = lin,
      pseudotime = .data[[pst_col]],
      bin = dplyr::ntile(.data[[pst_col]], n_bins)
    )
}

cell_bins <- lapply(paste0("lineage", 1:3), make_lineage_bins) |>
  bind_rows()

make_gene_bins <- function(lin) {
  cb <- cell_bins |> filter(.data$lineage == lin)
  split_cells <- split(cb$cellID, cb$bin)
  bind_rows(lapply(names(split_cells), function(bn) {
    cs <- split_cells[[bn]]
    tibble(
      gene = rownames(expr_mat),
      lineage = lin,
      bin = as.integer(bn),
      pseudotime_mean = mean(cb$pseudotime[match(cs, cb$cellID)]),
      expr_mean = Matrix::rowMeans(expr_mat[, cs, drop = FALSE])
    )
  }))
}

direct_binned <- lapply(paste0("lineage", 1:3), make_gene_bins) |>
  bind_rows()

plaur_metrics <- direct_binned |>
  filter(.data$gene == "PLAUR") |>
  group_by(.data$gene, .data$lineage) |>
  summarise(
    lineage_label = lineage_labels[first(.data$lineage)],
    lineage_short = recode(first(.data$lineage), lineage1 = "L1", lineage2 = "L2", lineage3 = "L3"),
    early_mean = mean(.data$expr_mean[.data$bin <= 16], na.rm = TRUE),
    late_mean = mean(.data$expr_mean[.data$bin >= 65], na.rm = TRUE),
    mid_mean = mean(.data$expr_mean[.data$bin >= 33 & .data$bin <= 48], na.rm = TRUE),
    min_expr = min(.data$expr_mean, na.rm = TRUE),
    max_expr = max(.data$expr_mean, na.rm = TRUE),
    dynamic_range = max_expr - min_expr,
    early_late_delta = late_mean - early_mean,
    abs_delta = abs(early_late_delta),
    mean_expr = mean(.data$expr_mean, na.rm = TRUE),
    waldStat = NA_real_,
    p_adj_BH = NA_real_,
    pattern = "user_prior",
    .groups = "drop"
  ) |>
  mutate(
    timing = if_else(.data$early_late_delta < 0, "early_high", "late_high"),
    endpoint_balance = .data$abs_delta / pmax(.data$dynamic_range, 1e-6),
    visual_score = .data$dynamic_range * 0.62 + .data$abs_delta * 0.20 + .data$endpoint_balance * 0.05,
    row_id = paste(.data$lineage, .data$gene, sep = "|"),
    display_gene = .data$gene
  )

candidate_pool_aug <- bind_rows(
  candidate_pool,
  plaur_metrics |> select(names(candidate_pool))
) |>
  distinct(.data$row_id, .keep_all = TRUE)

binned_pool <- direct_binned |>
  inner_join(
    candidate_pool_aug |>
      select(
        gene, lineage, row_id, display_gene, lineage_label, lineage_short,
        timing, pattern, waldStat, p_adj_BH, early_mean, late_mean,
        dynamic_range, early_late_delta, visual_score
      ),
    by = c("gene", "lineage")
  ) |>
  group_by(row_id) |>
  mutate(
    z = as.numeric(scale(expr_mean)),
    z = pmax(pmin(z, 1.55), -1.55)
  ) |>
  ungroup()

smooth_row <- function(df, n_bins = 720) {
  df <- df[order(df$bin), ]
  xout <- seq(min(df$bin), max(df$bin), length.out = n_bins)
  fit <- tryCatch(
    stats::loess(z ~ bin, data = df, span = 0.28, degree = 1, family = "symmetric"),
    error = function(e) NULL
  )
  if (!is.null(fit)) {
    approx_values <- as.numeric(predict(fit, newdata = data.frame(bin = xout)))
  } else {
    approx_values <- approx(x = df$bin, y = df$z, xout = xout, rule = 2)$y
  }
  if (anyNA(approx_values)) {
    fallback <- approx(x = df$bin, y = df$z, xout = xout, rule = 2)$y
    approx_values[is.na(approx_values)] <- fallback[is.na(approx_values)]
  }
  approx_values <- pmax(pmin(approx_values, 1.55), -1.55)
  tibble(
    row_id = df$row_id[1],
    lineage = df$lineage[1],
    lineage_label = df$lineage_label[1],
    lineage_short = df$lineage_short[1],
    gene = df$gene[1],
    pattern = df$pattern[1],
    smooth_bin = paste0("T", sprintf("%02d", seq_len(n_bins))),
    z = approx_values
  )
}

smooth_pool <- binned_pool |>
  group_by(row_id) |>
  group_modify(~ smooth_row(.x)) |>
  ungroup()

row_peak <- smooth_pool |>
  group_by(.data$row_id) |>
  summarise(
    peak_smooth_bin = which.max(.data$z),
    peak_z = max(.data$z, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    peak_zone = dplyr::case_when(
      .data$peak_smooth_bin <= 240 ~ "left",
      .data$peak_smooth_bin <= 480 ~ "middle",
      TRUE ~ "right"
    )
  )

select_balanced <- function(df, n_each_zone = 10) {
  zone_levels <- c("left", "middle", "right")
  picked <- lapply(zone_levels, function(zone) {
    df |>
      filter(.data$peak_zone == zone) |>
      arrange(desc(.data$visual_score), desc(.data$dynamic_range)) |>
      slice_head(n = n_each_zone)
  }) |>
    bind_rows()
  if (nrow(picked) < n_each_zone * length(zone_levels)) {
    fill_n <- n_each_zone * length(zone_levels) - nrow(picked)
    picked <- bind_rows(
      picked,
      df |>
        filter(!.data$row_id %in% picked$row_id) |>
        arrange(desc(.data$visual_score), desc(.data$dynamic_range)) |>
        slice_head(n = fill_n)
    )
  }
  picked |>
    distinct(.data$row_id, .keep_all = TRUE)
}

selected <- candidate_pool |>
  left_join(row_peak, by = "row_id") |>
  group_by(.data$lineage, .data$lineage_label, .data$lineage_short) |>
  group_modify(~ select_balanced(.x, n_each_zone = 10)) |>
  ungroup()

plaur_keep <- plaur_metrics |>
  left_join(row_peak, by = "row_id") |>
  filter(.data$lineage %in% c("lineage1", "lineage2")) |>
  arrange(.data$lineage)
for (i in seq_len(nrow(plaur_keep))) {
  plaur_row <- plaur_keep[i, ]
  if (!plaur_row$row_id %in% selected$row_id) {
    replace_lin <- plaur_row$lineage[1]
    replace_zone <- plaur_row$peak_zone[1]
    drop_candidates <- selected |>
      filter(.data$lineage == replace_lin, .data$peak_zone == replace_zone, .data$display_gene != "PLAUR")
    if (nrow(drop_candidates) == 0) {
      drop_candidates <- selected |>
        filter(.data$lineage == replace_lin, .data$display_gene != "PLAUR")
    }
    drop_row <- drop_candidates |>
      arrange(.data$visual_score) |>
      slice_head(n = 1) |>
      pull(.data$row_id)
    selected <- selected |>
      filter(!.data$row_id %in% drop_row) |>
      bind_rows(plaur_row |> select(names(selected)))
    log_msg("Forced PLAUR into", plaur_row$lineage_short[1], "zone", replace_zone, "replacing", drop_row)
  }
}

zone_counts <- selected |>
  count(.data$lineage_short, .data$peak_zone)
readr::write_csv(
  zone_counts,
  "figures/final_panels/source_data/supplement/S7A_peak_zone_balance.csv"
)

binned <- binned_pool |>
  semi_join(selected |> select(row_id), by = "row_id")
smooth_binned <- smooth_pool |>
  semi_join(selected |> select(row_id), by = "row_id")

row_order <- selected |>
  arrange(
    factor(.data$lineage, levels = c("lineage1", "lineage2", "lineage3")),
    factor(.data$peak_zone, levels = c("left", "middle", "right")),
    .data$peak_smooth_bin,
    desc(.data$dynamic_range),
    desc(.data$visual_score)
  ) |>
  pull(.data$row_id)

mat <- smooth_binned |>
  select(row_id, smooth_bin, z) |>
  pivot_wider(names_from = smooth_bin, values_from = z) |>
  tibble::column_to_rownames("row_id") |>
  as.matrix()
mat <- mat[row_order, , drop = FALSE]

row_meta <- selected |>
  mutate(row_id = factor(row_id, levels = row_order)) |>
  arrange(row_id) |>
  select(
    row_id, display_gene, lineage_label, lineage_short, timing,
    pattern, early_mean, late_mean, dynamic_range, early_late_delta,
    visual_score, peak_zone, peak_smooth_bin, peak_z, waldStat, p_adj_BH
  )

rownames(mat) <- row_meta$display_gene
row_split <- data.frame(
  Lineage = factor(row_meta$lineage_short, levels = c("L1", "L2", "L3"))
)

lineage_short_colors <- c("L1" = "#BC3C29", "L2" = "#20854E", "L3" = "#E18727")
left_anno <- rowAnnotation(
  Lineage = row_meta$lineage_short,
  col = list(Lineage = lineage_short_colors),
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 7),
  annotation_legend_param = list(
    Lineage = list(
      title = "Lineage",
      at = c("L1", "L2", "L3"),
      labels = c("L1: NPC-P to MES-I", "L2: NPC-P to MES-V", "L3: NPC-P to OPC-M"),
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7)
    )
  )
)

col_fun <- colorRamp2(c(-1.4, 0, 1.4), c("#3B4CC0", "white", "#B40426"))
ht <- Heatmap(
  mat,
  name = "Z",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_split,
  cluster_row_slices = FALSE,
  left_annotation = left_anno,
  row_names_gp = gpar(fontsize = 4.6),
  column_names_gp = gpar(fontsize = 0),
  show_column_names = FALSE,
  column_title = "S7A. Driver gene dynamics across three inferred lineages",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  row_title_gp = gpar(fontsize = 7, fontface = "bold"),
  heatmap_legend_param = list(title = "Scaled expression", title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7)),
  border = FALSE,
  use_raster = TRUE
)

pdf_path <- "figures/final_panels/supplement/S7A_three_lineage_driver_heatmap.pdf"
png_path <- "figures/final_panels/supplement/S7A_three_lineage_driver_heatmap.png"
pdf(pdf_path, width = 7.4, height = 10.8, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", padding = unit(c(3, 3, 3, 3), "mm"))
dev.off()

png(png_path, width = 7.4, height = 10.8, units = "in", res = 400, type = "cairo")
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", padding = unit(c(3, 3, 3, 3), "mm"))
dev.off()

write_csv(
  binned |> select(gene, lineage, lineage_label, lineage_short, pattern, bin, pseudotime_mean, expr_mean, z, waldStat, p_adj_BH),
  "figures/final_panels/source_data/supplement/S7A_three_lineage_driver_heatmap_original_bins.csv"
)
write_csv(
  smooth_binned |> left_join(selected |> select(row_id, waldStat, p_adj_BH), by = "row_id"),
  "figures/final_panels/source_data/supplement/S7A_three_lineage_driver_heatmap_smooth_bins.csv"
)
write_csv(
  row_meta,
  "figures/final_panels/source_data/supplement/S7A_three_lineage_driver_genes_selected.csv"
)

log_msg("Selected genes:", nrow(row_meta))
log_msg("Wrote", pdf_path)
log_msg("Wrote", png_path)
log_msg("STOP")
