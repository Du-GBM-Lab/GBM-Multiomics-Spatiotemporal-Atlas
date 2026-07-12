suppressPackageStartupMessages({
  .libPaths(c("<R_LIBS>", "<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(qs2)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(SingleCellExperiment)
  library(slingshot)
  library(patchwork)
})

# Phase 2B: multi-root Slingshot sensitivity + CytoTRACE2/pseudotime consistency.
# Read-only for source objects: no UMAP/PCA/cluster rerun and no qs2 writeback.

started_at <- Sys.time()
set.seed(42)

cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(cwd) == "scripts") {
  step_dir <- normalizePath(file.path(cwd, ".."), winslash = "/", mustWork = TRUE)
} else if (file.exists(file.path(cwd, "scripts", "_naming.R"))) {
  step_dir <- cwd
} else {
  step_dir <- normalizePath(file.path(cwd, "06_恶性细胞拟时序"), winslash = "/", mustWork = TRUE)
}
setwd(step_dir)

dir.create("tables", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/source_data", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

source(file.path("scripts", "_naming.R"))

log_file <- file.path("logs", "04_multi_root_sensitivity.log")
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_msg("Phase 2B multi-root sensitivity started")
log_msg("set.seed = 42")
log_msg("slingshot", as.character(packageVersion("slingshot")))
log_msg("SingleCellExperiment", as.character(packageVersion("SingleCellExperiment")))

obj_path <- normalizePath(file.path("..", "05_恶性细胞分亚群与Neftel对照", "outputs", "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"), winslash = "/", mustWork = TRUE)
coord_obj_path <- normalizePath(file.path("..", "05_恶性细胞分亚群与Neftel对照", "outputs", "GBM.malignant.subtyped.umap_candidates.qs2"), winslash = "/", mustWork = TRUE)
cytotrace_path <- file.path("outputs", "CytoTRACE2_per_cell.csv")
if (!file.exists(cytotrace_path)) stop("Missing cached CytoTRACE2 file: ", cytotrace_path)

log_msg("Input object:", obj_path)
log_msg("Coordinate object:", coord_obj_path)
log_msg("Coordinate reduction: umap_closer4")
log_msg("Cached CytoTRACE2:", cytotrace_path)

load_start <- Sys.time()
obj <- qs2::qs_read(obj_path)
coord_obj <- qs2::qs_read(coord_obj_path)
log_msg("Object load seconds:", round(as.numeric(difftime(Sys.time(), load_start, units = "secs")), 2))

md <- obj@meta.data |>
  rownames_to_column("cellID") |>
  recode_subtype()
cytotrace <- readr::read_csv(cytotrace_path, show_col_types = FALSE)
if (!all(c("cell_id", "CytoTRACE2_score") %in% colnames(cytotrace))) {
  stop("Cached CytoTRACE2 file lacks required columns.")
}
md <- md |>
  dplyr::left_join(cytotrace, by = c("cellID" = "cell_id"))

umap <- Seurat::Embeddings(coord_obj, "umap_closer4")
if (!setequal(rownames(umap), md$cellID)) stop("UMAP coordinate cell IDs do not match metadata.")
umap <- umap[md$cellID, 1:2, drop = FALSE]
colnames(umap) <- c("UMAP_1", "UMAP_2")

cell_df <- as.data.frame(umap) |>
  rownames_to_column("cellID") |>
  dplyr::left_join(
    md |> dplyr::select(cellID, subtype_k4, subtype_label_final, subtype_short, Pt_number, CytoTRACE2_score),
    by = "cellID"
  ) |>
  dplyr::mutate(subtype_short = factor(subtype_short, levels = subtype_naming_mapping$abbreviation))
if (nrow(cell_df) != 28213) stop("Unexpected cell count: ", nrow(cell_df))

short_levels <- subtype_naming_mapping$abbreviation
short_colors <- setNames(subtype_naming_mapping$color, subtype_naming_mapping$abbreviation)

make_base_sce <- function() {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(dummy = Matrix::Matrix(0, nrow = 1, ncol = nrow(cell_df), sparse = TRUE)),
    colData = md |> dplyr::select(cellID, subtype_k4, subtype_label_final, subtype_short)
  )
  colnames(sce) <- cell_df$cellID
  SingleCellExperiment::reducedDim(sce, "UMAP") <- as.matrix(cell_df[, c("UMAP_1", "UMAP_2")])
  sce
}

run_configs <- tibble(
  run_id = c("run1_NPCP", "run2_OPCM", "run3_MESV", "run4_auto"),
  root_cluster = c("NPC-P", "OPC-M", "MES-V", NA_character_),
  panel_label = c("Root: NPC-P", "Root: OPC-M", "Root: MES-V", "Root: auto")
)

run_slingshot_once <- function(run_id, root_cluster) {
  log_msg("Slingshot run started:", run_id, "root =", ifelse(is.na(root_cluster), "NULL", root_cluster))
  sce <- make_base_sce()
  cluster_labels <- factor(as.character(cell_df$subtype_short), levels = short_levels)
  run_start <- Sys.time()
  if (is.na(root_cluster)) {
    sce <- slingshot::slingshot(
      sce,
      clusterLabels = cluster_labels,
      reducedDim = "UMAP",
      start.clus = NULL,
      end.clus = NULL,
      approx_points = 150
    )
  } else {
    sce <- slingshot::slingshot(
      sce,
      clusterLabels = cluster_labels,
      reducedDim = "UMAP",
      start.clus = root_cluster,
      end.clus = NULL,
      approx_points = 150
    )
  }
  log_msg("Slingshot seconds:", run_id, round(as.numeric(difftime(Sys.time(), run_start, units = "secs")), 2))

  lineages <- slingshot::slingLineages(sce)
  pst <- as.data.frame(slingshot::slingPseudotime(sce))
  weights <- as.data.frame(slingshot::slingCurveWeights(sce))
  colnames(pst) <- paste0("lineage", seq_len(ncol(pst)), "_pst")
  colnames(weights) <- paste0("lineage", seq_len(ncol(weights)), "_weight")
  weight_mat <- as.matrix(weights)
  assignment <- if (ncol(weight_mat) > 0) paste0("lineage", max.col(weight_mat, ties.method = "first")) else rep(NA_character_, nrow(weight_mat))
  assignment[rowSums(weight_mat, na.rm = TRUE) == 0] <- NA_character_
  assigned_weight <- if (ncol(weight_mat) > 0) apply(weight_mat, 1, max, na.rm = TRUE) else rep(NA_real_, nrow(weight_mat))
  assigned_pst <- rep(NA_real_, nrow(pst))
  for (i in seq_len(nrow(pst))) {
    if (!is.na(assignment[i])) {
      pst_col <- paste0(assignment[i], "_pst")
      assigned_pst[i] <- pst[[pst_col]][i]
    }
  }

  per_cell <- bind_cols(
    cell_df |> dplyr::select(cellID, subtype_k4, subtype_label_final, subtype_short, Pt_number, CytoTRACE2_score),
    pst,
    weights,
    tibble(
      run_id = run_id,
      lineage_assignment = assignment,
      assigned_weight = assigned_weight,
      assigned_pseudotime = assigned_pst
    )
  )

  curve_coord <- lapply(seq_along(slingshot::slingCurves(sce)), function(i) {
    crv <- slingshot::slingCurves(sce)[[i]]
    tibble(
      run_id = run_id,
      lineage_id = paste0("lineage", i),
      point_id = seq_len(nrow(crv$s)),
      UMAP_1 = crv$s[, 1],
      UMAP_2 = crv$s[, 2],
      curve_order = seq_len(nrow(crv$s))
    )
  }) |>
    bind_rows()

  lineage_meta <- tibble(
    run_id = run_id,
    lineage_id = paste0("lineage", seq_along(lineages)),
    start_cluster = vapply(lineages, function(x) x[1], character(1)),
    end_cluster = vapply(lineages, function(x) x[length(x)], character(1)),
    lineage_path = vapply(lineages, paste, character(1), collapse = " -> ")
  )

  lineage_counts <- per_cell |>
    dplyr::filter(!is.na(lineage_assignment), assigned_weight > 0.5) |>
    dplyr::count(run_id, lineage_assignment, name = "n_cells_assigned") |>
    dplyr::rename(lineage_id = lineage_assignment)

  list(
    sce = sce,
    per_cell = per_cell,
    curve_coord = curve_coord,
    lineage_meta = lineage_meta |> dplyr::left_join(lineage_counts, by = c("run_id", "lineage_id")),
    n_lineages = length(lineages)
  )
}

results <- lapply(seq_len(nrow(run_configs)), function(i) {
  run_slingshot_once(run_configs$run_id[i], run_configs$root_cluster[i])
})
names(results) <- run_configs$run_id

per_cell_all <- bind_rows(lapply(results, `[[`, "per_cell"))
curve_all <- bind_rows(lapply(results, `[[`, "curve_coord"))
lineage_meta_all <- bind_rows(lapply(results, `[[`, "lineage_meta")) |>
  dplyr::left_join(run_configs, by = "run_id")

readr::write_csv(per_cell_all, file.path("tables", "multiroot_pseudotime_per_cell.csv"))
readr::write_csv(curve_all, file.path("tables", "multiroot_curve_coordinates.csv"))

lineage_summary <- lineage_meta_all |>
  dplyr::group_by(run_id, root_cluster, panel_label) |>
  dplyr::summarise(
    n_lineages = dplyr::n(),
    terminal_clusters = paste(end_cluster, collapse = ";"),
    lineage_paths = paste(paste0(lineage_id, ": ", lineage_path), collapse = " | "),
    n_cells_per_lineage = paste(paste0(lineage_id, "=", ifelse(is.na(n_cells_assigned), 0L, n_cells_assigned)), collapse = ";"),
    .groups = "drop"
  ) |>
  dplyr::arrange(run_id)
readr::write_csv(lineage_summary, file.path("tables", "multiroot_lineage_summary.csv"))

wide_pt <- per_cell_all |>
  dplyr::select(cellID, run_id, assigned_pseudotime) |>
  tidyr::pivot_wider(names_from = run_id, values_from = assigned_pseudotime)
run_ids <- run_configs$run_id
pt_corr <- tidyr::expand_grid(run_x = run_ids, run_y = run_ids) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    n_cells = sum(is.finite(wide_pt[[run_x]]) & is.finite(wide_pt[[run_y]])),
    spearman_rho = ifelse(
      n_cells >= 3,
      suppressWarnings(cor(wide_pt[[run_x]], wide_pt[[run_y]], method = "spearman", use = "pairwise.complete.obs")),
      NA_real_
    )
  ) |>
  dplyr::ungroup()
readr::write_csv(pt_corr, file.path("tables", "multiroot_pseudotime_correlation.csv"))

run1 <- results[["run1_NPCP"]]$per_cell
cyto_corr <- lapply(sort(unique(run1$lineage_assignment)), function(lin) {
  df <- run1 |>
    dplyr::filter(lineage_assignment == lin, assigned_weight > 0.5, !is.na(CytoTRACE2_score), is.finite(assigned_pseudotime))
  if (nrow(df) < 3) {
    return(tibble(lineage_id = lin, n_cells = nrow(df), spearman_rho = NA_real_, p_value = NA_real_, expected_direction = "negative", flag = "too_few_cells"))
  }
  ct <- suppressWarnings(cor.test(df$assigned_pseudotime, df$CytoTRACE2_score, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  tibble(
    lineage_id = lin,
    n_cells = nrow(df),
    spearman_rho = rho,
    p_value = ct$p.value,
    expected_direction = "negative",
    flag = dplyr::case_when(
      is.na(rho) ~ "NA",
      rho < 0 ~ "stemness_decreases",
      TRUE ~ "WARNING_positive_or_flat"
    )
  )
}) |>
  bind_rows()
readr::write_csv(cyto_corr, file.path("tables", "cytotrace2_pseudotime_correlation.csv"))

l3_robustness <- lineage_meta_all |>
  dplyr::group_by(run_id, root_cluster, panel_label) |>
  dplyr::summarise(
    OPC_M_reached_as_terminal = any(end_cluster == "OPC-M"),
    OPC_M_lineages = paste(lineage_id[end_cluster == "OPC-M"], collapse = ";"),
    n_cells_on_OPC_M_lineage = sum(ifelse(end_cluster == "OPC-M" & !is.na(n_cells_assigned), n_cells_assigned, 0L)),
    MES_I_terminal_present = any(end_cluster == "MES-I"),
    MES_I_lineages = paste(lineage_id[end_cluster == "MES-I"], collapse = ";"),
    n_cells_on_MES_I_lineage = sum(ifelse(end_cluster == "MES-I" & !is.na(n_cells_assigned), n_cells_assigned, 0L)),
    .groups = "drop"
  ) |>
  dplyr::arrange(run_id)
readr::write_csv(l3_robustness, file.path("tables", "L3_robustness_summary.csv"))

readr::write_csv(lineage_meta_all, file.path("figures/source_data", "04_lineage_meta.csv"))
readr::write_csv(curve_all, file.path("figures/source_data", "04_curve_coordinates.csv"))
readr::write_csv(pt_corr, file.path("figures/source_data", "04_panel_E_correlation_heatmap.csv"))
readr::write_csv(cyto_corr, file.path("figures/source_data", "05_cytotrace2_correlation.csv"))

base_theme <- theme_bw(base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    plot.title = element_text(size = 8),
    plot.subtitle = element_text(size = 6.2)
  )

plot_one_run <- function(run_id) {
  cfg <- run_configs |> dplyr::filter(run_id == .env$run_id)
  curves <- curve_all |> dplyr::filter(run_id == .env$run_id)
  labels <- curves |>
    dplyr::group_by(lineage_id) |>
    dplyr::slice_tail(n = 1) |>
    dplyr::ungroup()
  arrows <- curves |>
    dplyr::group_by(lineage_id) |>
    dplyr::slice_tail(n = 2) |>
    dplyr::mutate(row = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::select(lineage_id, row, UMAP_1, UMAP_2) |>
    tidyr::pivot_wider(id_cols = lineage_id, names_from = row, values_from = c(UMAP_1, UMAP_2))
  roots <- if (is.na(cfg$root_cluster)) {
    tibble(UMAP_1 = numeric(), UMAP_2 = numeric())
  } else {
    cell_df |>
      dplyr::filter(subtype_short == cfg$root_cluster) |>
      dplyr::summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
  }
  terminals <- lineage_meta_all |>
    dplyr::filter(run_id == .env$run_id) |>
    dplyr::select(end_cluster) |>
    dplyr::distinct() |>
    dplyr::left_join(
      cell_df |>
        dplyr::group_by(subtype_short) |>
        dplyr::summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop"),
      by = c("end_cluster" = "subtype_short")
    )
  ggplot(cell_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = subtype_short), size = 0.06, alpha = 0.18) +
    geom_path(data = curves, aes(UMAP_1, UMAP_2, group = lineage_id), inherit.aes = FALSE, color = "black", linewidth = 0.42, lineend = "round") +
    geom_segment(
      data = arrows,
      aes(x = UMAP_1_1, y = UMAP_2_1, xend = UMAP_1_2, yend = UMAP_2_2),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.07, "inches")),
      color = "black",
      linewidth = 0.42
    ) +
    geom_point(data = roots, aes(UMAP_1, UMAP_2), inherit.aes = FALSE, shape = 21, fill = "white", color = "black", size = 1.8, stroke = 0.45) +
    geom_point(data = terminals, aes(UMAP_1, UMAP_2), inherit.aes = FALSE, shape = 24, fill = "black", color = "black", size = 1.6) +
    geom_text(data = labels, aes(UMAP_1, UMAP_2, label = gsub("lineage", "L", lineage_id)), inherit.aes = FALSE, size = 2.1, vjust = -0.4) +
    scale_color_manual(values = short_colors, guide = "none") +
    coord_equal() +
    labs(title = cfg$panel_label, subtitle = paste0("Terminals: ", paste(unique(lineage_meta_all$end_cluster[lineage_meta_all$run_id == run_id]), collapse = ", ")), x = "UMAP 1", y = "UMAP 2") +
    base_theme
}

p_runs <- lapply(run_ids, plot_one_run)
names(p_runs) <- run_ids

p_heat <- ggplot(pt_corr, aes(run_x, run_y, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 2.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  coord_equal() +
  labs(title = "E. Cross-run pseudotime correlation", x = NULL, y = NULL) +
  theme_bw(base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    plot.title = element_text(size = 8)
  )

p_root <- (p_runs[[1]] | p_runs[[2]]) / (p_runs[[3]] | p_runs[[4]]) / p_heat +
  patchwork::plot_layout(heights = c(1, 1, 0.82)) +
  patchwork::plot_annotation(
    caption = "All runs use the same malignant-only umap_closer4 coordinates and subtype labels. White circle marks specified root; black triangles mark terminal subtype centroids.",
    theme = theme(plot.caption = element_text(size = 6, hjust = 0))
  )
pdf(file.path("figures", "04_root_sensitivity.pdf"), width = 9.0, height = 10.5, useDingbats = FALSE)
print(p_root)
dev.off()
ggsave(file.path("figures", "04_root_sensitivity_preview.png"), p_root, width = 9.0, height = 10.5, dpi = 180, bg = "white")
log_msg("Wrote root sensitivity figure:", file.path("figures", "04_root_sensitivity.pdf"))

run1_plot <- run1 |>
  dplyr::filter(!is.na(lineage_assignment), assigned_weight > 0.5, !is.na(CytoTRACE2_score), is.finite(assigned_pseudotime)) |>
  dplyr::mutate(lineage_label = factor(lineage_assignment, levels = paste0("lineage", 1:3), labels = paste0("L", 1:3)))
cyto_labels <- cyto_corr |>
  dplyr::mutate(
    lineage_label = factor(lineage_id, levels = paste0("lineage", 1:3), labels = paste0("L", 1:3)),
    label = paste0("rho = ", sprintf("%.2f", spearman_rho), "\np = ", format.pval(p_value, digits = 2, eps = 1e-300), "\nn = ", n_cells)
  )
readr::write_csv(run1_plot, file.path("figures/source_data", "05_cytotrace2_pseudotime_cells.csv"))

p_cyto <- ggplot(run1_plot, aes(assigned_pseudotime, CytoTRACE2_score, color = subtype_short)) +
  geom_point(size = 0.12, alpha = 0.28) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, linewidth = 0.55, color = "black", span = 0.65) +
  facet_wrap(~ lineage_label, scales = "free_x", nrow = 1) +
  geom_text(data = cyto_labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.05, vjust = 1.05, size = 2.1) +
  scale_color_manual(values = short_colors, name = NULL) +
  labs(
    title = "CytoTRACE2 vs primary pseudotime by lineage",
    subtitle = "Monotonic stemness loss is strongest for the OPC-M lineage (L3); MES-like lineages show non-monotonic ordering.",
    x = "Primary-run lineage pseudotime",
    y = "CytoTRACE2 score",
    caption = "L1/L2 trends should be interpreted cautiously because MES-like fine ordering is method-sensitive; monotonic stemness-loss claims are restricted to L3."
  ) +
  theme_bw(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 8),
    plot.subtitle = element_text(size = 6.4),
    plot.caption = element_text(size = 6, hjust = 0),
    strip.text = element_text(size = 7)
  )
pdf(file.path("figures", "05_cytotrace2_pseudotime.pdf"), width = 8.8, height = 3.4, useDingbats = FALSE)
print(p_cyto)
dev.off()
ggsave(file.path("figures", "05_cytotrace2_pseudotime_preview.png"), p_cyto, width = 8.8, height = 3.4, dpi = 180, bg = "white")
log_msg("Wrote CytoTRACE2 figure:", file.path("figures", "05_cytotrace2_pseudotime.pdf"))

log_msg("Lineage summary:", file.path("tables", "multiroot_lineage_summary.csv"))
log_msg("Pseudotime correlation:", file.path("tables", "multiroot_pseudotime_correlation.csv"))
log_msg("CytoTRACE2 correlation:", file.path("tables", "cytotrace2_pseudotime_correlation.csv"))
log_msg("L3 robustness:", file.path("tables", "L3_robustness_summary.csv"))
capture.output(print(lineage_summary)) |> paste(collapse = "\n") |> log_msg()
capture.output(print(cyto_corr)) |> paste(collapse = "\n") |> log_msg()

session_info <- capture.output(sessionInfo())
cat(session_info, sep = "\n", file = file.path("logs", "04_session_info.txt"))
log_msg("Session info:", file.path("logs", "04_session_info.txt"))
log_msg("Total elapsed seconds:", round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 2))
log_msg("STOP: Phase 2B complete. Await audit before Phase 3.")
