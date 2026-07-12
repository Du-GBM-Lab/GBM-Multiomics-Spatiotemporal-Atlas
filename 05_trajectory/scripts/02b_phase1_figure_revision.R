suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggridges)
  library(patchwork)
})

# Phase 1 figure revision.
# Read-only with respect to Slingshot: this script uses existing pseudotime,
# lineage weight, curve coordinate, UMAP, and connectivity CSV files only.

started_at <- Sys.time()

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

log_file <- file.path("logs", "02b_phase1_figure_revision.log")
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_msg("Revision script started")
log_msg("No Slingshot rerun; reading existing Phase 1 CSV outputs")

source_a <- readr::read_csv(file.path("figures/source_data", "02_panel_A.csv"), show_col_types = FALSE)
pseudotime <- readr::read_csv(file.path("tables", "slingshot_pseudotime_per_cell.csv"), show_col_types = FALSE)
curve_coord <- readr::read_csv(file.path("tables", "slingshot_curve_coordinates.csv"), show_col_types = FALSE)
connectivity <- readr::read_csv(file.path("tables", "paga_connectivity_matrix.csv"), show_col_types = FALSE)

short_colors <- setNames(subtype_naming_mapping$color, subtype_naming_mapping$abbreviation)
short_levels <- subtype_naming_mapping$abbreviation
lineage_ids <- sort(unique(curve_coord$lineage_id))
weight_threshold <- 0.5

required_cols <- c("cellID", "UMAP_1", "UMAP_2", "subtype_k4", "subtype_label_final", "subtype_short")
if (!all(required_cols %in% colnames(source_a))) {
  stop("02_panel_A.csv lacks required coordinate/subtype columns.")
}
if (!all(paste0(lineage_ids, "_pst") %in% colnames(pseudotime))) {
  stop("Pseudotime table lacks one or more lineage pseudotime columns.")
}
if (!all(paste0(lineage_ids, "_weight") %in% colnames(pseudotime))) {
  stop("Pseudotime table lacks one or more lineage weight columns.")
}

cell_df <- source_a |>
  left_join(pseudotime, by = c("cellID", "subtype_k4", "subtype_label_final", "subtype_short")) |>
  mutate(subtype_short = factor(subtype_short, levels = short_levels))
if (nrow(cell_df) != 28213) stop("Unexpected cell count after join: ", nrow(cell_df))

lineage_weight <- pseudotime |>
  select(cellID, subtype_k4, subtype_label_final, subtype_short, Pt_number, lineage_assignment, ends_with("_weight")) |>
  pivot_longer(
    cols = ends_with("_weight"),
    names_to = "lineage_id",
    values_to = "lineage_weight"
  ) |>
  mutate(
    lineage_id = sub("_weight$", "", lineage_id),
    assigned_by_weight = lineage_assignment == lineage_id & lineage_weight > weight_threshold
  )
readr::write_csv(lineage_weight, file.path("tables", "slingshot_lineage_weight_per_cell.csv"))

lineage_plot_long <- lineage_ids |>
  lapply(function(lin) {
    pst_col <- paste0(lin, "_pst")
    weight_col <- paste0(lin, "_weight")
    cell_df |>
      transmute(
        cellID,
        UMAP_1,
        UMAP_2,
        subtype_k4,
        subtype_label_final,
        subtype_short,
        lineage_id = lin,
        pseudotime = .data[[pst_col]],
        lineage_weight = .data[[weight_col]],
        assigned_by_weight = lineage_assignment == lin & .data[[weight_col]] > weight_threshold
      )
  }) |>
  bind_rows()

readr::write_csv(lineage_plot_long, file.path("figures/source_data", "02_panel_C.csv"))

ridge_source <- lineage_plot_long |>
  filter(assigned_by_weight, !is.na(pseudotime)) |>
  mutate(
    subtype_short = factor(subtype_short, levels = rev(short_levels)),
    lineage_panel = paste0(gsub("lineage", "L", lineage_id), " assigned cells")
  )
readr::write_csv(ridge_source, file.path("figures/source_data", "02_panel_E.csv"))

lineage_subtype_crosstab <- ridge_source |>
  count(lineage_id, subtype_k4, subtype_label_final, subtype_short, name = "n_cells") |>
  group_by(lineage_id) |>
  mutate(lineage_fraction = n_cells / sum(n_cells)) |>
  ungroup() |>
  arrange(lineage_id, subtype_short)
readr::write_csv(lineage_subtype_crosstab, file.path("tables", "lineage_subtype_crosstab.csv"))

branch_points <- curve_coord |>
  group_by(lineage_id) |>
  slice_min(curve_order, n = 1, with_ties = FALSE) |>
  ungroup()
branch_point_summary <- tibble(
  branch_point_id = "BP1",
  UMAP_1 = mean(branch_points$UMAP_1),
  UMAP_2 = mean(branch_points$UMAP_2),
  involved_lineages = paste(lineage_ids, collapse = ";"),
  estimation_rule = "centroid_of_first_curve_point_per_lineage",
  note = "Single trifurcation inferred from NPC-P; no further bifurcation detected. BP1 is a curve-start centroid, not an experimentally observed fate-branching event."
)
readr::write_csv(branch_point_summary, file.path("tables", "slingshot_branch_points.csv"))
readr::write_csv(
  branch_point_summary,
  file.path("figures/source_data", "02_panel_B_branch_points.csv")
)

curve_arrows <- curve_coord |>
  group_by(lineage_id) |>
  slice_tail(n = 2) |>
  mutate(row = row_number()) |>
  ungroup() |>
  select(lineage_id, row, UMAP_1, UMAP_2) |>
  pivot_wider(id_cols = lineage_id, names_from = row, values_from = c(UMAP_1, UMAP_2))
curve_labels <- curve_coord |>
  group_by(lineage_id) |>
  slice_tail(n = 1) |>
  ungroup()

node_df <- source_a |>
  count(subtype_short, subtype_label_final, name = "n_cells") |>
  mutate(
    angle = seq(0, 2 * pi, length.out = n() + 1)[-(n() + 1)],
    x = cos(angle),
    y = sin(angle),
    node_size = log1p(n_cells)
  )
edge_df <- connectivity |>
  filter(from_cluster < to_cluster, connectivity >= 0.1) |>
  left_join(node_df |> select(from_cluster = subtype_short, x_from = x, y_from = y), by = "from_cluster") |>
  left_join(node_df |> select(to_cluster = subtype_short, x_to = x, y_to = y), by = "to_cluster")

base_theme <- theme_bw(base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    plot.title = element_text(size = 8),
    plot.subtitle = element_text(size = 6.5)
  )

p_a <- ggplot(source_a, aes(UMAP_1, UMAP_2, color = subtype_short)) +
  geom_point(size = 0.12, alpha = 0.65) +
  scale_color_manual(values = short_colors, name = NULL) +
  coord_equal() +
  labs(title = "A. Malignant subtype UMAP", x = "UMAP 1", y = "UMAP 2") +
  base_theme

p_b <- ggplot(source_a, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = subtype_short), size = 0.08, alpha = 0.22) +
  geom_path(data = curve_coord, aes(group = lineage_id), color = "black", linewidth = 0.45, lineend = "round") +
  geom_segment(
    data = curve_arrows,
    aes(x = UMAP_1_1, y = UMAP_2_1, xend = UMAP_1_2, yend = UMAP_2_2),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.08, "inches")),
    linewidth = 0.45,
    color = "black"
  ) +
  geom_point(
    data = branch_point_summary,
    aes(UMAP_1, UMAP_2),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    color = "black",
    size = 2.1,
    stroke = 0.5
  ) +
  geom_text(data = branch_point_summary, aes(UMAP_1, UMAP_2, label = branch_point_id), inherit.aes = FALSE, size = 2.1, vjust = 1.6) +
  geom_text(data = curve_labels, aes(label = gsub("lineage", "L", lineage_id)), size = 2.2, vjust = -0.4, color = "black") +
  scale_color_manual(values = short_colors, guide = "none") +
  coord_equal() +
  labs(title = "B. Slingshot curves", subtitle = "Root: NPC-P; single trifurcation, no further bifurcation detected", x = "UMAP 1", y = "UMAP 2") +
  base_theme

p_c <- ggplot(lineage_plot_long, aes(UMAP_1, UMAP_2)) +
  geom_point(color = "grey88", size = 0.08, alpha = 0.45) +
  geom_point(
    data = lineage_plot_long |> filter(assigned_by_weight, !is.na(pseudotime)),
    aes(color = pseudotime),
    size = 0.12,
    alpha = 0.78
  ) +
  facet_wrap(~ lineage_id, nrow = 1, labeller = as_labeller(function(x) gsub("lineage", "L", x))) +
  scale_color_viridis_c(option = "viridis", name = "Pseudotime", na.value = "grey88") +
  coord_equal() +
  labs(title = "C. Lineage-specific pseudotime", subtitle = "Colored cells: assigned lineage and weight > 0.5", x = "UMAP 1", y = "UMAP 2") +
  base_theme +
  theme(strip.text = element_text(size = 7), legend.position = "right")

p_d <- ggplot() +
  geom_segment(data = edge_df, aes(x = x_from, y = y_from, xend = x_to, yend = y_to, linewidth = connectivity), color = "grey35", alpha = 0.75) +
  geom_point(data = node_df, aes(x, y, fill = subtype_short, size = node_size), shape = 21, color = "black", stroke = 0.35) +
  geom_text(data = node_df, aes(x, y, label = subtype_short), size = 2.4, vjust = -1.0) +
  scale_fill_manual(values = short_colors, guide = "none") +
  scale_size_continuous(range = c(4, 8), guide = "none") +
  scale_linewidth_continuous(range = c(0.3, 2.0), name = "Connectivity") +
  coord_equal(xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35)) +
  labs(title = "D. kNN-based connectivity", subtitle = "Edges with connectivity < 0.1 not shown", x = NULL, y = NULL) +
  base_theme

p_e <- ggplot(ridge_source, aes(x = pseudotime, y = subtype_short, fill = subtype_short)) +
  ggridges::geom_density_ridges(scale = 1.05, alpha = 0.72, linewidth = 0.25, color = "white") +
  facet_wrap(~ lineage_panel, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = short_colors, guide = "none") +
  labs(title = "E. Pseudotime density by lineage", subtitle = "Only cells assigned to each lineage with weight > 0.5 shown", x = "Lineage-specific pseudotime", y = NULL) +
  theme_bw(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 8),
    plot.subtitle = element_text(size = 6.5),
    axis.text = element_text(size = 6),
    strip.text = element_text(size = 7)
  )

caption <- paste(
  "Pseudotime represents an inferred transcriptional ordering and does not imply real-time progression or deterministic lineage commitment.",
  "Lineage-specific UMAP and ridge panels show cells assigned to each lineage with lineage weight > 0.5 only; other cells are grey or omitted."
)
p_all <- (p_a | p_b) / p_c / (p_d | p_e) +
  patchwork::plot_layout(heights = c(1, 0.82, 1.2)) +
  patchwork::plot_annotation(caption = caption, theme = theme(plot.caption = element_text(size = 6, hjust = 0)))

pdf(file.path("figures", "02_trajectory_primary.pdf"), width = 10.2, height = 11.2, useDingbats = FALSE)
print(p_all)
dev.off()
ggsave(file.path("figures", "02_trajectory_primary_preview.png"), p_all, width = 10.2, height = 11.2, dpi = 180, bg = "white")

log_msg("Wrote revised figure:", file.path("figures", "02_trajectory_primary.pdf"))
log_msg("Wrote preview:", file.path("figures", "02_trajectory_primary_preview.png"))
log_msg("Wrote branch points:", file.path("tables", "slingshot_branch_points.csv"))
log_msg("Wrote lineage weights:", file.path("tables", "slingshot_lineage_weight_per_cell.csv"))
log_msg("Wrote lineage-subtype crosstab:", file.path("tables", "lineage_subtype_crosstab.csv"))
log_msg("Lineages:", paste(lineage_ids, collapse = ", "))
log_msg("Weight threshold:", weight_threshold)
log_msg("Branch structure: single trifurcation from NPC-P; no further bifurcation detected")
log_msg("Total elapsed seconds:", round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 2))
log_msg("STOP: Phase 1 figure revision complete. Await audit before Phase 2.")
