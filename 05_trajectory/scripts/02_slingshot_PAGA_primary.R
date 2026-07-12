suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(Matrix)
})

# Phase 1 primary trajectory: Slingshot + PAGA/kNN connectivity.
# Read-only: does not modify or save the Seurat object.

started_at <- Sys.time()
set.seed(42)

cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(cwd) == "scripts") {
  step_dir <- normalizePath(file.path(cwd, ".."), winslash = "/", mustWork = TRUE)
} else if (file.exists(file.path(cwd, "scripts", "_naming.R"))) {
  step_dir <- cwd
} else {
  candidate <- file.path(cwd, "06_恶性细胞拟时序")
  if (!dir.exists(candidate)) candidate <- file.path(cwd, "06_鎭舵€х粏鑳炴嫙鏃跺簭")
  step_dir <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
}
setwd(step_dir)

dir.create("tables", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/source_data", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

source(file.path("scripts", "_naming.R"))

log_file <- file.path("logs", "02_slingshot_PAGA.log")
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

install_cran_if_needed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    log_msg("Installing CRAN package:", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

install_bioc_if_needed <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    log_msg("Installing Bioconductor package:", pkg)
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

install_bioc_if_needed("SingleCellExperiment")
install_bioc_if_needed("slingshot")
install_bioc_if_needed("DelayedMatrixStats")
install_cran_if_needed("FNN")
install_cran_if_needed("igraph")
install_cran_if_needed("ggridges")
install_cran_if_needed("patchwork")

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(slingshot)
  library(FNN)
  library(igraph)
  library(ggridges)
  library(patchwork)
})

log_msg("Pre-flight packages:")
log_msg("slingshot", as.character(packageVersion("slingshot")))
log_msg("SingleCellExperiment", as.character(packageVersion("SingleCellExperiment")))
if (requireNamespace("reticulate", quietly = TRUE)) {
  log_msg("reticulate", as.character(packageVersion("reticulate")))
} else {
  log_msg("reticulate unavailable")
}
log_msg("FNN", as.character(packageVersion("FNN")))
log_msg("igraph", as.character(packageVersion("igraph")))
log_msg("set.seed = 42")

obj_path <- file.path(
  "..",
  "05_恶性细胞分亚群与Neftel对照",
  "outputs",
  "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"
)
if (!file.exists(obj_path)) {
  obj_path <- file.path(
    "..",
    "05_鎭舵€х粏鑳炲垎浜氱兢涓嶯eftel瀵圭収",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"
  )
}
obj_path <- normalizePath(obj_path, winslash = "/", mustWork = TRUE)
coord_obj_path <- file.path(
  "..",
  "05_恶性细胞分亚群与Neftel对照",
  "outputs",
  "GBM.malignant.subtyped.umap_candidates.qs2"
)
if (!file.exists(coord_obj_path)) {
  coord_obj_path <- file.path(
    "..",
    "05_鎭舵€х粏鑳炲垎浜氱兢涓嶯eftel瀵圭収",
    "outputs",
    "GBM.malignant.subtyped.umap_candidates.qs2"
  )
}
coord_obj_path <- normalizePath(coord_obj_path, winslash = "/", mustWork = TRUE)
cytotrace_path <- file.path("outputs", "CytoTRACE2_per_cell.csv")
if (!file.exists(cytotrace_path)) stop("Missing cached CytoTRACE2 file: ", cytotrace_path)

log_msg("Input object:", obj_path)
log_msg("Coordinate object:", coord_obj_path)
log_msg("Coordinate reduction: umap_closer4")
log_msg("Cached CytoTRACE2:", cytotrace_path)
log_msg("Naming source:", file.path("scripts", "_naming.R"))
log_msg("Gene blacklist:", file.path("data", "gene_blacklist.txt"))

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
  left_join(cytotrace, by = c("cellID" = "cell_id"))
if (any(is.na(md$CytoTRACE2_score))) {
  log_msg("WARNING: CytoTRACE2 missing cells:", sum(is.na(md$CytoTRACE2_score)))
}

if (!"umap_closer4" %in% names(coord_obj@reductions)) stop("No malignant-only umap_closer4 reduction found in coordinate object.")
if (!"pca" %in% names(obj@reductions)) stop("No pre-existing pca reduction found.")
umap <- Seurat::Embeddings(coord_obj, reduction = "umap_closer4")
pca <- Seurat::Embeddings(obj, reduction = "pca")
overlap_cells <- intersect(rownames(umap), md$cellID)
log_msg("Coordinate cells:", nrow(umap))
log_msg("Metadata cells:", nrow(md))
log_msg("Coordinate-metadata overlap cells:", length(overlap_cells))
if (length(overlap_cells) != nrow(md)) {
  stop("umap_closer4 coordinates do not cover all metadata cells: overlap = ", length(overlap_cells), ", metadata = ", nrow(md))
}
if (!setequal(rownames(umap), md$cellID)) stop("umap_closer4 coordinate cell IDs do not match final labeled object cell IDs.")
if (!setequal(rownames(pca), md$cellID)) stop("PCA cell IDs do not match final labeled object cell IDs.")
log_msg("Using malignant-only UMAP dims:", paste(dim(umap), collapse = " x "))
log_msg("Using pre-existing PCA dims:", paste(dim(pca), collapse = " x "))

umap_df <- as.data.frame(umap[, 1:2, drop = FALSE])
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df <- umap_df |>
  rownames_to_column("cellID") |>
  left_join(md |> select(cellID, subtype_k4, subtype_label_source_metadata, subtype_label_original, subtype_label_final, subtype_short, Pt_number, CytoTRACE2_score), by = "cellID")
if (nrow(umap_df) != 28213) stop("Unexpected UMAP cell count: ", nrow(umap_df))
if (anyNA(umap_df$subtype_short)) stop("Missing subtype labels after coordinate/metadata join.")
if (anyNA(umap_df$UMAP_1) || anyNA(umap_df$UMAP_2)) stop("Missing umap_closer4 coordinates after join.")
umap_df$reduction <- "umap_closer4"

short_colors <- setNames(subtype_naming_mapping$color, subtype_naming_mapping$abbreviation)
short_levels <- subtype_naming_mapping$abbreviation
umap_df <- umap_df |>
  mutate(subtype_short = factor(subtype_short, levels = short_levels))

sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(dummy = Matrix::Matrix(0, nrow = 1, ncol = nrow(umap_df), sparse = TRUE)),
  colData = md |> select(cellID, subtype_k4, subtype_label_final, subtype_short)
)
colnames(sce) <- umap_df$cellID
SingleCellExperiment::reducedDim(sce, "UMAP") <- as.matrix(umap_df[, c("UMAP_1", "UMAP_2")])
cluster_labels <- factor(as.character(umap_df$subtype_short), levels = short_levels)

sl_start <- Sys.time()
sce <- slingshot::slingshot(
  sce,
  clusterLabels = cluster_labels,
  reducedDim = "UMAP",
  start.clus = "NPC-P",
  end.clus = NULL,
  approx_points = 150
)
log_msg("Slingshot seconds:", round(as.numeric(difftime(Sys.time(), sl_start, units = "secs")), 2))

pst <- as.data.frame(slingshot::slingPseudotime(sce))
weights <- as.data.frame(slingshot::slingCurveWeights(sce))
lineage_cols <- colnames(pst)
colnames(pst) <- paste0("lineage", seq_along(lineage_cols), "_pst")
colnames(weights) <- paste0("lineage", seq_along(lineage_cols), "_weight")

weight_mat <- as.matrix(weights)
lineage_assignment <- if (ncol(weight_mat) == 0) {
  rep(NA_character_, nrow(weight_mat))
} else {
  paste0("lineage", max.col(weight_mat, ties.method = "first"))
}
lineage_assignment[rowSums(weight_mat, na.rm = TRUE) == 0] <- NA_character_

pseudotime_per_cell <- bind_cols(
  umap_df |> select(cellID, subtype_label_final, subtype_k4, subtype_short, Pt_number),
  pst,
  weights,
  tibble(lineage_assignment = lineage_assignment)
)
readr::write_csv(pseudotime_per_cell, file.path("tables", "slingshot_pseudotime_per_cell.csv"))

lineages <- slingshot::slingLineages(sce)
log_msg("Slingshot lineage count:", length(lineages))
for (i in seq_along(lineages)) {
  log_msg("Lineage", i, "clusters:", paste(lineages[[i]], collapse = " -> "))
}

curve_coord <- lapply(seq_along(slingshot::slingCurves(sce)), function(i) {
  crv <- slingshot::slingCurves(sce)[[i]]
  tibble(
    lineage_id = paste0("lineage", i),
    point_id = seq_len(nrow(crv$s)),
    UMAP_1 = crv$s[, 1],
    UMAP_2 = crv$s[, 2],
    curve_order = seq_len(nrow(crv$s))
  )
}) |>
  bind_rows()
readr::write_csv(curve_coord, file.path("tables", "slingshot_curve_coordinates.csv"))

primary_pst_col <- "lineage1_pst"
lineage_summary <- lapply(seq_along(lineages), function(i) {
  pst_col <- paste0("lineage", i, "_pst")
  df <- pseudotime_per_cell |>
    filter(!is.na(.data[[pst_col]])) |>
    mutate(pseudotime_bin = ntile(.data[[pst_col]], 10))
  comp <- df |>
    dplyr::count(pseudotime_bin, subtype_label_final, subtype_short, name = "n_cells") |>
    group_by(pseudotime_bin) |>
    mutate(bin_fraction = n_cells / sum(n_cells)) |>
    ungroup() |>
    mutate(lineage_id = paste0("lineage", i))
  med <- df |>
    group_by(subtype_label_final, subtype_short) |>
    summarise(median_pseudotime = median(.data[[pst_col]], na.rm = TRUE), n_cells = n(), .groups = "drop") |>
    mutate(lineage_id = paste0("lineage", i), pseudotime_bin = NA_integer_, bin_fraction = NA_real_)
  bind_rows(
    comp |> mutate(median_pseudotime = NA_real_),
    med |> mutate(n_cells = n_cells)
  )
}) |>
  bind_rows()
lineage_meta <- tibble(
  lineage_id = paste0("lineage", seq_along(lineages)),
  start_cluster = vapply(lineages, function(x) x[1], character(1)),
  end_cluster = vapply(lineages, function(x) x[length(x)], character(1))
)
lineage_summary <- lineage_summary |>
  left_join(lineage_meta, by = "lineage_id")
readr::write_csv(lineage_summary, file.path("tables", "slingshot_lineage_summary.csv"))

scanpy_available <- FALSE
scanpy_version <- NA_character_
connectivity_method <- "kNN-based"
if (requireNamespace("reticulate", quietly = TRUE)) {
  scanpy_available <- tryCatch(reticulate::py_module_available("scanpy"), error = function(e) FALSE)
  if (scanpy_available) {
    scanpy_version <- tryCatch(as.character(reticulate::import("scanpy")$`__version__`), error = function(e) NA_character_)
  }
}
log_msg("scanpy available:", scanpy_available)
log_msg("scanpy version:", scanpy_version)

# Fallback kNN-based cluster connectivity on existing PCA.
connect_start <- Sys.time()
k <- 15
nn <- FNN::get.knn(pca[, seq_len(min(50, ncol(pca))), drop = FALSE], k = k)
cell_cluster <- setNames(as.character(umap_df$subtype_short), umap_df$cellID)
from <- rep(seq_len(nrow(nn$nn.index)), each = k)
to <- as.vector(t(nn$nn.index))
edges <- tibble(
  from_cluster = cell_cluster[rownames(pca)[from]],
  to_cluster = cell_cluster[rownames(pca)[to]]
) |>
  filter(from_cluster != to_cluster)
cluster_sizes <- table(cell_cluster)
cluster_pairs <- tidyr::expand_grid(
  from_cluster = short_levels,
  to_cluster = short_levels
)
edge_counts <- edges |>
  dplyr::count(from_cluster, to_cluster, name = "observed_edges")
connectivity <- cluster_pairs |>
  left_join(edge_counts, by = c("from_cluster", "to_cluster")) |>
  mutate(observed_edges = tidyr::replace_na(observed_edges, 0L)) |>
  rowwise() |>
  mutate(
    possible_edges = as.numeric(cluster_sizes[from_cluster]) * k,
    edge_fraction = observed_edges / possible_edges,
    expected_fraction = as.numeric(cluster_sizes[to_cluster]) / (sum(cluster_sizes) - 1),
    observed_expected_ratio = ifelse(expected_fraction > 0, edge_fraction / expected_fraction, NA_real_)
  ) |>
  ungroup()
connectivity_sym <- connectivity |>
  dplyr::select(from_cluster, to_cluster, edge_fraction, expected_fraction, observed_expected_ratio) |>
  left_join(
    connectivity |> dplyr::select(
      from_cluster = to_cluster,
      to_cluster = from_cluster,
      edge_fraction_rev = edge_fraction,
      expected_fraction_rev = expected_fraction,
      observed_expected_ratio_rev = observed_expected_ratio
    ),
    by = c("from_cluster", "to_cluster")
  ) |>
  mutate(
    observed_expected_ratio_sym = ifelse(
      from_cluster == to_cluster,
      NA_real_,
      rowMeans(cbind(observed_expected_ratio, observed_expected_ratio_rev), na.rm = TRUE)
    )
  )
max_ratio <- max(connectivity_sym$observed_expected_ratio_sym, na.rm = TRUE)
connectivity_sym <- connectivity_sym |>
  mutate(connectivity = ifelse(from_cluster == to_cluster, 1, observed_expected_ratio_sym / max_ratio)) |>
  left_join(subtype_naming_mapping |> dplyr::select(abbreviation, subtype_label_final), by = c("from_cluster" = "abbreviation")) |>
  dplyr::rename(from_label = subtype_label_final) |>
  left_join(subtype_naming_mapping |> dplyr::select(abbreviation, subtype_label_final), by = c("to_cluster" = "abbreviation")) |>
  dplyr::rename(to_label = subtype_label_final) |>
  mutate(method = connectivity_method)
readr::write_csv(connectivity_sym, file.path("tables", "paga_connectivity_matrix.csv"))
log_msg("Connectivity method:", connectivity_method)
log_msg("Connectivity seconds:", round(as.numeric(difftime(Sys.time(), connect_start, units = "secs")), 2))
log_msg("Connectivity edges threshold 0.1:", sum(connectivity_sym$from_cluster < connectivity_sym$to_cluster & connectivity_sym$connectivity >= 0.1))

source_a <- umap_df |> select(cellID, UMAP_1, UMAP_2, subtype_k4, subtype_label_final, subtype_short, reduction)
source_b <- bind_rows(
  source_a |> mutate(layer = "cell_background", lineage_id = NA_character_, point_id = NA_integer_),
  curve_coord |> mutate(cellID = NA_character_, subtype_k4 = NA, subtype_label_final = NA, subtype_short = NA, reduction = "umap_closer4", layer = "slingshot_curve")
) |>
  select(layer, cellID, lineage_id, point_id, UMAP_1, UMAP_2, subtype_k4, subtype_label_final, subtype_short, reduction)
source_c <- umap_df |>
  select(cellID, UMAP_1, UMAP_2, subtype_k4, subtype_label_final, subtype_short, reduction) |>
  left_join(pseudotime_per_cell |> select(cellID, all_of(primary_pst_col)), by = "cellID")
source_d <- connectivity_sym
source_e <- pseudotime_per_cell |>
  select(cellID, subtype_k4, subtype_label_final, subtype_short, all_of(primary_pst_col), lineage_assignment)

readr::write_csv(source_a, file.path("figures/source_data", "02_panel_A.csv"))
readr::write_csv(source_b, file.path("figures/source_data", "02_panel_B.csv"))
readr::write_csv(source_c, file.path("figures/source_data", "02_panel_C.csv"))
readr::write_csv(source_d, file.path("figures/source_data", "02_panel_D.csv"))
readr::write_csv(source_e, file.path("figures/source_data", "02_panel_E.csv"))

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
root_center <- source_a |>
  filter(subtype_short == "NPC-P") |>
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
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
  geom_point(data = root_center, aes(UMAP_1, UMAP_2), inherit.aes = FALSE, shape = 21, fill = "white", color = "black", size = 1.8, stroke = 0.45) +
  geom_text(data = curve_labels, aes(label = gsub("lineage", "L", lineage_id)), size = 2.2, vjust = -0.4, color = "black") +
  scale_color_manual(values = short_colors, guide = "none") +
  coord_equal() +
  labs(title = "B. Slingshot curves", subtitle = "Root: NPC-P", x = "UMAP 1", y = "UMAP 2") +
  base_theme

p_c <- ggplot(source_c, aes(UMAP_1, UMAP_2, color = .data[[primary_pst_col]])) +
  geom_point(size = 0.12, alpha = 0.72) +
  scale_color_viridis_c(option = "viridis", name = "Pseudotime", na.value = "grey88") +
  coord_equal() +
  labs(title = "C. Primary pseudotime", subtitle = "Lineage 1", x = "UMAP 1", y = "UMAP 2") +
  base_theme

node_df <- source_a |>
  dplyr::count(subtype_short, subtype_label_final, name = "n_cells") |>
  mutate(
    angle = seq(0, 2 * pi, length.out = n() + 1)[- (n() + 1)],
    x = cos(angle),
    y = sin(angle),
    node_size = log1p(n_cells)
  )
edge_df <- connectivity_sym |>
  filter(from_cluster < to_cluster, connectivity >= 0.1) |>
  left_join(node_df |> select(from_cluster = subtype_short, x_from = x, y_from = y), by = "from_cluster") |>
  left_join(node_df |> select(to_cluster = subtype_short, x_to = x, y_to = y), by = "to_cluster")
p_d <- ggplot() +
  geom_segment(data = edge_df, aes(x = x_from, y = y_from, xend = x_to, yend = y_to, linewidth = connectivity), color = "grey35", alpha = 0.75) +
  geom_point(data = node_df, aes(x, y, fill = subtype_short, size = node_size), shape = 21, color = "black", stroke = 0.35) +
  geom_text(data = node_df, aes(x, y, label = subtype_short), size = 2.4, vjust = -1.0) +
  scale_fill_manual(values = short_colors, guide = "none") +
  scale_size_continuous(range = c(4, 8), guide = "none") +
  scale_linewidth_continuous(range = c(0.3, 2.0), name = "Connectivity") +
  coord_equal(xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35)) +
  labs(title = "D. kNN-based connectivity", subtitle = "Threshold = 0.1", x = NULL, y = NULL) +
  base_theme

p_e <- source_e |>
  filter(!is.na(.data[[primary_pst_col]])) |>
  mutate(subtype_short = factor(subtype_short, levels = rev(short_levels))) |>
  ggplot(aes(x = .data[[primary_pst_col]], y = subtype_short, fill = subtype_short)) +
  ggridges::geom_density_ridges(scale = 1.15, alpha = 0.72, linewidth = 0.25, color = "white") +
  scale_fill_manual(values = short_colors, guide = "none") +
  labs(title = "E. Pseudotime density", x = "Lineage 1 pseudotime", y = NULL) +
  theme_bw(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 6)
  )

caption <- "Pseudotime represents an inferred transcriptional ordering and does not imply real-time progression or deterministic lineage commitment."
p_all <- (p_a | p_b) / (p_c | p_d) / p_e +
  patchwork::plot_annotation(caption = caption, theme = theme(plot.caption = element_text(size = 6, hjust = 0)))

pdf(file.path("figures", "02_trajectory_primary.pdf"), width = 8.4, height = 10.2, useDingbats = FALSE)
print(p_all)
dev.off()
log_msg("Wrote figure:", file.path("figures", "02_trajectory_primary.pdf"))
ggsave(file.path("figures", "02_trajectory_primary_preview.png"), p_all, width = 8.4, height = 10.2, dpi = 180, bg = "white")
log_msg("Wrote preview:", file.path("figures", "02_trajectory_primary_preview.png"))

session_info <- capture.output(sessionInfo())
cat(session_info, sep = "\n", file = file.path("logs", "02_session_info.txt"))
log_msg("Session info:", file.path("logs", "02_session_info.txt"))
log_msg("Total elapsed seconds:", round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 2))
log_msg("STOP: Phase 1 complete. Await audit before Phase 2.")
