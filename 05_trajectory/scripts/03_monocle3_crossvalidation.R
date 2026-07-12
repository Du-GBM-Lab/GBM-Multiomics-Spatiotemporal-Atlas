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
  library(patchwork)
  library(ggalluvial)
  library(SingleCellExperiment)
  library(monocle3)
  library(igraph)
})

# Phase 2A: Monocle3 cross-validation.
# Read-only for source objects: no UMAP / Seurat clustering rerun; no qs2 writeback.

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

log_file <- file.path("logs", "03_monocle3_crossvalidation.log")
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_msg("Phase 2A Monocle3 cross-validation started")
log_msg("set.seed = 42")
log_msg("monocle3", as.character(packageVersion("monocle3")))
log_msg("SingleCellExperiment", as.character(packageVersion("SingleCellExperiment")))
log_msg("leidenbase", ifelse(requireNamespace("leidenbase", quietly = TRUE), as.character(packageVersion("leidenbase")), "unavailable"))
log_msg("sf", ifelse(requireNamespace("sf", quietly = TRUE), as.character(packageVersion("sf")), "unavailable"))
log_msg("terra", ifelse(requireNamespace("terra", quietly = TRUE), as.character(packageVersion("terra")), "unavailable"))

obj_path <- normalizePath(file.path("..", "05_恶性细胞分亚群与Neftel对照", "outputs", "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"), winslash = "/", mustWork = TRUE)
coord_obj_path <- normalizePath(file.path("..", "05_恶性细胞分亚群与Neftel对照", "outputs", "GBM.malignant.subtyped.umap_candidates.qs2"), winslash = "/", mustWork = TRUE)
cytotrace_path <- file.path("outputs", "CytoTRACE2_per_cell.csv")
if (!file.exists(cytotrace_path)) stop("Missing cached CytoTRACE2 file: ", cytotrace_path)
slingshot_path <- file.path("tables", "slingshot_pseudotime_per_cell.csv")
if (!file.exists(slingshot_path)) stop("Missing Slingshot pseudotime file: ", slingshot_path)

log_msg("Input object:", obj_path)
log_msg("Coordinate object:", coord_obj_path)
log_msg("Coordinate reduction: umap_closer4")
log_msg("Cached CytoTRACE2:", cytotrace_path)
log_msg("Slingshot table:", slingshot_path)

get_counts <- function(obj) {
  assay <- "RNA"
  if (!assay %in% names(obj@assays)) assay <- DefaultAssay(obj)
  counts <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(counts)) {
    counts <- Seurat::GetAssayData(obj, assay = assay, slot = "counts")
  }
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  counts
}

load_start <- Sys.time()
obj <- qs2::qs_read(obj_path)
coord_obj <- qs2::qs_read(coord_obj_path)
log_msg("Object load seconds:", round(as.numeric(difftime(Sys.time(), load_start, units = "secs")), 2))

md <- obj@meta.data |>
  rownames_to_column("cellID") |>
  recode_subtype()
cytotrace <- readr::read_csv(cytotrace_path, show_col_types = FALSE)
if (!all(c("cell_id", "CytoTRACE2_score") %in% colnames(cytotrace))) {
  stop("Cached CytoTRACE2 file lacks cell_id / CytoTRACE2_score.")
}
md <- md |>
  left_join(cytotrace, by = c("cellID" = "cell_id"))
if (anyNA(md$CytoTRACE2_score)) log_msg("WARNING CytoTRACE2 missing cells:", sum(is.na(md$CytoTRACE2_score)))

if (!"umap_closer4" %in% names(coord_obj@reductions)) stop("Coordinate object lacks umap_closer4.")
umap <- Seurat::Embeddings(coord_obj, "umap_closer4")
if (!setequal(rownames(umap), md$cellID)) stop("UMAP coordinate cell IDs do not match metadata cell IDs.")
umap <- umap[md$cellID, 1:2, drop = FALSE]
colnames(umap) <- c("UMAP_1", "UMAP_2")
log_msg("UMAP cells:", nrow(umap))
log_msg("Metadata cells:", nrow(md))

counts <- get_counts(obj)
counts <- counts[, md$cellID, drop = FALSE]
log_msg("Counts dims:", paste(dim(counts), collapse = " x "))

cell_metadata <- md |>
  dplyr::select(cellID, subtype_k4, subtype_label_final, subtype_short, Pt_number, CytoTRACE2_score) |>
  as.data.frame()
rownames(cell_metadata) <- cell_metadata$cellID
gene_metadata <- data.frame(
  gene_short_name = rownames(counts),
  row.names = rownames(counts),
  stringsAsFactors = FALSE
)

cds <- monocle3::new_cell_data_set(
  expression_data = counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata,
  verbose = FALSE
)
SingleCellExperiment::reducedDim(cds, "UMAP") <- as.matrix(umap)
cluster_factor <- factor(cell_metadata$subtype_label_final, levels = subtype_naming_mapping$subtype_label_final)
names(cluster_factor) <- rownames(cell_metadata)
partition_vector <- setNames(rep(1, ncol(cds)), colnames(cds))
membership <- setNames(as.integer(factor(cell_metadata$subtype_label_final, levels = subtype_naming_mapping$subtype_label_final)), colnames(cds))
cds@clusters[["UMAP"]] <- list(
  cluster_result = list(
    g = NULL,
    relations = NULL,
    distMatrix = NULL,
    coord = NULL,
    edge_links = NULL,
    optim_res = list(membership = membership, modularity = NA_real_)
  ),
  partitions = partition_vector,
  clusters = cluster_factor
)
log_msg("Injected UMAP reduction and subtype_label_final cluster identity; cluster_cells not run.")

learn_start <- Sys.time()
cds <- monocle3::learn_graph(cds, use_partition = FALSE, close_loop = FALSE, verbose = FALSE)
log_msg("learn_graph seconds:", round(as.numeric(difftime(Sys.time(), learn_start, units = "secs")), 2))

root_candidates <- cell_metadata |>
  dplyr::filter(subtype_short == "NPC-P", !is.na(CytoTRACE2_score)) |>
  dplyr::arrange(dplyr::desc(CytoTRACE2_score))
root_pool <- root_candidates |>
  dplyr::slice_head(n = min(100, nrow(root_candidates))) |>
  dplyr::pull(cellID)
if (length(root_pool) == 0) stop("No NPC-P high-CytoTRACE2 root cells available.")
log_msg("Root cells: top NPC-P CytoTRACE2 cells n =", length(root_pool))
order_start <- Sys.time()
cds <- monocle3::order_cells(cds, reduction_method = "UMAP", root_cells = root_pool, verbose = FALSE)
log_msg("order_cells seconds:", round(as.numeric(difftime(Sys.time(), order_start, units = "secs")), 2))

mono_pst <- monocle3::pseudotime(cds)
mono_partition <- tryCatch(as.character(monocle3::partitions(cds, reduction_method = "UMAP")), error = function(e) rep(NA_character_, ncol(cds)))
mono_clusters <- tryCatch(as.character(monocle3::clusters(cds, reduction_method = "UMAP")), error = function(e) rep(NA_character_, ncol(cds)))
names(mono_partition) <- colnames(cds)
names(mono_clusters) <- colnames(cds)

slingshot <- readr::read_csv(slingshot_path, show_col_types = FALSE)
mono_cell <- cell_metadata |>
  rownames_to_column("cellID2") |>
  dplyr::select(-cellID2) |>
  dplyr::mutate(
    monocle3_pst = as.numeric(mono_pst[cellID]),
    monocle3_partition = mono_partition[cellID],
    monocle3_cluster = mono_clusters[cellID]
  ) |>
  left_join(
    slingshot |>
      dplyr::select(cellID, slingshot_assigned_lineage = lineage_assignment, dplyr::starts_with("lineage")),
    by = "cellID"
  )

g <- monocle3::principal_graph(cds)[["UMAP"]]
dp_mst <- monocle3::principal_graph_aux(cds)[["UMAP"]]$dp_mst
node_df <- as.data.frame(t(dp_mst))
colnames(node_df)[1:2] <- c("UMAP_1", "UMAP_2")
node_df <- node_df |>
  rownames_to_column("node_id")
if (!"node_id" %in% colnames(node_df) || any(is.na(node_df$node_id)) || any(node_df$node_id == "")) {
  node_df$node_id <- igraph::V(g)$name
}
node_df <- node_df |>
  dplyr::mutate(node_degree = igraph::degree(g)[node_id])
edge_df <- igraph::as_data_frame(g, what = "edges") |>
  as_tibble() |>
  dplyr::rename(from_node = "from", to_node = "to") |>
  dplyr::left_join(node_df |> dplyr::select(from_node = node_id, from_UMAP_1 = UMAP_1, from_UMAP_2 = UMAP_2), by = "from_node") |>
  dplyr::left_join(node_df |> dplyr::select(to_node = node_id, to_UMAP_1 = UMAP_1, to_UMAP_2 = UMAP_2), by = "to_node")

closest_vertex <- monocle3::principal_graph_aux(cds)[["UMAP"]]$pr_graph_cell_proj_closest_vertex
closest_idx <- as.integer(closest_vertex[colnames(cds), 1])
closest_node <- colnames(dp_mst)[closest_idx]
names(closest_node) <- colnames(cds)

leaf_nodes <- node_df |> dplyr::filter(node_degree == 1) |> dplyr::pull(node_id)
root_nodes <- unique(closest_node[root_pool])
root_node <- names(sort(table(root_nodes), decreasing = TRUE))[1]
leaf_nodes_for_terminal <- setdiff(leaf_nodes, root_node)
if (length(leaf_nodes_for_terminal) == 0) leaf_nodes_for_terminal <- leaf_nodes
dist_to_leaves <- igraph::distances(g, v = closest_node, to = leaf_nodes_for_terminal, weights = NA)
terminal_assignment <- leaf_nodes_for_terminal[max.col(-as.matrix(dist_to_leaves), ties.method = "first")]
names(terminal_assignment) <- names(closest_node)

mono_cell <- mono_cell |>
  dplyr::mutate(
    closest_graph_node = closest_node[cellID],
    monocle3_terminal = terminal_assignment[cellID],
    monocle3_terminal = ifelse(is.na(monocle3_terminal), "unassigned", monocle3_terminal)
  )

terminal_detection <- mono_cell |>
  dplyr::filter(monocle3_terminal != "unassigned") |>
  dplyr::count(monocle3_terminal, subtype_label_final, subtype_short, name = "n_cells") |>
  dplyr::group_by(monocle3_terminal) |>
  dplyr::mutate(
    terminal_fraction = n_cells / sum(n_cells),
    terminal_total_cells = sum(n_cells)
  ) |>
  dplyr::ungroup() |>
  dplyr::left_join(node_df |> dplyr::select(monocle3_terminal = node_id, terminal_UMAP_1 = UMAP_1, terminal_UMAP_2 = UMAP_2, terminal_degree = node_degree), by = "monocle3_terminal")
terminal_summary <- terminal_detection |>
  dplyr::group_by(monocle3_terminal, terminal_total_cells, terminal_UMAP_1, terminal_UMAP_2) |>
  dplyr::slice_max(terminal_fraction, n = 1, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::transmute(monocle3_terminal, terminal_total_cells, terminal_UMAP_1, terminal_UMAP_2, dominant_subtype = subtype_short, dominant_fraction = terminal_fraction)
log_msg("Terminal nodes detected:", nrow(terminal_summary))
log_msg("Terminal dominant subtypes:", paste(paste0(terminal_summary$monocle3_terminal, "=", terminal_summary$dominant_subtype), collapse = "; "))

method_corr <- lapply(paste0("lineage", 1:3), function(lin) {
  pst_col <- paste0(lin, "_pst")
  df <- mono_cell |>
    dplyr::filter(slingshot_assigned_lineage == lin, !is.na(.data[[pst_col]]), is.finite(monocle3_pst))
  if (nrow(df) < 3) {
    return(tibble(lineage_id = lin, n_cells = nrow(df), spearman_rho = NA_real_, p_value = NA_real_))
  }
  ct <- suppressWarnings(cor.test(df[[pst_col]], df$monocle3_pst, method = "spearman", exact = FALSE))
  tibble(lineage_id = lin, n_cells = nrow(df), spearman_rho = unname(ct$estimate), p_value = ct$p.value)
}) |>
  bind_rows()
readr::write_csv(method_corr, file.path("tables", "pseudotime_method_correlation.csv"))

graph_table <- bind_rows(
  node_df |>
    dplyr::transmute(record_type = "node", node_id, from_node = NA_character_, to_node = NA_character_, UMAP_1, UMAP_2, from_UMAP_1 = NA_real_, from_UMAP_2 = NA_real_, to_UMAP_1 = NA_real_, to_UMAP_2 = NA_real_, node_degree),
  edge_df |>
    dplyr::transmute(record_type = "edge", node_id = NA_character_, from_node, to_node, UMAP_1 = NA_real_, UMAP_2 = NA_real_, from_UMAP_1, from_UMAP_2, to_UMAP_1, to_UMAP_2, node_degree = NA_integer_)
)

readr::write_csv(
  mono_cell |>
    dplyr::select(cellID, subtype_label_final, subtype_short, subtype_k4, Pt_number, monocle3_pst, monocle3_partition, monocle3_cluster, closest_graph_node, monocle3_terminal, slingshot_assigned_lineage),
  file.path("tables", "monocle3_pseudotime_per_cell.csv")
)
readr::write_csv(graph_table, file.path("tables", "monocle3_principal_graph_edges.csv"))
readr::write_csv(terminal_detection, file.path("tables", "monocle3_terminal_detection.csv"))

umap_df <- as.data.frame(umap) |>
  rownames_to_column("cellID") |>
  dplyr::left_join(mono_cell |> dplyr::select(cellID, subtype_k4, subtype_label_final, subtype_short, monocle3_pst, monocle3_terminal, slingshot_assigned_lineage), by = "cellID") |>
  dplyr::mutate(subtype_short = factor(subtype_short, levels = subtype_naming_mapping$abbreviation))

plot_edges <- edge_df
plot_nodes <- node_df |>
  dplyr::mutate(
    node_role = case_when(
      node_id == root_node ~ "root",
      node_id %in% leaf_nodes_for_terminal ~ "terminal",
      TRUE ~ "internal"
    )
  )
alluvial_df <- mono_cell |>
  dplyr::count(slingshot_assigned_lineage, monocle3_terminal, name = "n_cells") |>
  dplyr::filter(!is.na(slingshot_assigned_lineage), !is.na(monocle3_terminal))
scatter_df <- lapply(paste0("lineage", 1:3), function(lin) {
  pst_col <- paste0(lin, "_pst")
  mono_cell |>
    dplyr::filter(slingshot_assigned_lineage == lin, !is.na(.data[[pst_col]]), is.finite(monocle3_pst)) |>
    dplyr::transmute(
      cellID,
      lineage_id = lin,
      subtype_short,
      slingshot_pst = .data[[pst_col]],
      monocle3_pst
    )
}) |>
  bind_rows() |>
  dplyr::mutate(lineage_id = factor(lineage_id, levels = paste0("lineage", 1:3), labels = paste0("L", 1:3)))
scatter_labels <- method_corr |>
  dplyr::mutate(
    lineage_id = factor(lineage_id, levels = paste0("lineage", 1:3), labels = paste0("L", 1:3)),
    label = paste0("rho = ", sprintf("%.2f", spearman_rho), "\np = ", format.pval(p_value, digits = 2, eps = 1e-300), "\nn = ", n_cells)
  )

readr::write_csv(umap_df, file.path("figures/source_data", "03_panel_A_B_cells.csv"))
readr::write_csv(plot_edges, file.path("figures/source_data", "03_panel_A_graph_edges.csv"))
readr::write_csv(plot_nodes, file.path("figures/source_data", "03_panel_A_graph_nodes.csv"))
readr::write_csv(scatter_df, file.path("figures/source_data", "03_panel_C.csv"))
readr::write_csv(alluvial_df, file.path("figures/source_data", "03_panel_D.csv"))

short_colors <- setNames(subtype_naming_mapping$color, subtype_naming_mapping$abbreviation)
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

p_a <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = subtype_short), size = 0.08, alpha = 0.35) +
  geom_segment(data = plot_edges, aes(x = from_UMAP_1, y = from_UMAP_2, xend = to_UMAP_1, yend = to_UMAP_2), inherit.aes = FALSE, color = "black", linewidth = 0.35) +
  geom_point(data = plot_nodes |> dplyr::filter(node_role == "internal"), aes(UMAP_1, UMAP_2), inherit.aes = FALSE, color = "black", size = 0.45) +
  geom_point(data = plot_nodes |> dplyr::filter(node_role == "root"), aes(UMAP_1, UMAP_2), inherit.aes = FALSE, shape = 21, fill = "red", color = "black", size = 2.0, stroke = 0.35) +
  geom_point(data = plot_nodes |> dplyr::filter(node_role == "terminal"), aes(UMAP_1, UMAP_2), inherit.aes = FALSE, shape = 17, color = "black", size = 1.8) +
  scale_color_manual(values = short_colors, name = NULL) +
  coord_equal() +
  labs(title = "A. Monocle3 principal graph", subtitle = "Root: high-CytoTRACE2 NPC-P cells", x = "UMAP 1", y = "UMAP 2") +
  base_theme

p_b <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = monocle3_pst)) +
  geom_point(size = 0.10, alpha = 0.72) +
  scale_color_viridis_c(option = "viridis", name = "Pseudotime", na.value = "grey88") +
  coord_equal() +
  labs(title = "B. Monocle3 pseudotime", x = "UMAP 1", y = "UMAP 2") +
  base_theme

p_c <- ggplot(scatter_df, aes(slingshot_pst, monocle3_pst, color = subtype_short)) +
  geom_point(size = 0.12, alpha = 0.35) +
  facet_wrap(~ lineage_id, scales = "free_x") +
  geom_text(data = scatter_labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.05, vjust = 1.05, size = 2.0) +
  scale_color_manual(values = short_colors, name = NULL) +
  labs(title = "C. Slingshot vs Monocle3 pseudotime", x = "Slingshot pseudotime", y = "Monocle3 pseudotime") +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), legend.position = "bottom", strip.text = element_text(size = 7), plot.title = element_text(size = 8))

p_d <- ggplot(alluvial_df, aes(axis1 = slingshot_assigned_lineage, axis2 = monocle3_terminal, y = n_cells)) +
  ggalluvial::geom_alluvium(aes(fill = slingshot_assigned_lineage), alpha = 0.72, width = 0.18) +
  ggalluvial::geom_stratum(width = 0.18, fill = "grey95", color = "grey30", linewidth = 0.25) +
  ggplot2::geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.2) +
  scale_x_discrete(limits = c("Slingshot lineage", "Monocle3 terminal"), expand = c(0.08, 0.08)) +
  labs(title = "D. Lineage-to-terminal correspondence", y = "Cells", x = NULL) +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), legend.position = "none", plot.title = element_text(size = 8), axis.text.y = element_text(size = 6))

p_all <- (p_a | p_b) / (p_c | p_d) +
  patchwork::plot_annotation(
    caption = "Monocle3 cross-validation uses the malignant-only umap_closer4 coordinates and does not rerun UMAP or Seurat clustering.",
    theme = theme(plot.caption = element_text(size = 6, hjust = 0))
  )
pdf(file.path("figures", "03_monocle3_vs_slingshot.pdf"), width = 10.2, height = 8.2, useDingbats = FALSE)
print(p_all)
dev.off()
ggsave(file.path("figures", "03_monocle3_vs_slingshot_preview.png"), p_all, width = 10.2, height = 8.2, dpi = 180, bg = "white")
log_msg("Wrote figure:", file.path("figures", "03_monocle3_vs_slingshot.pdf"))
log_msg("Wrote preview:", file.path("figures", "03_monocle3_vs_slingshot_preview.png"))

session_info <- capture.output(sessionInfo())
cat(session_info, sep = "\n", file = file.path("logs", "03_session_info.txt"))
log_msg("Session info:", file.path("logs", "03_session_info.txt"))
log_msg("Total elapsed seconds:", round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 2))
log_msg("STOP: Phase 2A complete. Await audit before Phase 2B.")
