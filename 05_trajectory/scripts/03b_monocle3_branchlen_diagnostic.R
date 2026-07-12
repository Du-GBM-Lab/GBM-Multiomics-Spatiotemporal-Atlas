suppressPackageStartupMessages({
  .libPaths(c("<R_LIBS>", "<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(qs2)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(SingleCellExperiment)
  library(monocle3)
  library(igraph)
})

# Phase 2A supplemental diagnostic:
# Monocle3 minimal_branch_len sensitivity.
# No Slingshot rerun, no UMAP rerun, no Seurat clustering rerun, no figure.

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
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

source(file.path("scripts", "_naming.R"))

log_file <- file.path("logs", "03b_diagnostic.log")
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_msg("Phase 2A supplemental minimal_branch_len diagnostic started")
log_msg("set.seed = 42")
log_msg("monocle3", as.character(packageVersion("monocle3")))
log_msg("sf", ifelse(requireNamespace("sf", quietly = TRUE), as.character(packageVersion("sf")), "unavailable"))

obj_path <- normalizePath(file.path("..", "05_恶性细胞分亚群与Neftel对照", "outputs", "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"), winslash = "/", mustWork = TRUE)
coord_obj_path <- normalizePath(file.path("..", "05_恶性细胞分亚群与Neftel对照", "outputs", "GBM.malignant.subtyped.umap_candidates.qs2"), winslash = "/", mustWork = TRUE)
cytotrace_path <- file.path("outputs", "CytoTRACE2_per_cell.csv")
slingshot_path <- file.path("tables", "slingshot_pseudotime_per_cell.csv")
if (!file.exists(cytotrace_path)) stop("Missing cached CytoTRACE2 file: ", cytotrace_path)
if (!file.exists(slingshot_path)) stop("Missing Slingshot pseudotime file: ", slingshot_path)

get_counts <- function(obj) {
  assay <- "RNA"
  if (!assay %in% names(obj@assays)) assay <- DefaultAssay(obj)
  counts <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) NULL
  )
  if (is.null(counts)) counts <- Seurat::GetAssayData(obj, assay = assay, slot = "counts")
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
md <- md |>
  dplyr::left_join(cytotrace, by = c("cellID" = "cell_id"))

umap <- Seurat::Embeddings(coord_obj, "umap_closer4")
if (!setequal(rownames(umap), md$cellID)) stop("UMAP cell IDs do not match metadata.")
umap <- umap[md$cellID, 1:2, drop = FALSE]
colnames(umap) <- c("UMAP_1", "UMAP_2")

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

base_cds <- monocle3::new_cell_data_set(
  expression_data = counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata,
  verbose = FALSE
)
SingleCellExperiment::reducedDim(base_cds, "UMAP") <- as.matrix(umap)
cluster_factor <- factor(cell_metadata$subtype_label_final, levels = subtype_naming_mapping$subtype_label_final)
names(cluster_factor) <- rownames(cell_metadata)
base_cds@clusters[["UMAP"]] <- list(
  cluster_result = list(
    g = NULL,
    relations = NULL,
    distMatrix = NULL,
    coord = NULL,
    edge_links = NULL,
    optim_res = list(
      membership = setNames(as.integer(factor(cell_metadata$subtype_label_final, levels = subtype_naming_mapping$subtype_label_final)), colnames(base_cds)),
      modularity = NA_real_
    )
  ),
  partitions = setNames(rep(1, ncol(base_cds)), colnames(base_cds)),
  clusters = cluster_factor
)

root_candidates <- cell_metadata |>
  dplyr::filter(subtype_short == "NPC-P", !is.na(CytoTRACE2_score)) |>
  dplyr::arrange(dplyr::desc(CytoTRACE2_score))
root_pool <- root_candidates |>
  dplyr::slice_head(n = min(100, nrow(root_candidates))) |>
  dplyr::pull(cellID)
if (length(root_pool) == 0) stop("No NPC-P root cells available.")
log_msg("Root cells: top NPC-P CytoTRACE2 cells n =", length(root_pool))

slingshot <- readr::read_csv(slingshot_path, show_col_types = FALSE)

run_one <- function(run_id, minimal_branch_len) {
  log_msg("Run started:", run_id, "minimal_branch_len =", minimal_branch_len)
  set.seed(42)
  cds <- base_cds
  run_start <- Sys.time()
  cds <- monocle3::learn_graph(
    cds,
    use_partition = FALSE,
    close_loop = FALSE,
    learn_graph_control = list(minimal_branch_len = minimal_branch_len),
    verbose = FALSE
  )
  cds <- monocle3::order_cells(cds, reduction_method = "UMAP", root_cells = root_pool, verbose = FALSE)
  log_msg("Run seconds:", run_id, round(as.numeric(difftime(Sys.time(), run_start, units = "secs")), 2))

  mono_pst <- monocle3::pseudotime(cds)
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

  closest_vertex <- monocle3::principal_graph_aux(cds)[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_idx <- as.integer(closest_vertex[colnames(cds), 1])
  closest_node <- colnames(dp_mst)[closest_idx]
  names(closest_node) <- colnames(cds)
  leaf_nodes <- node_df |>
    dplyr::filter(node_degree == 1) |>
    dplyr::pull(node_id)
  root_nodes <- unique(closest_node[root_pool])
  root_node <- names(sort(table(root_nodes), decreasing = TRUE))[1]
  terminal_nodes <- setdiff(leaf_nodes, root_node)
  if (length(terminal_nodes) == 0) terminal_nodes <- leaf_nodes
  dist_to_terminals <- igraph::distances(g, v = closest_node, to = terminal_nodes, weights = NA)
  terminal_assignment <- terminal_nodes[max.col(-as.matrix(dist_to_terminals), ties.method = "first")]
  names(terminal_assignment) <- names(closest_node)

  mono_cell <- cell_metadata |>
    dplyr::mutate(
      monocle3_pst = as.numeric(mono_pst[cellID]),
      closest_graph_node = closest_node[cellID],
      monocle3_terminal = terminal_assignment[cellID],
      monocle3_terminal = ifelse(is.na(monocle3_terminal), "unassigned", monocle3_terminal)
    ) |>
    dplyr::left_join(
      slingshot |>
        dplyr::select(cellID, slingshot_assigned_lineage = lineage_assignment, dplyr::starts_with("lineage")),
      by = "cellID"
    )

  terminal_detection <- mono_cell |>
    dplyr::filter(monocle3_terminal != "unassigned") |>
    dplyr::count(monocle3_terminal, subtype_label_final, subtype_short, name = "n_cells") |>
    dplyr::group_by(monocle3_terminal) |>
    dplyr::mutate(
      terminal_fraction = n_cells / sum(n_cells),
      terminal_total_cells = sum(n_cells)
    ) |>
    dplyr::ungroup()

  terminal_summary <- terminal_detection |>
    dplyr::group_by(monocle3_terminal, terminal_total_cells) |>
    dplyr::slice_max(terminal_fraction, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      run_id,
      minimal_branch_len,
      monocle3_terminal,
      terminal_total_cells,
      dominant_subtype = as.character(subtype_short),
      dominant_fraction = terminal_fraction,
      biological_terminal = dominant_subtype %in% c("OPC-M", "MES-V", "MES-I") & dominant_fraction >= 0.7
    )

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

  diagnostic <- tibble(
    run_id = run_id,
    minimal_branch_len = minimal_branch_len,
    n_terminals = nrow(terminal_summary),
    n_biological_terminals = sum(terminal_summary$biological_terminal),
    n_spurious_or_mixed_terminals = n_terminals - n_biological_terminals,
    L1_rho = method_corr$spearman_rho[method_corr$lineage_id == "lineage1"],
    L2_rho = method_corr$spearman_rho[method_corr$lineage_id == "lineage2"],
    L3_rho = method_corr$spearman_rho[method_corr$lineage_id == "lineage3"],
    L1_p = method_corr$p_value[method_corr$lineage_id == "lineage1"],
    L2_p = method_corr$p_value[method_corr$lineage_id == "lineage2"],
    L3_p = method_corr$p_value[method_corr$lineage_id == "lineage3"],
    L1_n = method_corr$n_cells[method_corr$lineage_id == "lineage1"],
    L2_n = method_corr$n_cells[method_corr$lineage_id == "lineage2"],
    L3_n = method_corr$n_cells[method_corr$lineage_id == "lineage3"]
  )

  list(
    diagnostic = diagnostic,
    terminal_summary = terminal_summary,
    terminal_detection = terminal_detection |> dplyr::mutate(run_id = run_id, minimal_branch_len = minimal_branch_len)
  )
}

runs <- list(
  run_default = 10,
  run_bl15 = 15,
  run_bl20 = 20
)
results <- lapply(names(runs), function(run_id) run_one(run_id, runs[[run_id]]))

diagnostic <- bind_rows(lapply(results, `[[`, "diagnostic"))
terminal_summary <- bind_rows(lapply(results, `[[`, "terminal_summary"))
terminal_detection <- bind_rows(lapply(results, `[[`, "terminal_detection"))

readr::write_csv(diagnostic, file.path("tables", "monocle3_branchlen_diagnostic.csv"))
readr::write_csv(terminal_summary, file.path("tables", "monocle3_branchlen_terminal_summary.csv"))
readr::write_csv(terminal_detection, file.path("tables", "monocle3_branchlen_terminal_detection.csv"))

log_msg("Diagnostic table:", file.path("tables", "monocle3_branchlen_diagnostic.csv"))
log_msg("Terminal summary:", file.path("tables", "monocle3_branchlen_terminal_summary.csv"))
log_msg("Terminal detection detail:", file.path("tables", "monocle3_branchlen_terminal_detection.csv"))
log_msg("Diagnostic rows:")
capture.output(print(diagnostic)) |> paste(collapse = "\n") |> log_msg()
log_msg("Total elapsed seconds:", round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 2))
log_msg("STOP: minimal_branch_len diagnostic complete. Await audit before Phase 2B.")
