# 05_恶性细胞分亚群与Neftel对照/10e_one_vs_rest_GSEA_wilcoxon_downsampled.R
# Downsampled formal Wilcoxon one-vs-rest DE and GSEA for Subtype1-4.
# This replaces the sparse-rank overview as the preferred input for subtype-level
# pathway figures while keeping full-cell sparse-rank outputs as sensitivity.

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(qs2)
  library(Seurat)
  library(dplyr)
  library(readr)
  library(clusterProfiler)
  library(msigdbr)
  library(future)
})

set.seed(42)
options(future.globals.maxSize = 8 * 1024^3)
future::plan("sequential")

params <- list(
  input_object = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.v2.qs2"
  ),
  tables_dir = file.path("05_恶性细胞分亚群与Neftel对照", "tables"),
  subtype_levels = paste0("Subtype", 1:4),
  max_cells_per_ident = 1000L,
  logfc_threshold_for_test = 0,
  min_pct_for_test = 0.1,
  min_gs_size = 10L,
  max_gs_size = 500L,
  p_cutoff = 1
)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

write_session_info <- function(path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  sink(con)
  print(sessionInfo())
  sink()
}

get_msigdb <- function(collection, subcollection = NULL) {
  args <- list(species = "Homo sapiens", collection = collection)
  if (!is.null(subcollection)) args$subcollection <- subcollection
  out <- tryCatch(
    do.call(msigdbr::msigdbr, args),
    error = function(e) {
      legacy_args <- list(species = "Homo sapiens", category = collection)
      if (!is.null(subcollection)) legacy_args$subcategory <- subcollection
      do.call(msigdbr::msigdbr, legacy_args)
    }
  )
  out |>
    dplyr::select(gs_name, gene_symbol) |>
    dplyr::distinct()
}

make_rank <- function(de_df, subtype) {
  stopifnot(all(c("gene", "avg_log2FC", "cluster") %in% colnames(de_df)))
  ranks <- de_df |>
    filter(
      cluster == subtype,
      !is.na(gene),
      !is.na(avg_log2FC),
      is.finite(avg_log2FC)
    ) |>
    group_by(gene) |>
    summarise(rank_value = avg_log2FC[which.max(abs(avg_log2FC))], .groups = "drop") |>
    arrange(desc(rank_value))
  gene_list <- ranks$rank_value
  names(gene_list) <- ranks$gene
  sort(gene_list, decreasing = TRUE)
}

run_gsea <- function(gene_list, term2gene, database, subtype) {
  gsea <- suppressWarnings(clusterProfiler::GSEA(
    geneList = gene_list,
    TERM2GENE = term2gene,
    minGSSize = params$min_gs_size,
    maxGSSize = params$max_gs_size,
    pvalueCutoff = params$p_cutoff,
    verbose = FALSE,
    seed = TRUE
  ))
  out <- as.data.frame(gsea)
  if (nrow(out) == 0) return(tibble())
  out |>
    as_tibble() |>
    mutate(
      subtype_k4 = subtype,
      comparison = paste0(subtype, "_vs_rest_downsampled"),
      database = database,
      direction = case_when(
        NES > 0 ~ "subtype_enriched",
        NES < 0 ~ "rest_enriched",
        TRUE ~ "neutral"
      ),
      .before = 1
    ) |>
    arrange(p.adjust, desc(abs(NES)))
}

msg("Loading object:", params$input_object)
obj <- qs2::qs_read(params$input_object)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}
obj@meta.data$subtype_k4 <- factor(as.character(obj@meta.data$subtype_k4), levels = params$subtype_levels)
Idents(obj) <- "subtype_k4"

md <- obj@meta.data
set.seed(42)
cells_use <- unlist(lapply(params$subtype_levels, function(st) {
  cells <- rownames(md)[md$subtype_k4 == st]
  sample(cells, size = min(length(cells), params$max_cells_per_ident))
}), use.names = FALSE)
obj_ds <- subset(obj, cells = cells_use)
Idents(obj_ds) <- "subtype_k4"
ds_counts <- as.data.frame(table(obj_ds@meta.data$subtype_k4))
colnames(ds_counts) <- c("subtype_k4", "n_cells_downsampled")
write_csv(ds_counts, file.path(params$tables_dir, "10e_one_vs_rest_wilcoxon_downsampled_cell_counts.csv"))
stopifnot(all(ds_counts$n_cells_downsampled <= params$max_cells_per_ident))
stopifnot(all(ds_counts$n_cells_downsampled > 0))

msg("Running FindAllMarkers on downsampled object")
markers <- FindAllMarkers(
  obj_ds,
  assay = "RNA",
  only.pos = FALSE,
  logfc.threshold = params$logfc_threshold_for_test,
  min.pct = params$min_pct_for_test,
  test.use = "wilcox",
  max.cells.per.ident = params$max_cells_per_ident,
  random.seed = 42,
  verbose = FALSE
)
if (!"avg_log2FC" %in% colnames(markers) && "avg_logFC" %in% colnames(markers)) {
  markers$avg_log2FC <- markers$avg_logFC
}
markers$BH_q <- p.adjust(markers$p_val, method = "BH")
markers <- markers |>
  rename(subtype_k4 = cluster) |>
  mutate(cluster = subtype_k4, .after = subtype_k4) |>
  arrange(subtype_k4, BH_q, desc(abs(avg_log2FC)))
write_csv(markers, file.path(params$tables_dir, "10e_one_vs_rest_FindAllMarkers_wilcoxon_downsampled.csv"))

msg("Loading MSigDB")
hallmark <- get_msigdb("H")
reactome <- get_msigdb("C2", "CP:REACTOME")

gsea_h_rows <- list()
gsea_r_rows <- list()
summary_rows <- list()

for (subtype in params$subtype_levels) {
  gene_list <- make_rank(markers, subtype)
  write_csv(
    tibble(
      subtype_k4 = subtype,
      gene = names(gene_list),
      rank_value_avg_log2FC_subtype_vs_rest = unname(gene_list)
    ),
    file.path(params$tables_dir, sprintf("10e_%s_vs_rest_wilcoxon_downsampled_ranked_gene_list.csv", subtype))
  )

  msg("GSEA", subtype, "Hallmark")
  gh <- run_gsea(gene_list, hallmark, "Hallmark", subtype)
  msg("GSEA", subtype, "Reactome")
  gr <- run_gsea(gene_list, reactome, "Reactome", subtype)
  gsea_h_rows[[subtype]] <- gh
  gsea_r_rows[[subtype]] <- gr

  sig <- bind_rows(gh, gr) |>
    filter(!is.na(p.adjust), p.adjust < 0.05, NES > 0)
  summary_rows[[subtype]] <- tibble(
    subtype_k4 = subtype,
    n_ranked_genes = length(gene_list),
    n_sig_subtype_enriched_pathways = nrow(sig),
    top_subtype_enriched_pathways = paste(head(sig$Description, 8), collapse = "; ")
  )
}

gsea_h_all <- bind_rows(gsea_h_rows)
gsea_r_all <- bind_rows(gsea_r_rows)
gsea_all <- bind_rows(gsea_h_all, gsea_r_all)
summary_df <- bind_rows(summary_rows) |>
  left_join(ds_counts, by = "subtype_k4")

write_csv(gsea_h_all, file.path(params$tables_dir, "10e_one_vs_rest_GSEA_hallmark_wilcoxon_downsampled.csv"))
write_csv(gsea_r_all, file.path(params$tables_dir, "10e_one_vs_rest_GSEA_reactome_wilcoxon_downsampled.csv"))
write_csv(gsea_all, file.path(params$tables_dir, "10e_one_vs_rest_GSEA_combined_wilcoxon_downsampled.csv"))
write_csv(summary_df, file.path(params$tables_dir, "10e_one_vs_rest_GSEA_summary_wilcoxon_downsampled.csv"))

sanity <- tibble(
  n_downsampled_cells = ncol(obj_ds),
  max_cells_per_ident = params$max_cells_per_ident,
  n_marker_rows = nrow(markers),
  n_hallmark_rows = nrow(gsea_h_all),
  n_reactome_rows = nrow(gsea_r_all),
  all_subtypes_have_markers = all(params$subtype_levels %in% unique(markers$subtype_k4)),
  all_subtypes_have_gsea = all(params$subtype_levels %in% unique(gsea_all$subtype_k4))
)
write_csv(sanity, file.path(params$tables_dir, "10e_one_vs_rest_GSEA_wilcoxon_downsampled_sanity_checks.csv"))

stopifnot(sanity$n_downsampled_cells <= length(params$subtype_levels) * params$max_cells_per_ident)
stopifnot(sanity$all_subtypes_have_markers)
stopifnot(sanity$all_subtypes_have_gsea)

write_session_info(file.path(params$tables_dir, "10e_one_vs_rest_GSEA_wilcoxon_downsampled_session_info.txt"))
msg("10e downsampled Wilcoxon one-vs-rest GSEA completed")
print(summary_df)
