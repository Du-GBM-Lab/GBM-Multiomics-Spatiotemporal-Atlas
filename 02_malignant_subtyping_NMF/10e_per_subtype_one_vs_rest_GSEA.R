# 05_恶性细胞分亚群与Neftel对照/10e_per_subtype_one_vs_rest_GSEA.R
# One-vs-rest RNA-assay sparse average-log2FC ranking and GSEA for Subtype1-4.
# This is a fast pathway overview, not a replacement for formal Wilcoxon DE.

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(qs2)
  library(Seurat)
  library(Matrix)
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
  min_pct_for_rank = 0.1,
  pseudocount = 1,
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

make_rank <- function(de_df) {
  ranks <- de_df |>
    filter(!is.na(gene), !is.na(avg_log2FC), is.finite(avg_log2FC)) |>
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
      comparison = paste0(subtype, "_vs_rest"),
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

compute_sparse_rank <- function(mat, cells_1, cells_2, subtype) {
  x1 <- mat[, cells_1, drop = FALSE]
  x2 <- mat[, cells_2, drop = FALSE]
  pct1 <- Matrix::rowMeans(x1 > 0)
  pct2 <- Matrix::rowMeans(x2 > 0)
  keep <- pct1 >= params$min_pct_for_rank | pct2 >= params$min_pct_for_rank
  x1 <- x1[keep, , drop = FALSE]
  x2 <- x2[keep, , drop = FALSE]
  pct1 <- pct1[keep]
  pct2 <- pct2[keep]
  avg1 <- Matrix::rowMeans(expm1(x1))
  avg2 <- Matrix::rowMeans(expm1(x2))
  avg_log2FC <- log2((avg1 + params$pseudocount) / (avg2 + params$pseudocount))
  tibble(
    subtype_k4 = subtype,
    comparison = paste0(subtype, "_vs_rest"),
    gene = rownames(x1),
    avg_log2FC = as.numeric(avg_log2FC),
    pct.1 = as.numeric(pct1),
    pct.2 = as.numeric(pct2),
    avg_exp_subtype = as.numeric(avg1),
    avg_exp_rest = as.numeric(avg2)
  ) |>
    arrange(desc(avg_log2FC))
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
stopifnot(all(params$subtype_levels %in% levels(Idents(obj))))
rna_data <- GetAssayData(obj, assay = "RNA", layer = "data")
stopifnot(inherits(rna_data, "dgCMatrix"))

msg("Loading MSigDB")
hallmark <- get_msigdb("H")
reactome <- get_msigdb("C2", "CP:REACTOME")

rank_rows <- list()
gsea_h_rows <- list()
gsea_r_rows <- list()
summary_rows <- list()

for (subtype in params$subtype_levels) {
  msg("Sparse average-log2FC rank", subtype, "vs rest")
  cells_1 <- rownames(md)[md$subtype_k4 == subtype]
  cells_2 <- rownames(md)[md$subtype_k4 != subtype]
  mk <- compute_sparse_rank(rna_data, cells_1, cells_2, subtype)
  write_csv(mk, file.path(params$tables_dir, sprintf("10e_%s_vs_rest_sparse_avg_log2FC_rank.csv", subtype)))
  rank_rows[[subtype]] <- mk

  gene_list <- make_rank(mk)
  write_csv(
    tibble(subtype_k4 = subtype, gene = names(gene_list), rank_value_avg_log2FC_subtype_vs_rest = unname(gene_list)),
    file.path(params$tables_dir, sprintf("10e_%s_vs_rest_ranked_gene_list.csv", subtype))
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
    n_cells_subtype = length(cells_1),
    n_cells_rest = length(cells_2),
    n_rank_input_genes = nrow(mk),
    n_ranked_genes = length(gene_list),
    n_sig_subtype_enriched_pathways = nrow(sig),
    top_subtype_enriched_pathways = paste(head(sig$Description, 8), collapse = "; ")
  )
}

rank_all <- bind_rows(rank_rows)
gsea_h_all <- bind_rows(gsea_h_rows)
gsea_r_all <- bind_rows(gsea_r_rows)
gsea_all <- bind_rows(gsea_h_all, gsea_r_all)
summary_df <- bind_rows(summary_rows)

write_csv(rank_all, file.path(params$tables_dir, "10e_one_vs_rest_sparse_avg_log2FC_rank_combined.csv"))
write_csv(gsea_h_all, file.path(params$tables_dir, "10e_one_vs_rest_GSEA_hallmark.csv"))
write_csv(gsea_r_all, file.path(params$tables_dir, "10e_one_vs_rest_GSEA_reactome.csv"))
write_csv(gsea_all, file.path(params$tables_dir, "10e_one_vs_rest_GSEA_combined.csv"))
write_csv(summary_df, file.path(params$tables_dir, "10e_one_vs_rest_GSEA_summary.csv"))

sanity <- tibble(
  n_subtypes = length(params$subtype_levels),
  n_rank_rows = nrow(rank_all),
  n_hallmark_rows = nrow(gsea_h_all),
  n_reactome_rows = nrow(gsea_r_all),
  all_subtypes_have_rank = all(params$subtype_levels %in% unique(rank_all$subtype_k4)),
  all_subtypes_have_gsea = all(params$subtype_levels %in% unique(gsea_all$subtype_k4))
)
write_csv(sanity, file.path(params$tables_dir, "10e_one_vs_rest_GSEA_sanity_checks.csv"))

stopifnot(sanity$n_subtypes == 4)
stopifnot(sanity$all_subtypes_have_rank)
stopifnot(sanity$all_subtypes_have_gsea)

write_session_info(file.path(params$tables_dir, "10e_one_vs_rest_GSEA_session_info.txt"))
msg("10e one-vs-rest GSEA completed")
print(summary_df)
