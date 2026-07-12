# 05_恶性细胞分亚群与Neftel对照/10e_S3_vs_S4_GSEA.R
# Rule v2 Evidence 3: S3 vs S4 ranked-logFC GSEA using Hallmark and Reactome.

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(dplyr)
  library(readr)
  library(clusterProfiler)
  library(msigdbr)
})

set.seed(42)

params <- list(
  tables_dir = file.path("05_恶性细胞分亚群与Neftel对照", "tables"),
  de_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10d_S3_vs_S4_FindMarkers_all.csv"
  ),
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
  stopifnot(all(c("gene", "avg_log2FC") %in% colnames(de_df)))
  ranks <- de_df |>
    filter(!is.na(gene), !is.na(avg_log2FC), is.finite(avg_log2FC)) |>
    group_by(gene) |>
    summarise(rank_value = avg_log2FC[which.max(abs(avg_log2FC))], .groups = "drop") |>
    arrange(desc(rank_value))
  gene_list <- ranks$rank_value
  names(gene_list) <- ranks$gene
  sort(gene_list, decreasing = TRUE)
}

run_gsea <- function(gene_list, term2gene, database) {
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
  if (nrow(out) == 0) {
    return(tibble())
  }
  out |>
    as_tibble() |>
    mutate(
      comparison = "Subtype3_vs_Subtype4",
      database = database,
      direction = case_when(
        NES > 0 ~ "Subtype3_enriched",
        NES < 0 ~ "Subtype4_enriched",
        TRUE ~ "neutral"
      ),
      .before = 1
    ) |>
    arrange(p.adjust, desc(abs(NES)))
}

msg("Reading DE table:", params$de_file)
de <- read_csv(params$de_file, show_col_types = FALSE)
gene_list <- make_rank(de)
stopifnot(length(gene_list) > 5000)
write_csv(
  tibble(gene = names(gene_list), rank_value_avg_log2FC_S3_vs_S4 = unname(gene_list)),
  file.path(params$tables_dir, "10e_S3_vs_S4_ranked_gene_list.csv")
)

msg("Loading MSigDB")
hallmark <- get_msigdb("H")
reactome <- get_msigdb("C2", "CP:REACTOME")

msg("Running Hallmark GSEA")
gsea_h <- run_gsea(gene_list, hallmark, "Hallmark")
msg("Running Reactome GSEA")
gsea_r <- run_gsea(gene_list, reactome, "Reactome")

write_csv(gsea_h, file.path(params$tables_dir, "10e_S3_vs_S4_GSEA_hallmark.csv"))
write_csv(gsea_r, file.path(params$tables_dir, "10e_S3_vs_S4_GSEA_reactome.csv"))

combined <- bind_rows(gsea_h, gsea_r)
write_csv(combined, file.path(params$tables_dir, "10e_S3_vs_S4_GSEA_combined.csv"))

sig <- combined |>
  filter(!is.na(p.adjust), p.adjust < 0.05)
s3_sig <- sig |> filter(NES > 0)
s4_sig <- sig |> filter(NES < 0)

decision <- tibble(
  evidence = "evidence3_pathway_functional",
  comparison = "Subtype3_vs_Subtype4",
  n_ranked_genes = length(gene_list),
  n_sig_pathways_total = nrow(sig),
  n_S3_enriched_sig_pathways = nrow(s3_sig),
  n_S4_enriched_sig_pathways = nrow(s4_sig),
  top_S3_pathways = paste(head(s3_sig$Description, 5), collapse = "; "),
  top_S4_pathways = paste(head(s4_sig$Description, 5), collapse = "; "),
  evidence3_satisfied = nrow(s3_sig) >= 3 & nrow(s4_sig) >= 3
)
write_csv(decision, file.path(params$tables_dir, "10e_evidence3_decision.csv"))

sanity <- tibble(
  n_de_genes = nrow(de),
  n_ranked_genes = length(gene_list),
  n_hallmark_results = nrow(gsea_h),
  n_reactome_results = nrow(gsea_r),
  n_sig_hallmark = sum(gsea_h$p.adjust < 0.05, na.rm = TRUE),
  n_sig_reactome = sum(gsea_r$p.adjust < 0.05, na.rm = TRUE),
  rank_max = max(gene_list),
  rank_min = min(gene_list)
)
write_csv(sanity, file.path(params$tables_dir, "10e_S3_vs_S4_GSEA_sanity_checks.csv"))

stopifnot(sanity$n_ranked_genes > 5000)
stopifnot(nrow(combined) > 0)

write_session_info(file.path(params$tables_dir, "10e_S3_vs_S4_GSEA_session_info.txt"))
msg("10e S3 vs S4 GSEA completed")
print(decision)
