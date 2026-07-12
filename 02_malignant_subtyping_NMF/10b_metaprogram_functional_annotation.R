# 05_恶性细胞分亚群与Neftel对照/10b_metaprogram_functional_annotation.R
# ORA-based functional annotation for NMF metaprograms.
#
# Main ORA: 10c 50-gene signatures used for UCell scoring.
# Sensitivity ORA: 10b strict consensus genes, variable gene-set sizes.
# Universe: union of per-patient HVGs used in 10a NMF.

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
  library(msigdbr)
})

set.seed(42)

params <- list(
  tables_dir = file.path("05_恶性细胞分亚群与Neftel对照", "tables"),
  hvg_file = file.path("05_恶性细胞分亚群与Neftel对照", "tables", "10a_per_patient_hvg.csv"),
  signature_file = file.path("05_恶性细胞分亚群与Neftel对照", "tables", "10c_metaprogram_signatures_used.csv"),
  metaprogram_summary_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10b_metaprogram_summary.csv"
  ),
  main_mps = c("MP01", "MP02", "MP03", "MP04", "MP05", "MP06"),
  forbidden_terms = c(
    "perivascular", "endothelial", "microglia", "macrophage",
    "pericyte", "smooth muscle", "fibroblast", "lymphocyte"
  )
)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

write_session_info <- function(path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  sink(con)
  print(sessionInfo())
  sink()
}

split_genes <- function(x) {
  genes <- trimws(unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE))
  genes[!is.na(genes) & genes != ""]
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

run_one_enricher <- function(genes, term2gene, universe, mp, gene_input, database) {
  genes <- intersect(unique(genes), universe)
  if (length(genes) < 5) {
    return(tibble(
      metaprogram_id = mp,
      gene_input = gene_input,
      database = database,
      n_genes_input = length(genes),
      ID = NA_character_,
      Description = NA_character_,
      GeneRatio = NA_character_,
      BgRatio = NA_character_,
      pvalue = NA_real_,
      p.adjust = NA_real_,
      qvalue = NA_real_,
      geneID = NA_character_,
      Count = NA_integer_
    ))
  }
  enr <- suppressWarnings(clusterProfiler::enricher(
    gene = genes,
    TERM2GENE = term2gene,
    universe = universe,
    pvalueCutoff = 1,
    qvalueCutoff = 1
  ))
  if (is.null(enr) || nrow(as.data.frame(enr)) == 0) {
    return(tibble(
      metaprogram_id = mp,
      gene_input = gene_input,
      database = database,
      n_genes_input = length(genes),
      ID = NA_character_,
      Description = NA_character_,
      GeneRatio = NA_character_,
      BgRatio = NA_character_,
      pvalue = NA_real_,
      p.adjust = NA_real_,
      qvalue = NA_real_,
      geneID = NA_character_,
      Count = NA_integer_
    ))
  }
  as.data.frame(enr) |>
    arrange(p.adjust, pvalue) |>
    slice_head(n = 10) |>
    as_tibble() |>
    mutate(
      metaprogram_id = mp,
      gene_input = gene_input,
      database = database,
      n_genes_input = length(genes),
      .before = 1
    )
}

clean_pathway_name <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^REACTOME_", "", x)
  x <- gsub("^GOBP_", "", x)
  x <- gsub("^GO_", "", x)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  x <- stringr::str_squish(x)
  x <- tools::toTitleCase(x)
  paste0(x, " program")
}

contains_forbidden <- function(x, forbidden_terms) {
  any(vapply(forbidden_terms, function(term) grepl(term, tolower(x), fixed = TRUE), logical(1)))
}

top_terms <- function(ora_df, mp, db, n = 3) {
  terms <- ora_df |>
    filter(metaprogram_id == mp, database == db, !is.na(Description)) |>
    arrange(p.adjust, pvalue) |>
    pull(Description)
  head(terms, n)
}

make_candidates <- function(ora_main_all, mp, forbidden_terms) {
  terms_ranked <- ora_main_all |>
    filter(metaprogram_id == mp, !is.na(Description)) |>
    arrange(p.adjust, pvalue) |>
    pull(Description) |>
    unique()
  candidates <- vapply(terms_ranked, clean_pathway_name, character(1))
  allowed <- candidates[!vapply(candidates, contains_forbidden, logical(1), forbidden_terms = forbidden_terms)]
  if (length(allowed) < 3) {
    allowed <- unique(c(allowed, paste0(mp, " pathway program")))
  }
  allowed <- head(rep(allowed, length.out = 3), 3)
  allowed
}

msg("Reading inputs")
hvg_tbl <- read_csv(params$hvg_file, show_col_types = FALSE)
stopifnot("hvg_gene" %in% colnames(hvg_tbl))
universe <- unique(hvg_tbl$hvg_gene)
stopifnot(length(universe) > 3000, length(universe) < 30000)

sig_tbl <- read_csv(params$signature_file, show_col_types = FALSE)
main_signatures <- sig_tbl |>
  filter(used_in_object) |>
  mutate(metaprogram_id = if_else(metaprogram == "MP05_cycling_supplement", "MP05", metaprogram)) |>
  group_by(metaprogram_id) |>
  summarise(genes = list(unique(gene)), n_genes = n_distinct(gene), .groups = "drop") |>
  filter(metaprogram_id %in% params$main_mps)
stopifnot(setequal(main_signatures$metaprogram_id, params$main_mps))
stopifnot(all(main_signatures$n_genes == 50))

mp_summary <- read_csv(params$metaprogram_summary_file, show_col_types = FALSE)
sensitivity_signatures <- mp_summary |>
  mutate(
    genes = lapply(top50_consensus_genes, split_genes),
    n_genes = lengths(genes)
  ) |>
  select(metaprogram_id, genes, n_genes, is_recurrent, is_cycling)
stopifnot(all(paste0("MP0", 1:7) %in% sensitivity_signatures$metaprogram_id))

msg("Loading MSigDB TERM2GENE tables")
hallmark <- get_msigdb("H")
reactome <- get_msigdb("C2", "CP:REACTOME")
gobp <- get_msigdb("C5", "GO:BP")
dbs <- list(Hallmark = hallmark, Reactome = reactome, GOBP = gobp)

run_ora_set <- function(signature_tbl, gene_input) {
  rows <- list()
  for (i in seq_len(nrow(signature_tbl))) {
    mp <- signature_tbl$metaprogram_id[i]
    genes <- signature_tbl$genes[[i]]
    for (db_name in names(dbs)) {
      msg("ORA", gene_input, mp, db_name, "n_genes", length(genes))
      rows[[length(rows) + 1]] <- run_one_enricher(
        genes = genes,
        term2gene = dbs[[db_name]],
        universe = universe,
        mp = mp,
        gene_input = gene_input,
        database = db_name
      )
    }
  }
  bind_rows(rows)
}

ora_main <- run_ora_set(main_signatures, "10c_50gene")
ora_sensitivity <- run_ora_set(sensitivity_signatures, "10b_strict_consensus")

write_db <- function(df, db_name, suffix) {
  out <- df |>
    filter(database == db_name)
  write_csv(out, file.path(params$tables_dir, sprintf("10b_metaprogram_ORA_%s_%s.csv", tolower(db_name), suffix)))
}

for (db_name in names(dbs)) {
  write_db(ora_main, db_name, "main")
  write_db(ora_sensitivity, db_name, "sensitivity")
}

sig_any <- bind_rows(ora_main, ora_sensitivity) |>
  filter(!is.na(p.adjust), p.adjust < 0.05) |>
  count(gene_input, metaprogram_id, name = "n_sig_pathways")
write_csv(sig_any, file.path(params$tables_dir, "10b_metaprogram_ORA_significant_pathway_counts.csv"))

naming_rows <- lapply(params$main_mps, function(mp) {
  candidates <- make_candidates(ora_main, mp, params$forbidden_terms)
  tibble(
    metaprogram_id = mp,
    gene_input = "10c_50gene",
    top3_hallmark = paste(top_terms(ora_main, mp, "Hallmark"), collapse = "; "),
    top3_reactome = paste(top_terms(ora_main, mp, "Reactome"), collapse = "; "),
    top3_gobp = paste(top_terms(ora_main, mp, "GOBP"), collapse = "; "),
    naming_candidate_1 = candidates[1],
    naming_candidate_2 = candidates[2],
    naming_candidate_3 = candidates[3],
    rationale_pathways = paste(unique(c(
      top_terms(ora_main, mp, "Hallmark", 2),
      top_terms(ora_main, mp, "Reactome", 2),
      top_terms(ora_main, mp, "GOBP", 2)
    )), collapse = "; ")
  )
}) |>
  bind_rows()

write_csv(naming_rows, file.path(params$tables_dir, "10b_metaprogram_proposed_naming_candidates.csv"))

candidate_long <- naming_rows |>
  pivot_longer(
    cols = starts_with("naming_candidate_"),
    names_to = "candidate_rank",
    values_to = "candidate_name"
  ) |>
  mutate(
    forbidden_hit = vapply(candidate_name, contains_forbidden, logical(1), forbidden_terms = params$forbidden_terms)
  ) |>
  filter(forbidden_hit)
write_csv(candidate_long, file.path(params$tables_dir, "10b_metaprogram_naming_forbidden_term_warnings.csv"))

warnings <- list()
for (mp in params$main_mps) {
  any_sig <- ora_main |>
    filter(metaprogram_id == mp, !is.na(p.adjust), p.adjust < 0.05) |>
    nrow() > 0
  if (!any_sig) warnings[[length(warnings) + 1]] <- sprintf("MP %s has no significant ORA pathway in main input", mp)
}
if (nrow(candidate_long) > 0) {
  warnings[[length(warnings) + 1]] <- sprintf("%d naming candidates contain forbidden terms", nrow(candidate_long))
}
warning_lines <- as.character(unlist(warnings))
if (length(warning_lines) == 0) {
  warning_lines <- character(0)
}
writeLines(warning_lines, file.path(params$tables_dir, "10b_metaprogram_ORA_warnings.txt"), useBytes = TRUE)

sanity <- tibble(
  universe_n_genes = length(universe),
  n_main_metaprograms = n_distinct(main_signatures$metaprogram_id),
  all_main_signature_size_50 = all(main_signatures$n_genes == 50),
  n_sensitivity_metaprograms = n_distinct(sensitivity_signatures$metaprogram_id),
  n_hallmark_terms = n_distinct(hallmark$gs_name),
  n_reactome_terms = n_distinct(reactome$gs_name),
  n_gobp_terms = n_distinct(gobp$gs_name),
  n_forbidden_candidate_warnings = nrow(candidate_long),
  n_warning_lines = length(warning_lines)
)
write_csv(sanity, file.path(params$tables_dir, "10b_metaprogram_ORA_sanity_checks.csv"))

stopifnot(sanity$n_main_metaprograms == 6)
stopifnot(sanity$all_main_signature_size_50)
stopifnot(sanity$n_sensitivity_metaprograms >= 7)
stopifnot(sanity$n_forbidden_candidate_warnings == 0)

write_session_info(file.path(params$tables_dir, "10b_metaprogram_ORA_session_info.txt"))

msg("10b metaprogram functional annotation completed")
print(sanity)
print(naming_rows)
