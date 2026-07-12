# 05_恶性细胞分亚群与Neftel对照/10b_program_to_metaprogram.R
# Collapse per-patient NMF programs into cross-patient metaprograms.
#
# Fixed decisions before execution:
# - Use only patient-level chosen k from 10a.
# - Align the 5 seed-specific W matrices per patient via Hungarian matching
#   on Pearson correlation between W columns, then average aligned W.
# - Select top 50 genes per program by raw W loading, not row-scaled W.
# - Program similarity = Jaccard similarity of top-50 gene sets.
# - Hierarchical clustering = average linkage on 1 - Jaccard.
# - Candidate metaprogram counts = 6:20; choose max silhouette, ties within
#   0.02 choose smaller k by parsimony.
# - Recurrent metaprogram = contributed by >=6 patients.
# - Cycling metaprogram = max Jaccard overlap with Neftel G1S/G2M > 0.2 OR
#   the maximal Neftel overlap belongs to G1S/G2M.

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(clue)
  library(cluster)
})

set.seed(42)

params <- list(
  project_dir = "05_恶性细胞分亚群与Neftel对照",
  nmf_dir = file.path("05_恶性细胞分亚群与Neftel对照", "outputs", "nmf", "per_patient"),
  output_dir = file.path("05_恶性细胞分亚群与Neftel对照", "outputs", "nmf"),
  tables_dir = file.path("05_恶性细胞分亚群与Neftel对照", "tables"),
  chosen_k_file = file.path("05_恶性细胞分亚群与Neftel对照", "tables", "10a_per_patient_chosen_k.csv"),
  signature_xlsx = file.path("05_恶性细胞分亚群与Neftel对照", "data", "Neftel2019_TableS2.xlsx"),
  seed_values = 1:5,
  top_n_genes = 50L,
  candidate_metaprogram_k = 6:20,
  recurrence_patient_n = 6L,
  cycling_overlap_threshold = 0.2
)

dir.create(params$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(params$tables_dir, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

write_session_info <- function(path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  sink(con)
  print(sessionInfo())
  sink()
}

read_neftel_signatures <- function(path) {
  raw <- readxl::read_excel(path, sheet = "Table S2", skip = 3, col_types = "text")
  names(raw) <- gsub("\\s+", "", names(raw))
  names(raw) <- gsub("/", "", names(raw))
  expected <- c("MES2", "MES1", "AC", "OPC", "NPC1", "NPC2", "G1S", "G2M")
  missing <- setdiff(expected, names(raw))
  if (length(missing) > 0) {
    stop("Missing signature columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  sigs <- lapply(expected, function(module) {
    genes <- trimws(as.character(raw[[module]]))
    unique(genes[!is.na(genes) & genes != ""])
  })
  names(sigs) <- expected
  sigs[c("MES1", "MES2", "AC", "OPC", "NPC1", "NPC2", "G1S", "G2M")]
}

align_W_to_ref <- function(W_target, W_ref) {
  cor_mat <- suppressWarnings(cor(W_ref, W_target, use = "pairwise.complete.obs"))
  cor_mat[is.na(cor_mat)] <- -1
  assignment <- as.integer(clue::solve_LSAP(cor_mat, maximum = TRUE))
  W_target[, assignment, drop = FALSE]
}

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

select_metaprogram_k <- function(selection_df) {
  max_sil <- max(selection_df$mean_silhouette, na.rm = TRUE)
  selection_df |>
    filter(mean_silhouette >= max_sil - 0.02) |>
    arrange(candidate_k) |>
    slice(1) |>
    pull(candidate_k)
}

msg("Reading chosen k:", params$chosen_k_file)
chosen_k <- read_csv(params$chosen_k_file, show_col_types = FALSE) |>
  mutate(patient = as.character(patient), chosen_k = as.integer(chosen_k)) |>
  arrange(patient)

neftel_sigs <- read_neftel_signatures(params$signature_xlsx)

program_top_rows <- list()
program_gene_sets <- list()
program_meta_rows <- list()

msg("Aligning seed-level programs and extracting top genes.")
for (i in seq_len(nrow(chosen_k))) {
  pt <- chosen_k$patient[i]
  k <- chosen_k$chosen_k[i]
  files <- file.path(params$nmf_dir, sprintf("%s_k%d_seed%d_nmf.rds", pt, k, params$seed_values))
  if (!all(file.exists(files))) {
    stop("Missing seed files for ", pt, " k=", k, call. = FALSE)
  }

  seed_objs <- lapply(files, readRDS)
  W_list <- lapply(seed_objs, function(x) {
    W <- x$W
    rownames(W) <- x$hvg
    colnames(W) <- paste0("program_", seq_len(ncol(W)))
    W
  })

  W_ref <- W_list[[1]]
  aligned <- list(W_ref)
  if (length(W_list) > 1) {
    for (s in 2:length(W_list)) {
      aligned[[s]] <- align_W_to_ref(W_list[[s]], W_ref)
    }
  }
  W_avg <- Reduce("+", aligned) / length(aligned)
  rownames(W_avg) <- rownames(W_ref)
  colnames(W_avg) <- paste0("program_", seq_len(k))

  for (program_idx in seq_len(k)) {
    program_id <- sprintf("%s_MP%02d", pt, program_idx)
    loading <- W_avg[, program_idx]
    ord <- order(loading, decreasing = TRUE)
    top_idx <- ord[seq_len(min(params$top_n_genes, length(ord)))]
    top_genes <- rownames(W_avg)[top_idx]
    top_loadings <- loading[top_idx]

    program_gene_sets[[program_id]] <- top_genes
    program_meta_rows[[program_id]] <- tibble(
      program_id = program_id,
      patient = pt,
      chosen_k = k,
      program_idx_within_patient = program_idx
    )
    program_top_rows[[program_id]] <- tibble(
      program_id = program_id,
      patient = pt,
      chosen_k = k,
      program_idx_within_patient = program_idx,
      gene_rank = seq_along(top_genes),
      gene = top_genes,
      loading = as.numeric(top_loadings)
    )
  }
  msg("Processed ", pt, " chosen_k=", k)
}

all_programs_top50 <- bind_rows(program_top_rows)
program_meta <- bind_rows(program_meta_rows)
write_csv(all_programs_top50, file.path(params$tables_dir, "10b_all_programs_top50_genes.csv"))

n_prog <- length(program_gene_sets)
expected_programs <- sum(chosen_k$chosen_k)
stopifnot(n_prog == expected_programs)

msg("Computing ", n_prog, " x ", n_prog, " Jaccard matrix.")
J <- matrix(0, n_prog, n_prog, dimnames = list(names(program_gene_sets), names(program_gene_sets)))
for (i in seq_len(n_prog)) {
  for (j in i:n_prog) {
    J[i, j] <- jaccard(program_gene_sets[[i]], program_gene_sets[[j]])
    J[j, i] <- J[i, j]
  }
}
diag(J) <- 1
saveRDS(J, file.path(params$output_dir, "10b_jaccard_matrix.rds"))

d <- as.dist(1 - J)
hc <- hclust(d, method = "average")
saveRDS(hc, file.path(params$output_dir, "10b_program_dendrogram.rds"))

selection_rows <- lapply(params$candidate_metaprogram_k, function(k) {
  cuts <- cutree(hc, k = k)
  sil <- mean(cluster::silhouette(cuts, d)[, 3])
  tibble(
    candidate_k = k,
    mean_silhouette = sil,
    gap_stat = NA_real_,
    gap_stat_note = "not_computed_for_jaccard_distance",
    n_clusters = length(unique(cuts))
  )
})
selection_df <- bind_rows(selection_rows)
best_k <- select_metaprogram_k(selection_df)
selection_df <- selection_df |>
  mutate(
    selected = candidate_k == best_k,
    selection_rule = "max_silhouette_tie_within_0.02_choose_smaller_k"
  )
write_csv(selection_df, file.path(params$tables_dir, "10b_metaprogram_cluster_selection.csv"))

assignment <- cutree(hc, k = best_k)
metaprogram_ids <- sprintf("MP%02d", assignment)
assignment_df <- tibble(
  program_id = names(assignment),
  metaprogram_id = metaprogram_ids
) |>
  left_join(program_meta, by = "program_id") |>
  arrange(metaprogram_id, patient, program_idx_within_patient)
write_csv(assignment_df, file.path(params$tables_dir, "10b_metaprogram_assignment.csv"))

freq_table_for_mp <- function(program_ids) {
  genes <- unlist(program_gene_sets[program_ids], use.names = FALSE)
  sort(table(genes), decreasing = TRUE)
}

overlap_row <- function(consensus_genes, top10_genes) {
  # Use top10 fallback if consensus is empty.
  genes_for_overlap <- if (length(consensus_genes) > 0) consensus_genes else top10_genes
  vals <- vapply(neftel_sigs, function(sig) jaccard(genes_for_overlap, sig), numeric(1))
  as.list(vals)
}

summary_rows <- list()
for (mp in sort(unique(assignment_df$metaprogram_id))) {
  rows <- assignment_df |> filter(metaprogram_id == mp)
  program_ids <- rows$program_id
  freq <- freq_table_for_mp(program_ids)
  n_programs <- length(program_ids)
  n_patients <- n_distinct(rows$patient)
  top10 <- names(freq)[seq_len(min(10, length(freq)))]
  consensus <- names(freq)[freq >= ceiling(0.5 * n_programs)]
  overlap_vals <- overlap_row(consensus, top10)
  overlap_numeric <- unlist(overlap_vals)
  max_overlap_module <- names(overlap_numeric)[which.max(overlap_numeric)]
  is_cycling <- isTRUE(max(overlap_numeric[c("G1S", "G2M")], na.rm = TRUE) > params$cycling_overlap_threshold) ||
    max_overlap_module %in% c("G1S", "G2M")

  summary_rows[[mp]] <- tibble(
    metaprogram_id = mp,
    n_programs = n_programs,
    n_patients_unique = n_patients,
    recurrence_pct = round(100 * n_patients / nrow(chosen_k), 2),
    top10_signature_genes = paste(top10, collapse = ","),
    top50_consensus_genes = paste(consensus, collapse = ","),
    n_top50_consensus_genes = length(consensus),
    max_neftel_overlap_module = max_overlap_module,
    is_cycling = is_cycling,
    is_recurrent = n_patients >= params$recurrence_patient_n
  ) |>
    bind_cols(as_tibble(overlap_vals) |> rename_with(~ paste0("neftel_overlap_", .x)))
}

meta_summary <- bind_rows(summary_rows) |>
  arrange(metaprogram_id)
write_csv(meta_summary, file.path(params$tables_dir, "10b_metaprogram_summary.csv"))

recurrent <- meta_summary |>
  filter(is_recurrent)
write_csv(recurrent, file.path(params$tables_dir, "10b_recurrent_metaprograms.csv"))

excluded_programs <- assignment_df |>
  left_join(meta_summary |> select(metaprogram_id, is_recurrent), by = "metaprogram_id") |>
  filter(!is_recurrent) |>
  arrange(patient, metaprogram_id, program_idx_within_patient)
write_csv(excluded_programs, file.path(params$tables_dir, "10b_excluded_patient_specific_programs.csv"))

excluded_by_patient <- excluded_programs |>
  count(patient, name = "n_excluded_programs") |>
  arrange(desc(n_excluded_programs), patient)
write_csv(excluded_by_patient, file.path(params$tables_dir, "10b_excluded_programs_by_patient.csv"))

pt9_recurrent_count <- assignment_df |>
  filter(patient == "Pt9") |>
  left_join(meta_summary |> select(metaprogram_id, is_recurrent), by = "metaprogram_id") |>
  summarise(n_pt9_programs_in_recurrent = sum(is_recurrent), n_pt9_programs_total = n()) 
write_csv(pt9_recurrent_count, file.path(params$tables_dir, "10b_pt9_recurrent_program_count.csv"))

pt9_cophenetic <- read_csv(
  file.path(params$tables_dir, "10a_per_patient_consensus_metrics.csv"),
  show_col_types = FALSE
) |>
  filter(patient == "Pt9") |>
  arrange(k)
write_csv(pt9_cophenetic, file.path(params$tables_dir, "10b_pt9_cophenetic_curve.csv"))

stopifnot(length(assignment) == expected_programs)
if (nrow(recurrent) < 6 || nrow(recurrent) > 15) {
  stop("Recurrent metaprogram count outside expected range 6-15: ", nrow(recurrent), call. = FALSE)
}
if (sum(meta_summary$is_cycling) < 1) {
  stop("No cycling metaprogram detected by fixed rule.", call. = FALSE)
}

write_session_info(file.path(params$tables_dir, "10b_session_info.txt"))

cat("\nCluster selection:\n")
print(selection_df)
cat("\nMetaprogram summary:\n")
print(meta_summary)
cat("\nRecurrent metaprograms:\n")
print(recurrent)
cat("\nPt9 recurrent program count:\n")
print(pt9_recurrent_count)
cat("\nExcluded programs by patient:\n")
print(excluded_by_patient)
msg("Done.")
