# 05_恶性细胞分亚群与Neftel对照/10c_metaprogram_subtype_association.R
# Score recurrent NMF metaprograms in all malignant cells and test subtype association.

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(qs2)
  library(Seurat)
  library(UCell)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(effsize)
})

set.seed(42)

params <- list(
  project_dir = "05_恶性细胞分亚群与Neftel对照",
  input_object = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.v2.qs2"
  ),
  tables_dir = file.path("05_恶性细胞分亚群与Neftel对照", "tables"),
  metaprogram_summary_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10b_metaprogram_summary.csv"
  ),
  all_program_top50_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10b_all_programs_top50_genes.csv"
  ),
  assignment_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10b_metaprogram_assignment.csv"
  ),
  subtype_levels = paste0("Subtype", 1:4),
  ucell_max_rank = 1500L
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

frequency_ranked_genes <- function(mp_id, all_program_top50, assignment_df) {
  programs <- assignment_df$program_id[assignment_df$metaprogram_id == mp_id]
  genes <- all_program_top50 |>
    filter(program_id %in% programs) |>
    count(gene, name = "n_programs_with_gene") |>
    arrange(desc(n_programs_with_gene), gene) |>
    pull(gene)
  genes
}

build_signature <- function(mp_id, consensus_string, all_program_top50, assignment_df, row_genes) {
  consensus <- split_genes(consensus_string)
  ranked <- frequency_ranked_genes(mp_id, all_program_top50, assignment_df)
  combined <- unique(c(consensus, ranked))
  genes_used <- intersect(combined, row_genes)
  genes_used <- genes_used[seq_len(min(50, length(genes_used)))]
  list(
    genes = genes_used,
    n_consensus = length(consensus),
    n_used = length(genes_used),
    signature_source = ifelse(length(consensus) >= 20, "consensus_ge20", "consensus_plus_frequency_ranked_topup")
  )
}

cliffs_delta_estimate <- function(x, y) {
  unname(effsize::cliff.delta(x, y)$estimate)
}

fast_permanova_two_group <- function(score_mat, group, n_perm = 999, seed = 42) {
  group <- factor(group)
  if (nlevels(group) != 2) {
    stop("fast_permanova_two_group requires exactly two groups.", call. = FALSE)
  }
  X <- as.matrix(score_mat)
  n <- nrow(X)
  p <- ncol(X)
  g <- nlevels(group)
  overall <- colMeans(X)
  total_ss <- sum(rowSums((X - matrix(overall, n, p, byrow = TRUE))^2))

  between_ss_for_group <- function(grp) {
    sum(vapply(levels(grp), function(lv) {
      idx <- grp == lv
      n_g <- sum(idx)
      mu_g <- colMeans(X[idx, , drop = FALSE])
      n_g * sum((mu_g - overall)^2)
    }, numeric(1)))
  }

  ss_between <- between_ss_for_group(group)
  ss_within <- total_ss - ss_between
  f_obs <- (ss_between / (g - 1)) / (ss_within / (n - g))
  r2 <- ss_between / total_ss

  set.seed(seed)
  f_perm <- replicate(n_perm, {
    grp_perm <- factor(sample(group), levels = levels(group))
    ss_b <- between_ss_for_group(grp_perm)
    ss_w <- total_ss - ss_b
    (ss_b / (g - 1)) / (ss_w / (n - g))
  })
  p_value <- (sum(f_perm >= f_obs) + 1) / (n_perm + 1)

  list(
    n = n,
    p = p,
    n_perm = n_perm,
    df_between = g - 1,
    df_within = n - g,
    ss_between = ss_between,
    ss_within = ss_within,
    total_ss = total_ss,
    F = f_obs,
    R2 = r2,
    p_value = p_value
  )
}

test_one_vs_rest <- function(df, subtype, mp_col) {
  x <- df[[mp_col]][df$subtype_k4 == subtype]
  y <- df[[mp_col]][df$subtype_k4 != subtype]
  wt <- suppressWarnings(wilcox.test(x, y))
  tibble(
    subtype_k4 = subtype,
    metaprogram = mp_col,
    n_subtype = sum(df$subtype_k4 == subtype),
    n_rest = sum(df$subtype_k4 != subtype),
    mean_subtype = mean(x),
    mean_rest = mean(y),
    median_subtype = median(x),
    median_rest = median(y),
    wilcox_p = wt$p.value,
    cliffs_delta = cliffs_delta_estimate(x, y)
  )
}

test_s3_s4 <- function(df, mp_col) {
  s34 <- df |>
    filter(subtype_k4 %in% c("Subtype3", "Subtype4"))
  s3 <- s34[[mp_col]][s34$subtype_k4 == "Subtype3"]
  s4 <- s34[[mp_col]][s34$subtype_k4 == "Subtype4"]
  wt <- suppressWarnings(wilcox.test(s3, s4))
  delta <- cliffs_delta_estimate(s3, s4)
  per_patient_delta <- s34 |>
    group_by(Pt_number) |>
    filter(any(subtype_k4 == "Subtype3") & any(subtype_k4 == "Subtype4")) |>
    summarise(
      delta_S3_minus_S4 = mean(.data[[mp_col]][subtype_k4 == "Subtype3"]) -
        mean(.data[[mp_col]][subtype_k4 == "Subtype4"]),
      .groups = "drop"
    )
  consistency <- mean(sign(per_patient_delta$delta_S3_minus_S4) == sign(delta)) * 100
  tibble(
    metaprogram = mp_col,
    n_S3 = length(s3),
    n_S4 = length(s4),
    mean_S3 = mean(s3),
    mean_S4 = mean(s4),
    median_S3 = median(s3),
    median_S4 = median(s4),
    delta_mean_S3_minus_S4 = mean(s3) - mean(s4),
    wilcox_p = wt$p.value,
    cliffs_delta_S3_vs_S4 = delta,
    within_patient_consistency_pct = consistency,
    n_patients_with_both = nrow(per_patient_delta)
  )
}

msg("Loading object:", params$input_object)
obj <- qs2::qs_read(params$input_object)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

md <- obj@meta.data
stopifnot(all(c("Pt_number", "subtype_k4", "G1S_UCell", "G2M_UCell") %in% colnames(md)))
md$subtype_k4 <- factor(as.character(md$subtype_k4), levels = params$subtype_levels)

mp_summary <- read_csv(params$metaprogram_summary_file, show_col_types = FALSE)
all_program_top50 <- read_csv(params$all_program_top50_file, show_col_types = FALSE)
assignment_df <- read_csv(params$assignment_file, show_col_types = FALSE)
noncycling_recurrent <- mp_summary |>
  filter(is_recurrent, !is_cycling) |>
  arrange(metaprogram_id)
cycling_supp <- mp_summary |>
  filter(metaprogram_id == "MP05")

stopifnot(nrow(noncycling_recurrent) == 5)
stopifnot(nrow(cycling_supp) == 1)

sig_rows <- list()
sig_list <- list()
for (i in seq_len(nrow(noncycling_recurrent))) {
  mp <- noncycling_recurrent$metaprogram_id[i]
  sig <- build_signature(
    mp,
    noncycling_recurrent$top50_consensus_genes[i],
    all_program_top50,
    assignment_df,
    rownames(obj)
  )
  genes_used <- sig$genes
  stopifnot(length(genes_used) >= 20)
  sig_list[[mp]] <- genes_used
  sig_rows[[mp]] <- tibble(
    metaprogram = mp,
    gene = genes_used,
    used_in_object = TRUE,
    signature_role = "recurrent_noncycling",
    n_consensus_genes = sig$n_consensus,
    n_signature_genes_used = sig$n_used,
    signature_source = sig$signature_source
  )
}

mp05_sig <- build_signature(
  "MP05",
  cycling_supp$top50_consensus_genes[1],
  all_program_top50,
  assignment_df,
  rownames(obj)
)
mp05_genes <- mp05_sig$genes
stopifnot(length(mp05_genes) >= 10)
sig_list[["MP05_cycling_supplement"]] <- mp05_genes
sig_rows[["MP05_cycling_supplement"]] <- tibble(
  metaprogram = "MP05_cycling_supplement",
  gene = mp05_genes,
  used_in_object = TRUE,
  signature_role = "cycling_supplement",
  n_consensus_genes = mp05_sig$n_consensus,
  n_signature_genes_used = mp05_sig$n_used,
  signature_source = mp05_sig$signature_source
)

write_csv(bind_rows(sig_rows), file.path(params$tables_dir, "10c_metaprogram_signatures_used.csv"))

msg("Running UCell scoring for ", length(sig_list), " metaprogram signatures.")
obj <- UCell::AddModuleScore_UCell(
  obj,
  features = sig_list,
  maxRank = params$ucell_max_rank,
  ncores = 1
)
md <- obj@meta.data

score_cols <- paste0(names(sig_list), "_UCell")
missing_scores <- setdiff(score_cols, colnames(md))
if (length(missing_scores) > 0) {
  stop("Missing UCell score columns: ", paste(missing_scores, collapse = ", "), call. = FALSE)
}

rename_map <- setNames(names(sig_list), score_cols)
for (old in names(rename_map)) {
  md[[rename_map[[old]]]] <- md[[old]]
}

noncycling_cols <- noncycling_recurrent$metaprogram_id
cycling_col <- "MP05_cycling_supplement"
all_mp_cols <- c(noncycling_cols, cycling_col)

out <- md |>
  mutate(
    cell_id = rownames(md),
    subtype_k4 = factor(as.character(subtype_k4), levels = params$subtype_levels)
  ) |>
  select(
    cell_id,
    Pt_number,
    subtype_k4,
    G1S_UCell,
    G2M_UCell,
    all_of(all_mp_cols)
  )

write_tsv(out, file.path(params$tables_dir, "10c_per_cell_metaprogram_scores.tsv"))

score_summary <- out |>
  pivot_longer(cols = all_of(all_mp_cols), names_to = "metaprogram", values_to = "score") |>
  group_by(subtype_k4, metaprogram) |>
  summarise(
    n_cells = n(),
    mean = mean(score),
    median = median(score),
    sd = sd(score),
    q25 = quantile(score, 0.25, names = FALSE),
    q75 = quantile(score, 0.75, names = FALSE),
    iqr = IQR(score),
    .groups = "drop"
  )
write_csv(score_summary, file.path(params$tables_dir, "10c_metaprogram_score_by_subtype.csv"))

one_vs_rest <- bind_rows(lapply(params$subtype_levels, function(st) {
  bind_rows(lapply(noncycling_cols, function(mp) test_one_vs_rest(out, st, mp)))
})) |>
  mutate(wilcox_BH_q = p.adjust(wilcox_p, method = "BH"))
write_csv(one_vs_rest, file.path(params$tables_dir, "10c_one_vs_rest_subtype_enrichment.csv"))

s3_s4 <- bind_rows(lapply(noncycling_cols, function(mp) test_s3_s4(out, mp))) |>
  mutate(wilcox_BH_q = p.adjust(wilcox_p, method = "BH")) |>
  arrange(desc(abs(cliffs_delta_S3_vs_S4)))
write_csv(s3_s4, file.path(params$tables_dir, "10c_S3_vs_S4_metaprogram_test.csv"))

s34 <- out |>
  filter(subtype_k4 %in% c("Subtype3", "Subtype4")) |>
  droplevels()
score_mat <- as.data.frame(s34[, noncycling_cols, drop = FALSE])
permanova <- fast_permanova_two_group(score_mat, s34$subtype_k4, n_perm = 999, seed = 42)
permanova_text <- c(
  "Fast PERMANOVA equivalent for two groups using Euclidean sums of squares",
  "Model: cell x recurrent non-cycling metaprogram score matrix ~ subtype_k4",
  "Groups: Subtype3 vs Subtype4",
  "Permutations: 999",
  sprintf("n = %s", permanova$n),
  sprintf("p = %s", permanova$p),
  sprintf("df_between = %s", permanova$df_between),
  sprintf("df_within = %s", permanova$df_within),
  sprintf("SS_between = %.10f", permanova$ss_between),
  sprintf("SS_within = %.10f", permanova$ss_within),
  sprintf("F = %.10f", permanova$F),
  sprintf("R2 = %.10f", permanova$R2),
  sprintf("p = %.10g", permanova$p_value)
)
writeLines(permanova_text, file.path(params$tables_dir, "10c_PERMANOVA_S3_vs_S4.txt"))

permanova_p <- permanova$p_value
univariate_pass <- s3_s4 |>
  filter(abs(cliffs_delta_S3_vs_S4) >= 0.15, wilcox_BH_q < 0.05)
evidence1_satisfied <- isTRUE(permanova_p < 0.05) && nrow(univariate_pass) >= 1

evidence1 <- tibble(
  evidence = "evidence1_NMF_program_level",
  PERMANOVA_p = permanova_p,
  any_univariate_passes = nrow(univariate_pass) >= 1,
  evidence_satisfied = evidence1_satisfied,
  metaprograms_passing_univariate = paste(univariate_pass$metaprogram, collapse = ",")
)
write_csv(evidence1, file.path(params$tables_dir, "10c_evidence1_decision.csv"))

cycling_summary <- out |>
  group_by(subtype_k4) |>
  summarise(
    n_cells = n(),
    G1S_UCell_mean = mean(G1S_UCell),
    G2M_UCell_mean = mean(G2M_UCell),
    cycling_pmax_UCell_mean = mean(pmax(G1S_UCell, G2M_UCell)),
    MP05_cycling_supplement_mean = mean(MP05_cycling_supplement),
    MP05_cycling_supplement_median = median(MP05_cycling_supplement),
    .groups = "drop"
  )
write_csv(cycling_summary, file.path(params$tables_dir, "10c_MP05_cycling_supplement_by_subtype.csv"))

sanity <- tibble(
  check = c(
    "n_cells",
    "n_noncycling_metaprograms",
    "cycling_supplement_present",
    "score_na_total",
    "score_min",
    "score_max",
    "S3_cells",
    "S4_cells"
  ),
  value = c(
    nrow(out),
    length(noncycling_cols),
    cycling_col %in% colnames(out),
    sum(is.na(out[, all_mp_cols])),
    min(as.matrix(out[, all_mp_cols])),
    max(as.matrix(out[, all_mp_cols])),
    sum(out$subtype_k4 == "Subtype3"),
    sum(out$subtype_k4 == "Subtype4")
  )
)
write_csv(sanity, file.path(params$tables_dir, "10c_sanity_checks.csv"))
write_session_info(file.path(params$tables_dir, "10c_session_info.txt"))

cat("\nSanity checks:\n")
print(sanity)
cat("\nMetaprogram score by subtype:\n")
print(score_summary)
cat("\nS3 vs S4 metaprogram test:\n")
print(s3_s4)
cat("\nPERMANOVA:\n")
cat(paste(permanova_text, collapse = "\n"), "\n")
cat("\nEvidence 1 decision:\n")
print(evidence1)
cat("\nCycling supplement summary:\n")
print(cycling_summary)
msg("Done.")
