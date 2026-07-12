# 05_恶性细胞分亚群与Neftel对照/10b_post_metaprogram_cell_audit.R
# Audit MP02-high and MP06-high malignant cells for CNV, doublet, glial-lineage,
# QC, patient concentration, and vascular/endothelial marker profiles.

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
  library(tidyr)
})

set.seed(42)

params <- list(
  input_object = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.v2.qs2"
  ),
  metaprogram_scores_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10c_per_cell_metaprogram_scores.tsv"
  ),
  tables_dir = file.path("05_恶性细胞分亚群与Neftel对照", "tables"),
  outputs_dir = file.path("05_恶性细胞分亚群与Neftel对照", "outputs"),
  mp_cols = c("MP02", "MP04", "MP06"),
  top_quantile = 0.90,
  cnv_score_col = "infercnv_immune_cnv_burden",
  cnv_correlation_col = "infercnv_immune_cnv_correlation",
  cnv_call_col = "infercnv_immune_call",
  malignant_status_col = "malignant_call_status_immune",
  malignant_bool_col = "is_malignant_for_downstream_immune",
  doublet_score_col = "doubletfinder_pANN",
  doublet_class_col = "doubletfinder_class",
  glial_markers = c("GFAP", "SOX2", "OLIG2", "S100B"),
  endothelial_markers = c("CLDN5", "FLT1", "PECAM1", "ENG", "PODXL", "A2M"),
  mural_markers = c("DCN", "COL3A1", "COL1A2", "FN1", "ACTA2", "TAGLN", "SPARCL1"),
  myeloid_markers = c(
    "C1QA", "C1QB", "C1QC", "TYROBP", "CD14", "FCGR3A",
    "CCL3", "CCL4", "CCL3L3", "CCL4L2", "MS4A7", "LAPTM5",
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1"
  ),
  neftel_cols = c(
    "MES1_UCell", "MES2_UCell", "AC_UCell", "OPC_UCell",
    "NPC1_UCell", "NPC2_UCell", "G1S_UCell", "G2M_UCell"
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

mean_sd <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), median = median(x, na.rm = TRUE))
}

safe_frac <- function(x) {
  if (length(x) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

extract_marker_expr <- function(obj, genes, cells) {
  genes_present <- intersect(genes, rownames(obj))
  out <- matrix(NA_real_, nrow = length(cells), ncol = length(genes), dimnames = list(cells, genes))
  if (length(genes_present) == 0) return(as.data.frame(out))
  expr <- GetAssayData(obj, assay = "RNA", layer = "data")[genes_present, cells, drop = FALSE]
  expr_df <- as.data.frame(as.matrix(t(expr)))
  out[, genes_present] <- as.matrix(expr_df[, genes_present, drop = FALSE])
  as.data.frame(out)
}

summarise_high_rest_numeric <- function(df, high_flag, cols, mp_name) {
  rows <- list()
  for (col in cols) {
    if (!col %in% colnames(df)) next
    high <- suppressWarnings(as.numeric(df[[col]][high_flag]))
    rest <- suppressWarnings(as.numeric(df[[col]][!high_flag]))
    rows[[length(rows) + 1]] <- tibble(
      metaprogram = mp_name,
      metric = col,
      high_mean = mean(high, na.rm = TRUE),
      high_sd = sd(high, na.rm = TRUE),
      high_median = median(high, na.rm = TRUE),
      rest_mean = mean(rest, na.rm = TRUE),
      rest_sd = sd(rest, na.rm = TRUE),
      rest_median = median(rest, na.rm = TRUE),
      high_minus_rest_mean = high_mean - rest_mean
    )
  }
  bind_rows(rows)
}

top_patient_summary <- function(df, high_flag, mp_name) {
  high_pt <- df |>
    filter(high_flag) |>
    count(Pt_number, name = "n_high") |>
    arrange(desc(n_high), Pt_number) |>
    mutate(frac_high = n_high / sum(n_high))
  tibble(
    metaprogram = mp_name,
    top1_patient = high_pt$Pt_number[1],
    top1_frac = high_pt$frac_high[1],
    top2_patient = paste(high_pt$Pt_number[seq_len(min(2, nrow(high_pt)))], collapse = ","),
    top2_frac = sum(high_pt$frac_high[seq_len(min(2, nrow(high_pt)))]),
    n_patients_high = nrow(high_pt)
  )
}

subtype_summary <- function(df, high_flag, mp_name) {
  df |>
    filter(high_flag) |>
    count(subtype_k4, name = "n_high") |>
    mutate(
      metaprogram = mp_name,
      frac_high = n_high / sum(n_high),
      .before = 1
    )
}

decide_disposition <- function(mp_name, summary_row, overall_stats, is_mp06 = FALSE) {
  cnv_available <- !is.na(summary_row$cnv_score_high_mean)
  pANN_available <- !is.na(summary_row$doublet_pANN_high_mean)

  cnv_low <- cnv_available &&
    summary_row$cnv_score_high_mean < (overall_stats$cnv_mean - overall_stats$cnv_sd)
  glial_low <- !is.na(summary_row$glial_max_high_mean) &&
    summary_row$glial_max_high_mean < 0.1
  doublet_high <- pANN_available &&
    summary_row$doublet_pANN_high_mean > overall_stats$pANN_95th
  patient_concentrated <- !is.na(summary_row$top2_patient_frac) &&
    summary_row$top2_patient_frac > 0.70

  if (is_mp06 && cnv_low && glial_low) {
    disposition <- "likely_contamination_flag"
    reason <- "CNV score lower than overall mean-1SD and glial marker retention low"
  } else if (doublet_high) {
    disposition <- "likely_doublet_flag"
    reason <- "DoubletFinder pANN mean exceeds overall malignant 95th percentile"
  } else if (patient_concentrated) {
    disposition <- "patient_specific_concern_flag"
    reason <- "Top two patients contribute >70% of high cells"
  } else if (!cnv_available || !pANN_available) {
    disposition <- "retain_with_available_evidence_limitations"
    reason <- "Key CNV or doublet metrics unavailable; no automatic removal rule triggered"
  } else {
    disposition <- "retain_with_biology_note"
    reason <- "CNV/doublet/patient-concentration removal flags not triggered"
  }

  tibble(
    metaprogram = mp_name,
    disposition = disposition,
    reason = reason,
    cnv_low_flag = cnv_low,
    glial_low_flag = glial_low,
    doublet_high_flag = doublet_high,
    patient_concentrated_flag = patient_concentrated
  )
}

msg("Loading object:", params$input_object)
obj <- qs2::qs_read(params$input_object)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

md <- obj@meta.data
md$cell_id <- rownames(md)
md$subtype_k4 <- factor(as.character(md$subtype_k4), levels = paste0("Subtype", 1:4))

mp_scores <- read_tsv(params$metaprogram_scores_file, show_col_types = FALSE)
stopifnot(all(c("cell_id", params$mp_cols) %in% colnames(mp_scores)))
stopifnot(identical(sort(md$cell_id), sort(mp_scores$cell_id)))

md <- md |>
  left_join(mp_scores[, c("cell_id", params$mp_cols)], by = "cell_id")

marker_genes <- unique(c(
  params$glial_markers,
  params$endothelial_markers,
  params$mural_markers,
  params$myeloid_markers
))
marker_expr <- extract_marker_expr(obj, marker_genes, md$cell_id)
marker_expr$cell_id <- rownames(marker_expr)
md <- md |>
  left_join(marker_expr, by = "cell_id")

glial_present <- intersect(params$glial_markers, colnames(md))
endo_present <- intersect(params$endothelial_markers, colnames(md))
mural_present <- intersect(params$mural_markers, colnames(md))
myeloid_present <- intersect(params$myeloid_markers, colnames(md))
md$glial_marker_max <- apply(md[, glial_present, drop = FALSE], 1, max, na.rm = TRUE)
md$glial_marker_mean <- rowMeans(md[, glial_present, drop = FALSE], na.rm = TRUE)
md$endothelial_marker_max <- apply(md[, endo_present, drop = FALSE], 1, max, na.rm = TRUE)
md$endothelial_marker_mean <- rowMeans(md[, endo_present, drop = FALSE], na.rm = TRUE)
md$mural_marker_max <- apply(md[, mural_present, drop = FALSE], 1, max, na.rm = TRUE)
md$mural_marker_mean <- rowMeans(md[, mural_present, drop = FALSE], na.rm = TRUE)
md$myeloid_marker_max <- apply(md[, myeloid_present, drop = FALSE], 1, max, na.rm = TRUE)
md$myeloid_marker_mean <- rowMeans(md[, myeloid_present, drop = FALSE], na.rm = TRUE)

numeric_metrics <- unique(c(
  params$cnv_score_col,
  params$cnv_correlation_col,
  params$doublet_score_col,
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt",
  "glial_marker_max",
  "glial_marker_mean",
  "endothelial_marker_max",
  "endothelial_marker_mean",
  "mural_marker_max",
  "mural_marker_mean",
  "myeloid_marker_max",
  "myeloid_marker_mean",
  params$neftel_cols,
  params$glial_markers,
  params$endothelial_markers,
  params$mural_markers,
  params$myeloid_markers
))
numeric_metrics <- intersect(numeric_metrics, colnames(md))

overall_stats <- tibble(
  cnv_mean = if (params$cnv_score_col %in% colnames(md)) mean(md[[params$cnv_score_col]], na.rm = TRUE) else NA_real_,
  cnv_sd = if (params$cnv_score_col %in% colnames(md)) sd(md[[params$cnv_score_col]], na.rm = TRUE) else NA_real_,
  pANN_95th = if (params$doublet_score_col %in% colnames(md)) quantile(md[[params$doublet_score_col]], 0.95, na.rm = TRUE, names = FALSE) else NA_real_
)

summary_rows <- list()
decision_rows <- list()
profile_files <- character()

for (mp in params$mp_cols) {
  msg("Auditing", mp)
  threshold <- quantile(md[[mp]], params$top_quantile, na.rm = TRUE, names = FALSE)
  high_flag <- md[[mp]] >= threshold
  high_md <- md[high_flag, , drop = FALSE]

  profile_cols <- unique(c(
    "cell_id", "Pt_number", "subtype_k4", mp,
    params$cnv_score_col, params$cnv_correlation_col, params$cnv_call_col,
    params$malignant_status_col, params$malignant_bool_col,
    params$doublet_score_col, params$doublet_class_col,
    "nCount_RNA", "nFeature_RNA", "percent.mt",
    "glial_marker_max", "glial_marker_mean",
    "endothelial_marker_max", "endothelial_marker_mean",
    "mural_marker_max", "mural_marker_mean",
    "myeloid_marker_max", "myeloid_marker_mean",
    params$neftel_cols,
    marker_genes
  ))
  profile_cols <- intersect(profile_cols, colnames(high_md))
  profile <- high_md[, profile_cols, drop = FALSE] |>
    arrange(desc(.data[[mp]]))

  profile_file <- file.path(params$tables_dir, sprintf("10b_audit_%shigh_cell_profile.csv", mp))
  write_csv(profile, profile_file)
  profile_files <- c(profile_files, profile_file)

  if (mp == "MP06") {
    write_csv(
      profile,
      file.path(params$outputs_dir, "10b_audit_MP06_high_cells_for_review.csv")
    )
  }

  metric_summary <- summarise_high_rest_numeric(md, high_flag, numeric_metrics, mp)
  write_csv(metric_summary, file.path(params$tables_dir, sprintf("10b_audit_%s_metric_high_vs_rest.csv", mp)))

  patient_s <- top_patient_summary(md, high_flag, mp)
  subtype_s <- subtype_summary(md, high_flag, mp)
  write_csv(subtype_s, file.path(params$tables_dir, sprintf("10b_audit_%s_subtype_distribution.csv", mp)))

  cnv_vals <- mean_sd(high_md[[params$cnv_score_col]])
  cnv_rest <- mean_sd(md[[params$cnv_score_col]][!high_flag])
  pANN_vals <- mean_sd(high_md[[params$doublet_score_col]])
  pANN_rest <- mean_sd(md[[params$doublet_score_col]][!high_flag])
  glial_vals <- mean_sd(high_md$glial_marker_max)
  glial_rest <- mean_sd(md$glial_marker_max[!high_flag])
  endo_vals <- mean_sd(high_md$endothelial_marker_max)
  mural_vals <- mean_sd(high_md$mural_marker_max)
  myeloid_vals <- mean_sd(high_md$myeloid_marker_max)

  summary_row <- tibble(
    metaprogram = mp,
    score_threshold_top10 = threshold,
    n_high = nrow(high_md),
    n_rest = sum(!high_flag),
    cnv_score_high_mean = cnv_vals["mean"],
    cnv_score_high_sd = cnv_vals["sd"],
    cnv_score_rest_mean = cnv_rest["mean"],
    cnv_score_rest_sd = cnv_rest["sd"],
    cnv_correlation_high_mean = mean(high_md[[params$cnv_correlation_col]], na.rm = TRUE),
    cnv_correlation_rest_mean = mean(md[[params$cnv_correlation_col]][!high_flag], na.rm = TRUE),
    doublet_pANN_high_mean = pANN_vals["mean"],
    doublet_pANN_high_sd = pANN_vals["sd"],
    doublet_pANN_rest_mean = pANN_rest["mean"],
    doublet_pANN_rest_sd = pANN_rest["sd"],
    doublet_pANN_overall_95th = overall_stats$pANN_95th,
    nFeature_high_mean = mean(high_md$nFeature_RNA, na.rm = TRUE),
    nFeature_rest_mean = mean(md$nFeature_RNA[!high_flag], na.rm = TRUE),
    percent_mt_high_mean = mean(high_md$percent.mt, na.rm = TRUE),
    percent_mt_rest_mean = mean(md$percent.mt[!high_flag], na.rm = TRUE),
    glial_max_high_mean = glial_vals["mean"],
    glial_max_rest_mean = glial_rest["mean"],
    endothelial_max_high_mean = endo_vals["mean"],
    mural_max_high_mean = mural_vals["mean"],
    myeloid_max_high_mean = myeloid_vals["mean"],
    malignant_downstream_high_frac = if (params$malignant_bool_col %in% colnames(high_md)) {
      mean(as.logical(high_md[[params$malignant_bool_col]]), na.rm = TRUE)
    } else {
      NA_real_
    }
  ) |>
    bind_cols(patient_s |> select(top1_patient, top1_frac, top2_patient, top2_frac, n_patients_high)) |>
    rename(top2_patient_frac = top2_frac)

  summary_rows[[mp]] <- summary_row
  decision_rows[[mp]] <- decide_disposition(mp, summary_row, overall_stats, is_mp06 = mp == "MP06")
}

summary_df <- bind_rows(summary_rows)
decision_df <- bind_rows(decision_rows) |>
  left_join(summary_df, by = "metaprogram")

write_csv(summary_df, file.path(params$tables_dir, "10b_audit_summary.csv"))
write_csv(decision_df, file.path(params$tables_dir, "10b_audit_decision.csv"))

sanity <- tibble(
  n_cells = nrow(md),
  n_mp02_high = sum(md$MP02 >= quantile(md$MP02, params$top_quantile, na.rm = TRUE)),
  n_mp06_high = sum(md$MP06 >= quantile(md$MP06, params$top_quantile, na.rm = TRUE)),
  cnv_score_available = params$cnv_score_col %in% colnames(md),
  doublet_pANN_available = params$doublet_score_col %in% colnames(md),
  n_glial_markers_present = length(glial_present),
  n_endothelial_markers_present = length(endo_present),
  n_mural_markers_present = length(mural_present),
  n_myeloid_markers_present = length(myeloid_present),
  profile_files = paste(profile_files, collapse = ";")
)
write_csv(sanity, file.path(params$tables_dir, "10b_audit_sanity_checks.csv"))
write_session_info(file.path(params$tables_dir, "10b_audit_session_info.txt"))

stopifnot(nrow(summary_df) == length(params$mp_cols))
stopifnot(nrow(decision_df) == length(params$mp_cols))
stopifnot(all(sanity$n_glial_markers_present >= 3))
stopifnot(all(sanity$n_endothelial_markers_present >= 4))
stopifnot(all(sanity$n_mural_markers_present >= 5))
stopifnot(all(sanity$n_myeloid_markers_present >= 10))

msg("10b post-metaprogram audit completed")
print(summary_df)
print(decision_df |> select(metaprogram, disposition, reason))
