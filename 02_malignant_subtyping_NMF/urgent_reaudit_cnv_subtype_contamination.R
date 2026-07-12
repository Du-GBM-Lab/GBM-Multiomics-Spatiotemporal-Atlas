suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(dplyr)
  library(readr)
})

# URGENT re-audit:
#   Check whether Subtype3/Subtype4 look CNV-high comparable to Subtype1/2
#   or suspiciously CNV-low, and summarize patient composition and doublet risk.
#
# Inputs:
#   outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2
#   tables/10b_audit_summary.csv
#   tables/10b_audit_decision.csv
#   tables/03c_subtype_k4_patient_composition.csv
#
# Outputs:
#   tables/urgent_reaudit_cnv_by_subtype.csv
#   tables/urgent_reaudit_cnv_relative_to_S1S2.csv
#   tables/urgent_reaudit_patient_subtype_counts.csv
#   tables/urgent_reaudit_patient_subtype_dominance.csv
#   tables/urgent_reaudit_niche_audit_summary_existing.csv
#   tables/urgent_reaudit_threshold_summary.csv
#   tables/urgent_reaudit_session_info.txt

step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir) {
  setwd(file.path(getwd(), step_dir))
}

obj_path <- file.path("outputs", "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2")
out_cnv <- file.path("tables", "urgent_reaudit_cnv_by_subtype.csv")
out_rel <- file.path("tables", "urgent_reaudit_cnv_relative_to_S1S2.csv")
out_patient_counts <- file.path("tables", "urgent_reaudit_patient_subtype_counts.csv")
out_patient_dominance <- file.path("tables", "urgent_reaudit_patient_subtype_dominance.csv")
out_niche <- file.path("tables", "urgent_reaudit_niche_audit_summary_existing.csv")
out_threshold <- file.path("tables", "urgent_reaudit_threshold_summary.csv")
out_session <- file.path("tables", "urgent_reaudit_session_info.txt")

obj <- qs2::qs_read(obj_path)
md <- obj@meta.data

required_cols <- c(
  "subtype_k4", "Pt_number", "doubletfinder_pANN", "doubletfinder_class",
  "infercnv_immune_cnv_burden", "infercnv_immune_cnv_burden_z",
  "infercnv_immune_cnv_correlation", "infercnv_immune_cnv_correlation_ref_z",
  "infercnv_immune_call", "malignant_call_status_immune",
  "is_malignant_for_downstream_immune"
)
missing_cols <- setdiff(required_cols, colnames(md))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

subtype_levels <- paste0("Subtype", 1:4)
md$subtype_k4 <- factor(as.character(md$subtype_k4), levels = subtype_levels)
stopifnot(!any(is.na(md$subtype_k4)))

summarise_num <- function(x) {
  tibble(
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}

cnv_by_subtype <- md %>%
  as_tibble(rownames = "cell_id") %>%
  group_by(subtype_k4) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(Pt_number),
    cnv_burden_mean = mean(infercnv_immune_cnv_burden, na.rm = TRUE),
    cnv_burden_median = median(infercnv_immune_cnv_burden, na.rm = TRUE),
    cnv_burden_q25 = quantile(infercnv_immune_cnv_burden, 0.25, na.rm = TRUE),
    cnv_burden_q75 = quantile(infercnv_immune_cnv_burden, 0.75, na.rm = TRUE),
    cnv_burden_z_mean = mean(infercnv_immune_cnv_burden_z, na.rm = TRUE),
    cnv_burden_z_median = median(infercnv_immune_cnv_burden_z, na.rm = TRUE),
    cnv_burden_z_q25 = quantile(infercnv_immune_cnv_burden_z, 0.25, na.rm = TRUE),
    cnv_burden_z_q75 = quantile(infercnv_immune_cnv_burden_z, 0.75, na.rm = TRUE),
    cnv_correlation_mean = mean(infercnv_immune_cnv_correlation, na.rm = TRUE),
    cnv_correlation_median = median(infercnv_immune_cnv_correlation, na.rm = TRUE),
    cnv_correlation_q25 = quantile(infercnv_immune_cnv_correlation, 0.25, na.rm = TRUE),
    cnv_correlation_q75 = quantile(infercnv_immune_cnv_correlation, 0.75, na.rm = TRUE),
    cnv_correlation_ref_z_mean = mean(infercnv_immune_cnv_correlation_ref_z, na.rm = TRUE),
    cnv_correlation_ref_z_median = median(infercnv_immune_cnv_correlation_ref_z, na.rm = TRUE),
    cnv_correlation_ref_z_q25 = quantile(infercnv_immune_cnv_correlation_ref_z, 0.25, na.rm = TRUE),
    cnv_correlation_ref_z_q75 = quantile(infercnv_immune_cnv_correlation_ref_z, 0.75, na.rm = TRUE),
    pANN_mean = mean(doubletfinder_pANN, na.rm = TRUE),
    pANN_median = median(doubletfinder_pANN, na.rm = TRUE),
    pANN_q75 = quantile(doubletfinder_pANN, 0.75, na.rm = TRUE),
    doublet_rate_pct = mean(doubletfinder_class != "Singlet", na.rm = TRUE) * 100,
    high_conf_malignant_pct = mean(is_malignant_for_downstream_immune == TRUE, na.rm = TRUE) * 100,
    .groups = "drop"
  )

readr::write_csv(cnv_by_subtype, out_cnv)

s1s2 <- cnv_by_subtype %>%
  filter(subtype_k4 %in% c("Subtype1", "Subtype2")) %>%
  summarise(
    s1s2_cnv_burden_z_mean = weighted.mean(cnv_burden_z_mean, n_cells),
    s1s2_cnv_correlation_ref_z_mean = weighted.mean(cnv_correlation_ref_z_mean, n_cells),
    s1s2_cnv_burden_mean = weighted.mean(cnv_burden_mean, n_cells),
    s1s2_cnv_correlation_mean = weighted.mean(cnv_correlation_mean, n_cells)
  )

relative <- cnv_by_subtype %>%
  bind_cols(s1s2[rep(1, nrow(cnv_by_subtype)), ]) %>%
  mutate(
    burden_z_ratio_vs_S1S2 = cnv_burden_z_mean / s1s2_cnv_burden_z_mean,
    corr_ref_z_ratio_vs_S1S2 = cnv_correlation_ref_z_mean / s1s2_cnv_correlation_ref_z_mean,
    burden_ratio_vs_S1S2 = cnv_burden_mean / s1s2_cnv_burden_mean,
    corr_ratio_vs_S1S2 = cnv_correlation_mean / s1s2_cnv_correlation_mean,
    burden_z_drop_pct_vs_S1S2 = (1 - burden_z_ratio_vs_S1S2) * 100,
    corr_ref_z_drop_pct_vs_S1S2 = (1 - corr_ref_z_ratio_vs_S1S2) * 100
  ) %>%
  select(
    subtype_k4, n_cells,
    cnv_burden_z_mean, s1s2_cnv_burden_z_mean, burden_z_ratio_vs_S1S2, burden_z_drop_pct_vs_S1S2,
    cnv_correlation_ref_z_mean, s1s2_cnv_correlation_ref_z_mean, corr_ref_z_ratio_vs_S1S2, corr_ref_z_drop_pct_vs_S1S2,
    cnv_burden_mean, burden_ratio_vs_S1S2,
    cnv_correlation_mean, corr_ratio_vs_S1S2,
    pANN_mean, doublet_rate_pct, high_conf_malignant_pct
  )

readr::write_csv(relative, out_rel)

patient_counts <- md %>%
  as_tibble(rownames = "cell_id") %>%
  count(Pt_number, subtype_k4, name = "n_cells") %>%
  group_by(Pt_number) %>%
  mutate(patient_total = sum(n_cells), pct_within_patient = n_cells / patient_total * 100) %>%
  ungroup()

readr::write_csv(patient_counts, out_patient_counts)

patient_dominance <- patient_counts %>%
  group_by(subtype_k4) %>%
  arrange(desc(n_cells), .by_group = TRUE) %>%
  mutate(rank_within_subtype = row_number()) %>%
  summarise(
    subtype_total = sum(n_cells),
    n_patients_contributing = n_distinct(Pt_number),
    top_patient = Pt_number[which.max(n_cells)],
    top_patient_cells = max(n_cells),
    top_patient_pct_of_subtype = max(n_cells) / subtype_total * 100,
    patients_gt30pct_within_patient = paste(Pt_number[pct_within_patient > 30], collapse = ";"),
    n_patients_gt30pct_within_patient = sum(pct_within_patient > 30),
    .groups = "drop"
  )

readr::write_csv(patient_dominance, out_patient_dominance)

niche_summary_path <- file.path("tables", "10b_audit_summary.csv")
niche_decision_path <- file.path("tables", "10b_audit_decision.csv")
if (file.exists(niche_summary_path)) {
  niche_summary <- readr::read_csv(niche_summary_path, show_col_types = FALSE)
  niche_decision <- if (file.exists(niche_decision_path)) {
    readr::read_csv(niche_decision_path, show_col_types = FALSE)
  } else {
    tibble()
  }
  readr::write_csv(
    list(summary = niche_summary, decision = niche_decision) %>%
      bind_rows(.id = "source_table"),
    out_niche
  )
}

threshold_summary <- tibble::tibble(
  malignant_definition = "immune-reference two-axis inferCNV high-confidence",
  cnv_burden_z_threshold = "> 3",
  cnv_correlation_ref_z_threshold = "> 3",
  included_downstream_call = "malignant_like_CNV_high_confidence / is_malignant_for_downstream_immune == TRUE",
  excluded_from_main_subtyping = "malignant_like_CNV_burden_only",
  final_object_cells = nrow(md),
  high_conf_malignant_pct_in_final_object = mean(md$is_malignant_for_downstream_immune == TRUE, na.rm = TRUE) * 100,
  note = "All cells in final malignant subtype object should be high-confidence malignant by immune-reference inferCNV two-axis call."
)
readr::write_csv(threshold_summary, out_threshold)

writeLines(capture.output(sessionInfo()), out_session)

cat("CNV by subtype:\n")
print(cnv_by_subtype)
cat("\nRelative to weighted S1/S2 baseline:\n")
print(relative)
cat("\nPatient dominance:\n")
print(patient_dominance)
cat("\nThreshold summary:\n")
print(threshold_summary)
cat("\nOutputs:\n")
cat(out_cnv, "\n", out_rel, "\n", out_patient_counts, "\n", out_patient_dominance, "\n", out_niche, "\n", out_threshold, "\n", sep = "")
