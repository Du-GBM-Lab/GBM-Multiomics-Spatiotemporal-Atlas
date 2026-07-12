#!/usr/bin/env Rscript
# =============================================================================
# R9 batch 1 | A2 D02 vascular-fail slice stratification
#
# Sensitivity analysis only. D02 definition is frozen from script 44 output.
# No reclustering, no resolution/k tuning, no D02 boundary edits.
#
# Goal:
#   Explain the 6 slices where D02 vascular enrichment is not FDR-significant
#   using pre-specified slice-level variables:
#     1) IDH_status
#     2) D02_vascular_spot_n
#     3) slice_vascular_spot_total
#
# Vascular spot counts are slice-level dominant-RCTD counts:
#   dominant RCTD label = max weight per spot; vascular = Endothelial or Mural cells.
# This is a count-based slice metric, not a vascular-score spot subset filter.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
})

base_dir <- getwd()
gate_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A2_D02_gate_review")
out_dir <- file.path(gate_dir, "slice_stratification")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

enrich_path <- file.path(gate_dir, "A2_D02_perslice_enrichment_vs_nonD02.csv")
wide_path <- file.path(gate_dir, "A2_D02_perslice_gate_wide.csv")
prox_summary_path <- file.path(gate_dir, "A2_D02_region_v3B_proximity_summary.csv")
null_audit_path <- file.path(gate_dir, "A2_D02_region_v3B_null_audit.csv")
weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
labels_path <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A2_unbiased_domain/A2_domain_labels_per_spot.csv")

enrich <- fread(enrich_path)
wide <- fread(wide_path)
prox_summary <- fread(prox_summary_path)
null_audit <- fread(null_audit_path)
weights <- as.data.table(qs2::qs_read(weights_path))
labels <- fread(labels_path)

vascular_fail <- enrich[feature == "vascular" & !pass_positive_fdr_feature]
vascular_negative <- enrich[feature == "vascular" & log2FC <= 0]
vascular_positive_nonsig <- enrich[feature == "vascular" & log2FC > 0 & !pass_positive_fdr_feature]

fail_list <- vascular_fail[, .(
  slice,
  failure_class = fifelse(log2FC > 0, "positive_but_FDR_ge_0.05", "negative_or_zero"),
  n_region,
  slice_spots,
  region_fraction,
  vascular_mean_D02 = mean_region,
  vascular_mean_background = mean_background_nonD02,
  vascular_log2FC = log2FC,
  vascular_fdr = fdr_bh_by_feature
)]
fail_list[, frozen_order := .I]

cell_cols <- setdiff(names(weights), c("spot_id", "slice", "image", "x", "y"))
if (!all(c("Endothelial", "Mural cells") %in% cell_cols)) {
  stop("RCTD weights do not contain Endothelial and Mural cells columns.")
}
weights[, dominant_RCTD_label := cell_cols[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = cell_cols]
weights[, dominant_is_vascular := dominant_RCTD_label %in% c("Endothelial", "Mural cells")]

spot_meta <- merge(
  labels[, .(spot_id, slice, domain_id, D02_region = domain_id == "D02")],
  weights[, .(spot_id, slice, dominant_RCTD_label, dominant_is_vascular,
              vascular_weight = Endothelial + `Mural cells`)],
  by = c("spot_id", "slice"),
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(spot_meta$dominant_is_vascular)) stop("dominant vascular flag missing after label-weight merge.")

slice_meta <- spot_meta[, .(
  slice_spots = .N,
  IDH_status = ifelse(grepl("IDHMutant", slice), "IDH-Mutant", "IDH-WT"),
  slice_vascular_spot_total = sum(dominant_is_vascular),
  slice_vascular_spot_fraction = mean(dominant_is_vascular),
  D02_size = sum(D02_region),
  D02_vascular_spot_n = sum(D02_region & dominant_is_vascular),
  D02_vascular_spot_fraction = fifelse(sum(D02_region) > 0, sum(D02_region & dominant_is_vascular) / sum(D02_region), NA_real_),
  slice_vascular_weight_mean = mean(vascular_weight),
  D02_vascular_weight_mean = mean(vascular_weight[D02_region])
), by = slice]

slice_meta <- merge(
  slice_meta,
  wide[, .(
    slice,
    D02_vascular_log2FC = `log2FC_vascular`,
    D02_vascular_FDR = `fdr_bh_by_feature_vascular`,
    D02_vascular_positive_FDR = `pass_positive_fdr_feature_vascular`,
    D02_MESV_positive_FDR = `pass_positive_fdr_feature_MES-V`,
    D02_MESV_vascular_positive_FDR = pass_MESV_vascular
  )],
  by = "slice",
  all.x = TRUE,
  sort = FALSE
)
slice_meta[, vascular_fail := !D02_vascular_positive_FDR]
slice_meta[, vascular_fail_class := fifelse(
  !vascular_fail, "pass",
  fifelse(D02_vascular_log2FC > 0, "positive_but_FDR_ge_0.05", "negative_or_zero")
)]

fisher_or_na <- function(tab) {
  if (any(dim(tab) < 2)) return(NA_real_)
  suppressWarnings(fisher.test(tab)$p.value)
}
wilcox_or_na <- function(x, group) {
  if (length(unique(group)) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(x ~ group, exact = FALSE)$p.value)
}

strat_tests <- list()
strat_tests[[1]] <- data.table(
  variable = "IDH_status",
  priority = 1L,
  test = "Fisher exact: vascular_fail by IDH_status",
  p_value = fisher_or_na(table(slice_meta$vascular_fail, slice_meta$IDH_status)),
  fail_distribution = paste(capture.output(print(table(slice_meta[vascular_fail == TRUE]$IDH_status))), collapse = " | "),
  pass_distribution = paste(capture.output(print(table(slice_meta[vascular_fail == FALSE]$IDH_status))), collapse = " | ")
)
for (var in c("D02_vascular_spot_n", "slice_vascular_spot_total")) {
  strat_tests[[length(strat_tests) + 1L]] <- data.table(
    variable = var,
    priority = ifelse(var == "D02_vascular_spot_n", 2L, 3L),
    test = "Wilcoxon: vascular_fail vs pass",
    p_value = wilcox_or_na(slice_meta[[var]], slice_meta$vascular_fail),
    fail_distribution = paste0(
      "median=", median(slice_meta[vascular_fail == TRUE][[var]]),
      "; range=", min(slice_meta[vascular_fail == TRUE][[var]]), "-", max(slice_meta[vascular_fail == TRUE][[var]])
    ),
    pass_distribution = paste0(
      "median=", median(slice_meta[vascular_fail == FALSE][[var]]),
      "; range=", min(slice_meta[vascular_fail == FALSE][[var]]), "-", max(slice_meta[vascular_fail == FALSE][[var]])
    )
  )
}
strat_tests <- rbindlist(strat_tests, use.names = TRUE)

# Pre-specified, non-outcome-fitted power thresholds for exploratory reporting.
# D02_vascular_spot_n threshold uses a minimal count of 10 vascular-dominant spots
# inside D02 as a conservative small-count warning.
# Whole-slice vascular threshold uses the lower quartile of slice vascular dominant
# counts as a descriptive low-vascular stratum.
thresholds <- data.table(
  threshold_name = c("IDH_WT_only", "D02_vascular_spot_n_ge_10", "slice_vascular_spot_total_ge_Q1"),
  priority = c(1L, 2L, 3L),
  variable = c("IDH_status", "D02_vascular_spot_n", "slice_vascular_spot_total"),
  threshold = c("IDH-WT", ">=10", paste0(">=", as.integer(quantile(slice_meta$slice_vascular_spot_total, 0.25, type = 1)))),
  rationale = c(
    "pre-specified biological subgroup; IDH-mutant slices may have different vascular biology",
    "pre-specified minimum count warning for region-level vascular test power; not fitted to outcome",
    "pre-specified descriptive low-whole-slice-vascular stratum using Q1 of the independent slice-level count"
  )
)

eval_subset <- function(name, keep, reason) {
  sub <- slice_meta[keep]
  data.table(
    threshold_name = name,
    N = nrow(sub),
    n_pass_vascular_FDR = sum(sub$D02_vascular_positive_FDR),
    n_positive = sum(sub$D02_vascular_log2FC > 0),
    pass_rate = sum(sub$D02_vascular_positive_FDR) / nrow(sub),
    positive_rate = sum(sub$D02_vascular_log2FC > 0) / nrow(sub),
    excluded_slices = paste(slice_meta[!keep]$slice, collapse = ";"),
    exclusion_reason = reason
  )
}
q1_vascular_total <- as.integer(quantile(slice_meta$slice_vascular_spot_total, 0.25, type = 1))
subset_results <- rbindlist(list(
  eval_subset("IDH_WT_only", slice_meta$IDH_status == "IDH-WT", "IDH-Mutant"),
  eval_subset("D02_vascular_spot_n_ge_10", slice_meta$D02_vascular_spot_n >= 10, "D02_vascular_spot_n < 10"),
  eval_subset("slice_vascular_spot_total_ge_Q1", slice_meta$slice_vascular_spot_total >= q1_vascular_total, "slice_vascular_spot_total below Q1")
), use.names = TRUE)

sens_thresholds <- rbindlist(list(
  data.table(variable = "D02_vascular_spot_n", threshold = c(5L, 10L, 15L, 20L)),
  data.table(variable = "slice_vascular_spot_total", threshold = as.integer(quantile(slice_meta$slice_vascular_spot_total, c(0.10, 0.25, 0.33, 0.50), type = 1)))
), use.names = TRUE)
sens_results <- sens_thresholds[, {
  keep <- slice_meta[[variable]] >= threshold
  .(
    N = sum(keep),
    n_excluded = sum(!keep),
    n_pass_vascular_FDR = sum(slice_meta[keep]$D02_vascular_positive_FDR),
    pass_rate = sum(slice_meta[keep]$D02_vascular_positive_FDR) / sum(keep),
    excluded_slices = paste(slice_meta[!keep]$slice, collapse = ";")
  )
}, by = .(variable, threshold)]

excluded_reasons <- rbindlist(list(
  slice_meta[IDH_status == "IDH-Mutant", .(slice, reason = "IDH-Mutant", variable = "IDH_status")],
  slice_meta[D02_vascular_spot_n < 10, .(slice, reason = "D02_vascular_spot_n < 10", variable = "D02_vascular_spot_n")],
  slice_meta[slice_vascular_spot_total < q1_vascular_total, .(slice, reason = paste0("slice_vascular_spot_total < Q1(", q1_vascular_total, ")"), variable = "slice_vascular_spot_total")]
), use.names = TRUE)

null_check <- null_audit[, .(
  n_rows = .N,
  n_preserved = sum(count_preserved_all_perm),
  all_preserved = all(count_preserved_all_perm)
)]

neuron_check <- prox_summary[
  sensitivity == "all_slices" & top_rule == "top10" &
    ((response == "D02_region" & predictor == "neuron_control") |
       (response == "neuron_control" & predictor == "D02_region")),
  .(k, response, predictor, median_delta_crossK, n_crossK_positive, median_p_crossK,
    n_slice_crossK_p_lt_0p05, median_delta_nn_closer, n_nn_closer_positive, median_p_nn)
]

outcome <- "outcome_2_retreat_to_supplement_plus_main_text_one_sentence"
outcome_reason <- paste(
  "The six vascular-FDR-fail slices do not concentrate in the pre-specified IDH-Mutant subgroup.",
  "The two IDH-Mutant slices split one fail and one pass.",
  "Small D02 vascular-dominant count explains only part of the failures and also excludes several passing slices.",
  "No pre-specified stratum yields a clean 15-17/18-equivalent promotion path."
)

decision <- data.table(
  outcome = outcome,
  outcome_label = "退补充 + 正文一句引用",
  candidate_for_main_text_overview_via_stratified_sensitivity = FALSE,
  reason = outcome_reason,
  allowed_wording = "D02 may be cited as an unbiased supplementary landscape/domain signal consistent with the locked MES-like vascular proximity result.",
  forbidden_wording = "No causal proof, physical contact, single-cell colocalization, interaction, recruitment, vascular drives MES, TAM recruitment, or myeloid/TAM reinterpretation."
)

params <- data.table(
  parameter = c(
    "D02_definition",
    "vascular_fail_definition",
    "IDH_status_source",
    "vascular_spot_count_source",
    "D02_vascular_spot_n_threshold",
    "slice_vascular_spot_total_threshold",
    "QC_metadata",
    "priority_order",
    "boundary"
  ),
  value = c(
    "Frozen A2 D02 labels from A2_domain_labels_per_spot.csv and script 44 output; no reclustering or boundary edits.",
    "D02 vascular enrichment pass_positive_fdr_feature == FALSE in A2_D02_perslice_enrichment_vs_nonD02.csv.",
    "Parsed from slice ID; slices containing IDHMutant are IDH-Mutant, all others IDH-WT.",
    "Dominant-RCTD label per spot from R9_A2_RCTD_weights_allslices_long.qs2; vascular = Endothelial or Mural cells.",
    ">=10 vascular-dominant spots inside D02, pre-specified small-count warning.",
    paste0(">= Q1 of slice_vascular_spot_total, Q1=", q1_vascular_total, "."),
    "No per-slice UMI/gene QC source table was used in this stratification.",
    "1 IDH_status; 2 D02_vascular_spot_n; 3 slice_vascular_spot_total.",
    "Exploratory slice-level sensitivity only; no spot subset filtering by vascular score; no main figure promoted."
  )
)

fwrite(fail_list, file.path(out_dir, "A2_D02_vascular_fail_frozen_slices.csv"))
fwrite(slice_meta, file.path(out_dir, "A2_D02_slice_level_metadata_and_gate.csv"))
fwrite(strat_tests, file.path(out_dir, "A2_D02_prespecified_stratification_tests.csv"))
fwrite(thresholds, file.path(out_dir, "A2_D02_stratification_thresholds.csv"))
fwrite(subset_results, file.path(out_dir, "A2_D02_stratified_pass_rates.csv"))
fwrite(sens_results, file.path(out_dir, "A2_D02_threshold_sensitivity.csv"))
fwrite(excluded_reasons, file.path(out_dir, "A2_D02_excluded_slice_reasons_by_threshold.csv"))
fwrite(neuron_check, file.path(out_dir, "A2_D02_neuron_control_stratification_check.csv"))
fwrite(null_check, file.path(out_dir, "A2_D02_null_count_preservation_check.csv"))
fwrite(decision, file.path(out_dir, "A2_D02_stratification_decision.csv"))
fwrite(params, file.path(out_dir, "A2_D02_stratification_parameters.csv"))

cat("\n== Step 1 frozen vascular-fail slices:\n")
print(fail_list)
cat("\n== Pre-specified slice-level metadata:\n")
print(slice_meta[order(vascular_fail, IDH_status, slice)])
cat("\n== Stratification tests:\n")
print(strat_tests)
cat("\n== Stratified pass rates:\n")
print(subset_results)
cat("\n== Threshold sensitivity:\n")
print(sens_results)
cat("\n== Decision:\n")
print(decision)
cat("\n[STOP A2 D02 vascular-fail slice stratification] Review outcome before any figure-role change.\n")
