#!/usr/bin/env Rscript
# =============================================================================
# R9 Batch1 | A1 toroidal fail-slice identity check
#
# Nature:
#   Interpretability check only. The A1 program toroidal test from script 47 is
#   frozen: no rerun, no statistic/null/neighborhood change, no slice picking.
#   This script only reads existing per-slice results and already-recorded batch
#   1 weak-slice metadata from A2/D02 stratification. No new weak-slice variable
#   or threshold is introduced.
#
# Binary decision:
#   If A1 failing slices are contained in pre-existing independent weak-slice
#   records, mark candidate for stratified sensitivity discussion. Otherwise,
#   keep A1 as strong supplement/source-data with at most one main-text citation.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- getwd()
a1_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A1_program_colocalization_toroidal")
a2_strat_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A2_D02_gate_review/slice_stratification")
out_dir <- file.path(a1_dir, "slice_check")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

a1_perslice_path <- file.path(a1_dir, "A1_program_toroidal_perslice.csv")
a1_gate_path <- file.path(a1_dir, "A1_program_toroidal_gate.csv")
a2_fail_path <- file.path(a2_strat_dir, "A2_D02_vascular_fail_frozen_slices.csv")
a2_meta_path <- file.path(a2_strat_dir, "A2_D02_slice_level_metadata_and_gate.csv")
a2_threshold_path <- file.path(a2_strat_dir, "A2_D02_stratification_thresholds.csv")

required_files <- c(a1_perslice_path, a1_gate_path, a2_fail_path, a2_meta_path, a2_threshold_path)
missing <- required_files[!file.exists(required_files)]
if (length(missing)) stop("Missing required file(s):\n", paste(missing, collapse = "\n"))

params <- data.table(
  parameter = c(
    "scope",
    "frozen_A1_input",
    "independent_weak_slice_sources",
    "no_new_strata",
    "decision_rule"
  ),
  value = c(
    "read-only interpretability check; no A1 toroidal rerun or parameter change",
    a1_perslice_path,
    paste(c(a2_fail_path, a2_meta_path, a2_threshold_path), collapse = "; "),
    "TRUE",
    "candidate only if A1 non-pass slices are contained in pre-existing independent weak-slice records; otherwise strong supplement"
  )
)
fwrite(params, file.path(out_dir, "A1_toroidal_fail_slice_check_parameters.csv"))

a1 <- fread(a1_perslice_path)
a1_gate <- fread(a1_gate_path)
a2_fail <- fread(a2_fail_path)
a2_meta <- fread(a2_meta_path)
a2_thresholds <- fread(a2_threshold_path)

primary <- a1[pair == "MESV_vascular"]
if (nrow(primary) != 18) stop("Expected 18 MESV_vascular rows in A1 per-slice output.")

negative_or_zero <- primary[positive == FALSE, .(
  slice,
  A1_failure_class = "A1_direction_negative_or_zero",
  obs_stat,
  null_median,
  delta_vs_null_median,
  emp_p_greater,
  fdr_bh_by_pair,
  positive,
  pass_positive_fdr
)]

positive_nonsig <- primary[positive == TRUE & pass_positive_fdr == FALSE, .(
  slice,
  A1_failure_class = "A1_positive_but_FDR_ge_0.05",
  obs_stat,
  null_median,
  delta_vs_null_median,
  emp_p_greater,
  fdr_bh_by_pair,
  positive,
  pass_positive_fdr
)]

all_nonpass <- primary[pass_positive_fdr == FALSE, .(
  slice,
  A1_failure_class = fifelse(positive, "A1_positive_but_FDR_ge_0.05", "A1_direction_negative_or_zero"),
  obs_stat,
  null_median,
  delta_vs_null_median,
  emp_p_greater,
  fdr_bh_by_pair,
  positive,
  pass_positive_fdr
)]
all_nonpass[, frozen_order := .I]

frozen_lists <- rbindlist(list(negative_or_zero, positive_nonsig), use.names = TRUE)
fwrite(frozen_lists, file.path(out_dir, "A1_toroidal_frozen_failure_slices.csv"))
fwrite(all_nonpass, file.path(out_dir, "A1_toroidal_all_nonpass_slices.csv"))

a2_neg <- a2_fail[failure_class == "negative_or_zero", slice]
a2_fail_slices <- a2_fail$slice

threshold_values <- a2_thresholds[, .(threshold_name, variable, threshold)]
d02_min_count <- as.numeric(gsub("[^0-9.]+", "", threshold_values[threshold_name == "D02_vascular_spot_n_ge_10", threshold]))
whole_vasc_q1 <- as.numeric(gsub("[^0-9.]+", "", threshold_values[threshold_name == "slice_vascular_spot_total_ge_Q1", threshold]))
if (!is.finite(d02_min_count)) d02_min_count <- 10
if (!is.finite(whole_vasc_q1)) whole_vasc_q1 <- 5

weak_meta <- copy(a2_meta)
weak_meta[, weak_D02_vascular_fail_6 := slice %in% a2_fail_slices]
weak_meta[, weak_D02_vascular_negative_or_zero := slice %in% a2_neg]
weak_meta[, weak_D02_vascular_spot_n_below_existing_min := D02_vascular_spot_n < d02_min_count]
weak_meta[, weak_slice_vascular_total_below_existing_Q1 := slice_vascular_spot_total < whole_vasc_q1]
weak_meta[, any_existing_weak_record := weak_D02_vascular_fail_6 |
            weak_D02_vascular_negative_or_zero |
            weak_D02_vascular_spot_n_below_existing_min |
            weak_slice_vascular_total_below_existing_Q1]
weak_meta[, existing_weak_reasons := paste(c(
  if (weak_D02_vascular_fail_6) "A2_D02_vascular_FDR_fail_6" else NULL,
  if (weak_D02_vascular_negative_or_zero) "A2_D02_vascular_negative_or_zero" else NULL,
  if (weak_D02_vascular_spot_n_below_existing_min) paste0("D02_vascular_spot_n<", d02_min_count) else NULL,
  if (weak_slice_vascular_total_below_existing_Q1) paste0("slice_vascular_spot_total<", whole_vasc_q1) else NULL
), collapse = ";"), by = slice]

fwrite(weak_meta, file.path(out_dir, "A1_toroidal_existing_weak_slice_metadata.csv"))

overlap_all <- merge(
  all_nonpass,
  weak_meta[, .(
    slice, IDH_status, slice_vascular_spot_total, D02_vascular_spot_n,
    weak_D02_vascular_fail_6,
    weak_D02_vascular_negative_or_zero,
    weak_D02_vascular_spot_n_below_existing_min,
    weak_slice_vascular_total_below_existing_Q1,
    any_existing_weak_record,
    existing_weak_reasons
  )],
  by = "slice",
  all.x = TRUE,
  sort = FALSE
)
overlap_all <- overlap_all[match(all_nonpass$slice, slice)]
fwrite(overlap_all, file.path(out_dir, "A1_toroidal_failure_overlap_existing_weak_records.csv"))

overlap_summary <- rbindlist(list(
  data.table(
    failure_set = "direction_negative_or_zero_only",
    n_slices = nrow(negative_or_zero),
    slices = paste(negative_or_zero$slice, collapse = ";"),
    n_in_A2_D02_vascular_fail_6 = sum(negative_or_zero$slice %in% a2_fail_slices),
    n_in_A2_D02_negative_or_zero = sum(negative_or_zero$slice %in% a2_neg),
    n_with_any_existing_weak_record = sum(negative_or_zero$slice %in% weak_meta[any_existing_weak_record == TRUE, slice])
  ),
  data.table(
    failure_set = "all_A1_positive_FDR_nonpass",
    n_slices = nrow(all_nonpass),
    slices = paste(all_nonpass$slice, collapse = ";"),
    n_in_A2_D02_vascular_fail_6 = sum(all_nonpass$slice %in% a2_fail_slices),
    n_in_A2_D02_negative_or_zero = sum(all_nonpass$slice %in% a2_neg),
    n_with_any_existing_weak_record = sum(all_nonpass$slice %in% weak_meta[any_existing_weak_record == TRUE, slice])
  )
), use.names = TRUE)
fwrite(overlap_summary, file.path(out_dir, "A1_toroidal_failure_overlap_summary.csv"))

strata <- list(
  exclude_A2_D02_vascular_fail_6 = setdiff(unique(primary$slice), a2_fail_slices),
  keep_D02_vascular_spot_n_ge_10 = weak_meta[D02_vascular_spot_n >= d02_min_count, slice],
  keep_slice_vascular_spot_total_ge_Q1 = weak_meta[slice_vascular_spot_total >= whole_vasc_q1, slice],
  exclude_any_existing_weak_record = weak_meta[any_existing_weak_record == FALSE, slice]
)

stratified_rates <- rbindlist(lapply(names(strata), function(nm) {
  keep <- strata[[nm]]
  pd <- primary[slice %in% keep]
  data.table(
    stratum = nm,
    N = nrow(pd),
    excluded_N = 18L - nrow(pd),
    excluded_slices = paste(setdiff(primary$slice, keep), collapse = ";"),
    n_positive = sum(pd$positive),
    n_pass_positive_fdr = sum(pd$pass_positive_fdr),
    pass_rate = ifelse(nrow(pd) > 0, sum(pd$pass_positive_fdr) / nrow(pd), NA_real_)
  )
}), use.names = TRUE)
fwrite(stratified_rates, file.path(out_dir, "A1_toroidal_existing_strata_pass_rates.csv"))

all_nonpass_explained <- all(overlap_all$any_existing_weak_record)
direction_fail_explained <- all(negative_or_zero$slice %in% weak_meta[any_existing_weak_record == TRUE, slice])
candidate_main <- all_nonpass_explained

decision <- data.table(
  item = c(
    "direction_negative_or_zero_explained_by_existing_weak_records",
    "all_A1_FDR_nonpass_explained_by_existing_weak_records",
    "existing_strata_pass_rates_reported_for_judgment",
    "binary_decision"
  ),
  value = c(
    paste0(sum(negative_or_zero$slice %in% weak_meta[any_existing_weak_record == TRUE, slice]), "/", nrow(negative_or_zero)),
    paste0(sum(overlap_all$any_existing_weak_record), "/", nrow(overlap_all)),
    paste(stratified_rates[, paste0(stratum, "=", n_pass_positive_fdr, "/", N)], collapse = "; "),
    ifelse(candidate_main, "candidate_for_main_text_via_stratified_sensitivity", "strong_supplement_main_text_one_sentence")
  )
)
fwrite(decision, file.path(out_dir, "A1_toroidal_fail_slice_check_decision.csv"))

stop_lines <- c(
  "# R9 Batch1 A1 toroidal fail-slice check STOP",
  "",
  paste0("- Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "- Nature: interpretability check only. A1 toroidal definition is frozen; no rerun, no statistic/null/neighborhood change, no score thresholding, no new weak-slice variable.",
  "- Sources: script 47 A1 per-slice output plus script 45 D02 stratification outputs only.",
  "",
  "## Frozen A1 Failure Lists",
  paste0("- Direction negative/zero slices: ", paste(negative_or_zero$slice, collapse = "; "), "."),
  paste0("- Positive but BH FDR>=0.05 slices: ", paste(positive_nonsig$slice, collapse = "; "), "."),
  paste0("- All A1 FDR non-pass slices: ", paste(all_nonpass$slice, collapse = "; "), "."),
  "",
  "## Overlap With Existing Weak-Slice Records",
  paste0("- Direction negative/zero overlap with any existing weak record: ",
         sum(negative_or_zero$slice %in% weak_meta[any_existing_weak_record == TRUE, slice]), "/", nrow(negative_or_zero), "."),
  paste0("- All A1 FDR non-pass overlap with any existing weak record: ",
         sum(overlap_all$any_existing_weak_record), "/", nrow(overlap_all), "."),
  paste0("- Overlap table: ", normalizePath(file.path(out_dir, "A1_toroidal_failure_overlap_existing_weak_records.csv"), winslash = "/")),
  "",
  "## Existing Strata Pass Rates",
  paste(apply(stratified_rates[, .(stratum, N, n_pass_positive_fdr, excluded_N)], 1, function(x) {
    paste0("- ", x[["stratum"]], ": ", x[["n_pass_positive_fdr"]], "/", x[["N"]], " pass; excluded ", x[["excluded_N"]], " slices.")
  }), collapse = "\n"),
  "",
  "## Decision",
  paste0("- Binary decision: ", decision[item == "binary_decision", value], "."),
  if (candidate_main) {
    "- Candidate caveat: all A1 non-pass slices are covered by existing independent weak-slice records; use only the reported pre-existing strata for judgment, and do not treat this as changing the frozen 14/18 primary result."
  } else {
    "- Placement: A1 remains strong supplement/source-data; main text at most one sentence that A1 program-level direction is consistent with v3-B."
  },
  "- Forbidden: no physical/single-cell colocalization, interaction, recruitment, drives, causality, myeloid/TAM reinterpretation."
)
writeLines(stop_lines, file.path(base_dir, "docs/R9_batch1_A1_toroidal_fail_slice_check_STOP.md"), useBytes = TRUE)

message("Done. Outputs written to: ", out_dir)
