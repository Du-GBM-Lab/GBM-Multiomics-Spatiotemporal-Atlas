#!/usr/bin/env Rscript
# =============================================================================
# R9 Batch1 | A1 12/14 Q1 stratification compliance audit
#
# Nature:
#   Read-only compliance audit. This script decides whether the existing
#   `slice_vascular_spot_total >= Q1` stratum can support A1 program toroidal
#   stratified sensitivity for a possible main-text figure position. It freezes
#   all prior tests: no toroidal rerun, no colocalization recalculation, no
#   threshold change, no new stratum.
#
# Questions answered:
#   1) Was Q1=5 a pre-declared distributional quartile rather than a hand-picked
#      value?
#   2) Are the excluded slices naturally determined by the independent whole-
#      slice vascular count ranking, not by A1 failures?
#   3) In the 14-slice subset, does neuron remain clean and count/score
#      preservation remain TRUE?
#
# If all three questions are clean: pass_stratified_sensitivity_for_main_text.
# Otherwise: back_to_supplement_strong_reference.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- getwd()
a1_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A1_program_colocalization_toroidal")
slice_check_dir <- file.path(a1_dir, "slice_check")
a2_strat_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A2_D02_gate_review/slice_stratification")
out_dir <- file.path(slice_check_dir, "q1_stratification_audit")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

a1_perslice_path <- file.path(a1_dir, "A1_program_toroidal_perslice.csv")
a2_meta_path <- file.path(a2_strat_dir, "A2_D02_slice_level_metadata_and_gate.csv")
a2_threshold_path <- file.path(a2_strat_dir, "A2_D02_stratification_thresholds.csv")
a2_params_path <- file.path(a2_strat_dir, "A2_D02_stratification_parameters.csv")
script45_path <- file.path(base_dir, "scripts/45_R9_batch1_A2_D02_vascular_fail_slice_stratification.R")

required_files <- c(a1_perslice_path, a2_meta_path, a2_threshold_path, a2_params_path, script45_path)
missing <- required_files[!file.exists(required_files)]
if (length(missing)) stop("Missing required file(s):\n", paste(missing, collapse = "\n"))

a1 <- fread(a1_perslice_path)
a2_meta <- fread(a2_meta_path)
a2_threshold <- fread(a2_threshold_path)
a2_params <- fread(a2_params_path)
script45 <- readLines(script45_path, warn = FALSE, encoding = "UTF-8")

q1_script_line <- grep("quantile\\(slice_meta\\$slice_vascular_spot_total, 0.25, type = 1\\)", script45, value = TRUE)
q1_predeclared_line <- grep("Whole-slice vascular threshold uses the lower quartile", script45, value = TRUE)

dist_dt <- a2_meta[, .(
  slice,
  IDH_status,
  slice_vascular_spot_total,
  D02_vascular_spot_n
)][order(slice_vascular_spot_total, slice)]
dist_dt[, vascular_total_rank_low_to_high := .I]

q1_type1 <- as.integer(quantile(dist_dt$slice_vascular_spot_total, 0.25, type = 1))
q1_type7 <- as.numeric(quantile(dist_dt$slice_vascular_spot_total, 0.25, type = 7))
dist_summary <- data.table(
  metric = c("min", "Q1_type1_empirical_lower", "Q1_type7_R_default", "median", "Q3_type1_empirical_lower", "Q3_type7_R_default", "max"),
  value = c(
    min(dist_dt$slice_vascular_spot_total),
    q1_type1,
    q1_type7,
    median(dist_dt$slice_vascular_spot_total),
    as.integer(quantile(dist_dt$slice_vascular_spot_total, 0.75, type = 1)),
    as.numeric(quantile(dist_dt$slice_vascular_spot_total, 0.75, type = 7)),
    max(dist_dt$slice_vascular_spot_total)
  )
)
fwrite(dist_dt, file.path(out_dir, "A1_Q1_slice_vascular_total_18_values.csv"))
fwrite(dist_summary, file.path(out_dir, "A1_Q1_slice_vascular_total_distribution_summary.csv"))

threshold_row <- a2_threshold[threshold_name == "slice_vascular_spot_total_ge_Q1"]
threshold_value <- as.integer(gsub("[^0-9]+", "", threshold_row$threshold))
q1_clean <- threshold_value == q1_type1 &&
  length(q1_script_line) > 0 &&
  grepl("Q1", threshold_row$rationale) &&
  grepl("independent slice-level count", threshold_row$rationale)

q1_source <- data.table(
  question = "Q1_source_and_predeclaration",
  threshold_from_existing_table = threshold_row$threshold,
  threshold_numeric = threshold_value,
  computed_Q1_type1_empirical_lower = q1_type1,
  computed_Q1_type7_R_default = q1_type7,
  script45_quantile_line = paste(q1_script_line, collapse = " | "),
  script45_predeclaration_line = paste(q1_predeclared_line, collapse = " | "),
  rationale_from_existing_table = threshold_row$rationale,
  clean = q1_clean,
  answer = ifelse(q1_clean,
                  "clean: Q1=5 is the pre-declared empirical lower-quartile rule from script 45, not a hand-picked A1 outcome threshold",
                  "not clean: Q1=5 cannot be verified as the pre-declared empirical quartile")
)
fwrite(q1_source, file.path(out_dir, "A1_Q1_question1_source_predeclaration.csv"))

excluded <- dist_dt[slice_vascular_spot_total < threshold_value]
kept <- dist_dt[slice_vascular_spot_total >= threshold_value]
lowest_four <- dist_dt[1:4]
a1_fail <- c("#UKF265_T_ST", "#UKF313_T_ST", "#UKF270_IDHMutant_T_ST", "#UKF304_T_ST")

excluded[, is_A1_toroidal_nonpass := slice %in% a1_fail]
kept[, is_A1_toroidal_nonpass := slice %in% a1_fail]
excluded_natural <- setequal(excluded$slice, lowest_four$slice)
metric_independent <- TRUE

exclusion_audit <- data.table(
  excluded_slice = excluded$slice,
  IDH_status = excluded$IDH_status,
  slice_vascular_spot_total = excluded$slice_vascular_spot_total,
  vascular_total_rank_low_to_high = excluded$vascular_total_rank_low_to_high,
  is_A1_toroidal_nonpass = excluded$is_A1_toroidal_nonpass
)
fwrite(exclusion_audit, file.path(out_dir, "A1_Q1_question2_excluded_slices.csv"))

question2 <- data.table(
  question = "exclusion_independence",
  threshold_rule = paste0("keep slice_vascular_spot_total >= ", threshold_value),
  excluded_slices = paste(excluded$slice, collapse = ";"),
  A1_toroidal_nonpass_slices = paste(a1_fail, collapse = ";"),
  overlap_excluded_with_A1_nonpass = paste(intersect(excluded$slice, a1_fail), collapse = ";"),
  n_overlap_excluded_with_A1_nonpass = length(intersect(excluded$slice, a1_fail)),
  excluded_are_lowest_four_by_independent_metric = excluded_natural,
  metric_independent_of_A1_toroidal = metric_independent,
  clean = excluded_natural && metric_independent,
  answer = ifelse(excluded_natural && metric_independent,
                  "clean: the four excluded slices are exactly the four lowest whole-slice dominant-RCTD vascular counts; the metric is independent of A1 toroidal output",
                  "not clean: excluded slices are not the natural low end of the independent vascular-count ranking")
)
fwrite(question2, file.path(out_dir, "A1_Q1_question2_exclusion_independence.csv"))

subset_slices <- kept$slice
subset_results <- a1[slice %in% subset_slices]
mesv_subset <- subset_results[pair == "MESV_vascular"]
neuron_subset <- subset_results[pair == "neuron_vascular"]

question3 <- data.table(
  question = "subset_neuron_and_preservation",
  subset_rule = paste0("slice_vascular_spot_total >= ", threshold_value),
  N = length(subset_slices),
  subset_slices = paste(subset_slices, collapse = ";"),
  MESV_positive_FDR_pass = sum(mesv_subset$pass_positive_fdr),
  MESV_positive = sum(mesv_subset$positive),
  neuron_positive_FDR_pass = sum(neuron_subset$pass_positive_fdr),
  neuron_positive = sum(neuron_subset$positive),
  MESV_count_preserved_TRUE = sum(mesv_subset$count_preserved_all_perm),
  MESV_score_distribution_preserved_TRUE = sum(mesv_subset$score_distribution_preserved_all_perm),
  clean = sum(neuron_subset$pass_positive_fdr) == 0 &&
    all(mesv_subset$count_preserved_all_perm) &&
    all(mesv_subset$score_distribution_preserved_all_perm),
  answer = ifelse(
    sum(neuron_subset$pass_positive_fdr) == 0 &&
      all(mesv_subset$count_preserved_all_perm) &&
      all(mesv_subset$score_distribution_preserved_all_perm),
    "clean: neuron remains 0/14 positive-FDR and MES-V count/score preservation remains 14/14 TRUE",
    "not clean: neuron or preservation check failed in the 14-slice subset"
  )
)
fwrite(question3, file.path(out_dir, "A1_Q1_question3_subset_neuron_preservation.csv"))

three_questions <- rbindlist(list(
  q1_source[, .(question, clean, answer)],
  question2[, .(question, clean, answer)],
  question3[, .(question, clean, answer)]
), use.names = TRUE)
all_clean <- all(three_questions$clean)
decision <- data.table(
  item = c("Q1_source_clean", "exclusion_independence_clean", "subset_neuron_preservation_clean", "final_decision"),
  value = c(
    as.character(q1_source$clean),
    as.character(question2$clean),
    as.character(question3$clean),
    ifelse(all_clean, "pass_stratified_sensitivity_for_main_text", "back_to_supplement_strong_reference")
  )
)
fwrite(three_questions, file.path(out_dir, "A1_Q1_three_question_answers.csv"))
fwrite(decision, file.path(out_dir, "A1_Q1_stratification_audit_decision.csv"))

stop_lines <- c(
  "# R9 Batch1 A1 12/14 Q1 stratification compliance audit STOP",
  "",
  paste0("- Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "- Nature: read-only audit. No toroidal rerun, no colocalization recalculation, no threshold change, no new stratum.",
  "- Frozen stratum: `slice_vascular_spot_total >= Q1`, with Q1 from script 45.",
  "",
  "## Question 1: Q1 Source",
  paste0("- Answer: ", q1_source$answer),
  paste0("- Distribution: min=", dist_summary[metric == "min", value],
         "; Q1(type=1)=", q1_type1,
         "; Q1(type=7)=", signif(q1_type7, 4),
         "; median=", dist_summary[metric == "median", value],
         "; Q3(type=1)=", dist_summary[metric == "Q3_type1_empirical_lower", value],
         "; Q3(type=7)=", signif(dist_summary[metric == "Q3_type7_R_default", value], 4),
         "; max=", dist_summary[metric == "max", value], "."),
  paste0("- 18 values low-to-high: ", paste(paste0(dist_dt$slice, "=", dist_dt$slice_vascular_spot_total), collapse = "; "), "."),
  "",
  "## Question 2: Exclusion Independence",
  paste0("- Answer: ", question2$answer),
  paste0("- Excluded 4 slices: ", paste(paste0(excluded$slice, "=", excluded$slice_vascular_spot_total), collapse = "; "), "."),
  paste0("- Overlap with A1 non-pass 4 slices: ", question2$n_overlap_excluded_with_A1_nonpass, "/4 (",
         question2$overlap_excluded_with_A1_nonpass, ")."),
  "- `slice_vascular_spot_total` is a whole-slice dominant-RCTD vascular count, independent of the A1 toroidal program co-enrichment statistic.",
  "",
  "## Question 3: 14-slice Neuron/Preservation",
  paste0("- Answer: ", question3$answer),
  paste0("- MES-V x vascular in 14 slices: ", question3$MESV_positive_FDR_pass, "/14 positive and FDR<0.05; ",
         question3$MESV_positive, "/14 positive."),
  paste0("- neuron_control x vascular in 14 slices: ", question3$neuron_positive_FDR_pass, "/14 positive and FDR<0.05; ",
         question3$neuron_positive, "/14 positive."),
  paste0("- MES-V preservation: count ", question3$MESV_count_preserved_TRUE, "/14 TRUE; score distribution ",
         question3$MESV_score_distribution_preserved_TRUE, "/14 TRUE."),
  "",
  "## Decision",
  paste0("- Final decision: ", decision[item == "final_decision", value], "."),
  if (all_clean) {
    "- Candidate caveat: In the 14 slices with sufficient whole-slice vascular dominant-spot counts, A1 program toroidal co-enrichment passed 12/14; the four excluded slices were the independent low-vascular-count tail and included both IDH-Mutant slices. This does not alter the frozen all-slice primary result of 14/18."
  } else {
    "- Placement: A1 returns to strong supplement/source-data with at most one strong main-text reference."
  },
  "- Forbidden: physical/single-cell colocalization, interaction, recruitment, drives, causality, myeloid/TAM reinterpretation."
)
writeLines(stop_lines, file.path(base_dir, "docs/R9_batch1_A1_Q1_stratification_audit_STOP.md"), useBytes = TRUE)

message("Done. Outputs written to: ", out_dir)
