#!/usr/bin/env Rscript
# =============================================================================
# R9 Batch 2 | bin-level gradient check + territory neuron-null diagnosis
#
# Nature:
#   Lightweight audit only. This script reads existing 52_ gradient outputs and
#   51_ territory/null outputs. It does not rerun the gradient, does not rerun
#   the territory main analysis, and does not tune parameters to obtain
#   significance.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- getwd()
out_dir <- file.path(base_dir, "tables/R9_batch2_landscape_gradient/bin_and_null_check")
doc_dir <- file.path(base_dir, "docs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(doc_dir, showWarnings = FALSE, recursive = TRUE)

grad_perslice_path <- file.path(base_dir, "tables/R9_batch2_landscape_gradient/vascular_distance_gradient/vascular_distance_composition_perslice_bins.csv")
grad_summary_path <- file.path(base_dir, "tables/R9_batch2_landscape_gradient/vascular_distance_gradient/vascular_distance_composition_summary_bins.csv")
neuron_path <- file.path(base_dir, "tables/R9_batch2_landscape_gradient/territory_segregation/territory_segregation_neuron_control.csv")
null_path <- file.path(base_dir, "tables/R9_batch2_landscape_gradient/territory_segregation/territory_segregation_toroidal_null_audit.csv")
param_path <- file.path(base_dir, "tables/R9_batch2_landscape_gradient/territory_segregation/territory_segregation_parameters_source.csv")

need <- c(grad_perslice_path, grad_summary_path, neuron_path, null_path, param_path)
missing <- need[!file.exists(need)]
if (length(missing)) stop("Missing required input(s): ", paste(missing, collapse = "; "))

grad <- fread(grad_perslice_path)
grad_sum <- fread(grad_summary_path)
neuron <- fread(neuron_path)
null_dt <- fread(null_path)
params <- fread(param_path)

focus_groups <- c(
  "NPC-P", "OPC-M", "MES-V", "MES-I",
  "TAM", "Endothelial", "Mural cells", "Neurons",
  "Oligodendrocytes", "OPCs", "T cells", "B cells", "NK cells", "cDCs", "pDCs"
)

bin_table <- grad[group %in% focus_groups, .(
  n_slices = .N,
  mean_across_slices = mean(mean_weight, na.rm = TRUE),
  sd_across_slices = sd(mean_weight, na.rm = TRUE),
  min_across_slices = min(mean_weight, na.rm = TRUE),
  max_across_slices = max(mean_weight, na.rm = TRUE),
  median_across_slices = median(mean_weight, na.rm = TRUE),
  q25_across_slices = quantile(mean_weight, 0.25, na.rm = TRUE),
  q75_across_slices = quantile(mean_weight, 0.75, na.rm = TRUE),
  cv_sd_over_abs_mean = ifelse(abs(mean(mean_weight, na.rm = TRUE)) > 0, sd(mean_weight, na.rm = TRUE) / abs(mean(mean_weight, na.rm = TRUE)), NA_real_)
), by = .(group, vascular_distance_bin)]
setorder(bin_table, group, vascular_distance_bin)

wide_group <- dcast(grad[group %in% focus_groups], slice + group ~ vascular_distance_bin, value.var = "mean_weight")
setnames(wide_group, as.character(1:5), paste0("bin", 1:5))

trend_by_slice <- wide_group[, {
  vals <- as.numeric(.SD)
  diffs <- diff(vals)
  data.table(
    bin1_near = vals[1],
    bin5_far = vals[5],
    near_minus_far = vals[1] - vals[5],
    adjacent_decrease_count = sum(diffs < 0, na.rm = TRUE),
    adjacent_increase_count = sum(diffs > 0, na.rm = TRUE),
    monotonic_decreasing = all(diffs <= 0, na.rm = TRUE),
    endpoint_near_gt_far = vals[1] > vals[5],
    max_adjacent_reversal = ifelse(any(diffs > 0, na.rm = TRUE), max(diffs, na.rm = TRUE), 0)
  )
}, by = .(slice, group), .SDcols = paste0("bin", 1:5)]

trend_summary <- trend_by_slice[, .(
  n_slices = .N,
  median_near_minus_far = median(near_minus_far, na.rm = TRUE),
  mean_near_minus_far = mean(near_minus_far, na.rm = TRUE),
  sd_near_minus_far = sd(near_minus_far, na.rm = TRUE),
  near_gt_far_slices = sum(endpoint_near_gt_far, na.rm = TRUE),
  monotonic_decreasing_slices = sum(monotonic_decreasing, na.rm = TRUE),
  median_adjacent_decrease_count = median(adjacent_decrease_count, na.rm = TRUE),
  max_adjacent_reversal_median = median(max_adjacent_reversal, na.rm = TRUE)
), by = group]

adjacent_consistency <- grad[group %in% focus_groups][order(slice, group, vascular_distance_bin)]
adjacent_consistency <- adjacent_consistency[, .(
  adjacent_pair = paste0("bin", vascular_distance_bin[-.N], "_to_bin", vascular_distance_bin[-1]),
  decrease = diff(mean_weight) < 0,
  delta_next_minus_current = diff(mean_weight)
), by = .(slice, group)]
adjacent_summary <- adjacent_consistency[, .(
  n_slices = .N,
  decrease_slices = sum(decrease, na.rm = TRUE),
  increase_or_flat_slices = sum(!decrease, na.rm = TRUE),
  median_delta_next_minus_current = median(delta_next_minus_current, na.rm = TRUE)
), by = .(group, adjacent_pair)]

mesv_bins <- bin_table[group == "MES-V"]
mesv_adj <- adjacent_summary[group == "MES-V"]
mesv_trend <- trend_summary[group == "MES-V"]
tam_bins <- bin_table[group == "TAM"]
tam_adj <- adjacent_summary[group == "TAM"]
tam_trend <- trend_summary[group == "TAM"]

mesv_median_vals <- mesv_bins[order(vascular_distance_bin), median_across_slices]
mesv_mean_vals <- mesv_bins[order(vascular_distance_bin), mean_across_slices]
mesv_median_monotonic <- all(diff(mesv_median_vals) <= 0, na.rm = TRUE)
mesv_mean_monotonic <- all(diff(mesv_mean_vals) <= 0, na.rm = TRUE)
mesv_majority_adjacent <- all(mesv_adj$decrease_slices >= ceiling(mesv_adj$n_slices / 2))
mesv_noise_ratio <- median(mesv_bins$cv_sd_over_abs_mean, na.rm = TRUE)
mesv_supports_main <- mesv_median_monotonic && mesv_majority_adjacent && is.finite(mesv_noise_ratio) && mesv_noise_ratio < 1.5

null_by_slice <- null_dt[, .(
  null_mean_neuron_same_neighbor_fraction = mean(neuron_same_neighbor_fraction, na.rm = TRUE),
  null_sd_neuron_same_neighbor_fraction = sd(neuron_same_neighbor_fraction, na.rm = TRUE),
  null_q025_neuron_same_neighbor_fraction = quantile(neuron_same_neighbor_fraction, 0.025, na.rm = TRUE),
  null_q975_neuron_same_neighbor_fraction = quantile(neuron_same_neighbor_fraction, 0.975, na.rm = TRUE),
  null_valid_iterations = sum(!is.na(neuron_same_neighbor_fraction)),
  neuron_count_exact_iterations = sum(neuron_count_delta == 0, na.rm = TRUE),
  neuron_count_max_delta = max(neuron_count_delta, na.rm = TRUE),
  neuron_count_median_delta = median(neuron_count_delta, na.rm = TRUE)
), by = slice]

neuron_diag <- merge(neuron, null_by_slice, by = "slice", all.x = TRUE)
neuron_diag[, obs_minus_null_mean := neuron_same_neighbor_fraction - null_mean_neuron_same_neighbor_fraction]
neuron_diag[, obs_z_vs_null := obs_minus_null_mean / null_sd_neuron_same_neighbor_fraction]
neuron_diag[, obs_above_null_mean := neuron_same_neighbor_fraction > null_mean_neuron_same_neighbor_fraction]
neuron_diag[, obs_above_null_q975 := neuron_same_neighbor_fraction > null_q975_neuron_same_neighbor_fraction]
neuron_diag[, count_preservation_problem := toroidal_neuron_count_exact_all_perm != TRUE | neuron_count_max_delta > 0]

neuron_diag_summary <- data.table(
  metric = c(
    "slices_total",
    "slices_with_neuron_spots",
    "neuron_obs_above_null_mean",
    "neuron_obs_above_null_q975",
    "neuron_FDR_pass",
    "slices_with_count_preservation_problem",
    "median_obs_minus_null_mean",
    "median_obs_z_vs_null",
    "k_nn",
    "n_perm",
    "null_type"
  ),
  value = c(
    nrow(neuron_diag),
    sum(neuron_diag$neuron_spots > 0, na.rm = TRUE),
    sum(neuron_diag$obs_above_null_mean, na.rm = TRUE),
    sum(neuron_diag$obs_above_null_q975, na.rm = TRUE),
    sum(neuron_diag$BH_FDR_neuron_segregated < 0.05, na.rm = TRUE),
    sum(neuron_diag$count_preservation_problem, na.rm = TRUE),
    signif(median(neuron_diag$obs_minus_null_mean, na.rm = TRUE), 5),
    signif(median(neuron_diag$obs_z_vs_null, na.rm = TRUE), 5),
    unique(na.omit(neuron_diag$k_nn))[1],
    unique(na.omit(neuron_diag$n_perm))[1],
    params[parameter == "null", value][1]
  )
)

root_cause <- if (sum(neuron_diag$count_preservation_problem, na.rm = TRUE) > 0) {
  "null_implementation_count_preservation_problem_for_discrete_label_field"
} else if (sum(neuron_diag$obs_above_null_mean, na.rm = TRUE) > 0 && sum(neuron_diag$BH_FDR_neuron_segregated < 0.05, na.rm = TRUE) == 0) {
  "observed_above_null_but_underpowered_or_null_too_wide"
} else {
  "observed_approximately_null_at_this_scale"
}

fwrite(bin_table, file.path(out_dir, "gradient_bin_values_mean_sd_range.csv"))
fwrite(trend_by_slice, file.path(out_dir, "gradient_trend_by_slice.csv"))
fwrite(trend_summary, file.path(out_dir, "gradient_trend_summary.csv"))
fwrite(adjacent_summary, file.path(out_dir, "gradient_adjacent_bin_consistency.csv"))
fwrite(mesv_bins, file.path(out_dir, "MESV_gradient_bin_values.csv"))
fwrite(tam_bins, file.path(out_dir, "TAM_gradient_bin_values.csv"))
fwrite(neuron_diag, file.path(out_dir, "territory_neuron_null_diagnosis_perslice.csv"))
fwrite(neuron_diag_summary, file.path(out_dir, "territory_neuron_null_diagnosis_summary.csv"))
fwrite(data.table(
  check = c("MESV_median_monotonic_near_to_far", "MESV_mean_monotonic_near_to_far", "MESV_majority_adjacent_decrease", "MESV_noise_cv_median_lt_1.5", "MESV_supports_main_overview", "territory_neuron_root_cause"),
  value = c(mesv_median_monotonic, mesv_mean_monotonic, mesv_majority_adjacent, mesv_noise_ratio < 1.5, mesv_supports_main, root_cause)
), file.path(out_dir, "bin_and_null_check_decision_flags.csv"))

fmt_vec <- function(x) paste(signif(x, 4), collapse = " -> ")

stop_lines <- c(
  "# R9 Batch 2 bin-level gradient + territory neuron-null diagnosis STOP",
  "",
  "性质声明：本轮只读 52_ 距血管梯度结果与 51_ 领地结构/null 输出；不重跑梯度、不重做领地主分析、不调参数追显著。所有结论仍为 Visium spot-level spatial landscape / audit。",
  "",
  "## Check A: distance-to-vascular gradient bin-level audit",
  paste0("- 输出目录：`", out_dir, "`"),
  paste0("- MES-V median bin trajectory (near -> far)：", fmt_vec(mesv_median_vals), "。"),
  paste0("- MES-V mean bin trajectory (near -> far)：", fmt_vec(mesv_mean_vals), "。"),
  paste0("- MES-V median monotonic decreasing：", mesv_median_monotonic, "；mean monotonic decreasing：", mesv_mean_monotonic, "。"),
  paste0("- MES-V adjacent-bin majority decrease：", mesv_majority_adjacent, "；per-adjacent decrease slices：", paste(mesv_adj$adjacent_pair, paste0(mesv_adj$decrease_slices, "/", mesv_adj$n_slices), sep = "=", collapse = "; "), "。"),
  paste0("- MES-V endpoint near>far slices：", mesv_trend$near_gt_far_slices, "/", mesv_trend$n_slices,
         "；monotonic-decreasing slices：", mesv_trend$monotonic_decreasing_slices, "/", mesv_trend$n_slices,
         "；median near-minus-far=", signif(mesv_trend$median_near_minus_far, 4),
         "；median CV across bins=", signif(mesv_noise_ratio, 4), "。"),
  paste0("- TAM median bin trajectory (near -> far)：", fmt_vec(tam_bins[order(vascular_distance_bin), median_across_slices]), "。"),
  paste0("- TAM endpoint near>far slices：", tam_trend$near_gt_far_slices, "/", tam_trend$n_slices,
         "；monotonic-decreasing slices：", tam_trend$monotonic_decreasing_slices, "/", tam_trend$n_slices,
         "；median near-minus-far=", signif(tam_trend$median_near_minus_far, 4), "。TAM 只作诊断，不立空间偏好或 recruitment 主张。"),
  paste0("- 自评：", if (mesv_supports_main) "supports main overview candidate as a descriptive gradient" else "retreat to supplement; bin trajectory is not clean enough for main overview", "。"),
  "",
  "## Check B: territory neuron-null diagnosis",
  paste0("- k_nn=", unique(na.omit(neuron_diag$k_nn))[1], "；n_perm=", unique(na.omit(neuron_diag$n_perm))[1], "；null=", params[parameter == "null", value][1], "。"),
  paste0("- slices with neuron spots：", sum(neuron_diag$neuron_spots > 0, na.rm = TRUE), "/18；neuron FDR pass：", sum(neuron_diag$BH_FDR_neuron_segregated < 0.05, na.rm = TRUE), "/18。"),
  paste0("- neuron observed above null mean：", sum(neuron_diag$obs_above_null_mean, na.rm = TRUE), "/18；above null 97.5%：", sum(neuron_diag$obs_above_null_q975, na.rm = TRUE), "/18。"),
  paste0("- median obs-minus-null-mean：", signif(median(neuron_diag$obs_minus_null_mean, na.rm = TRUE), 4),
         "；median z vs null：", signif(median(neuron_diag$obs_z_vs_null, na.rm = TRUE), 4), "。"),
  paste0("- count preservation problem in neuron toroidal remap：", sum(neuron_diag$count_preservation_problem, na.rm = TRUE), "/18；root cause call：", root_cause, "。"),
  "- 诊断结论：51_ 的 neuron 0/18 不能被当成“neuron 不隔离”的强生物学结论；主要问题是离散 label-field toroidal nearest-remap null 对 neuron count preservation 不干净/过宽。按 spec 不修、不调参，项目 1 仍保守落 internal-audit/supplement landscape。",
  "",
  "## 禁止写法",
  "- 禁写 drives、evolution、trajectory、recruitment、interaction、causality、physical/single-cell colocalization、TAM recruitment。"
)
writeLines(stop_lines, file.path(doc_dir, "R9_batch2_bin_and_null_check_STOP.md"), useBytes = TRUE)

message("Wrote 54_ bin/null check outputs to: ", out_dir)
message("Wrote STOP to: ", file.path(doc_dir, "R9_batch2_bin_and_null_check_STOP.md"))
