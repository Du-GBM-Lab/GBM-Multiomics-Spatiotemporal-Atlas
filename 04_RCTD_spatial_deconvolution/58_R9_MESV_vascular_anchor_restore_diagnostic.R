#!/usr/bin/env Rscript
# =============================================================================
# R9 | MES-V vascular anchor restore diagnostic
#
# Purpose:
#   Protect the locked v3-B / Figure-A MES-V vascular anchor after 57_ failed
#   under a stricter exploratory niche-neighbor gate. This script distinguishes:
#   (1) harmless stricter gate penalty in 57_, versus
#   (2) a true non-reproducibility problem in the locked source data.
#
# This is a diagnostic only. It does not change the locked v3-B result and does
# not seek a new method.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
})

base_dir <- getwd()
out_dir <- file.path(base_dir, "tables/R9_MESV_niche_neighbors/anchor_restore_diagnostic")
doc_dir <- file.path(base_dir, "docs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(doc_dir, showWarnings = FALSE, recursive = TRUE)

weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
threshold_path <- file.path(base_dir, "tables/C3_niche_preflight/R9_C3_1_neighborhood_niche_scores.csv")
locked_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_C3_4_v3B_MESV_primary_crossK_forest_source.csv")
audit57_path <- file.path(base_dir, "tables/R9_MESV_niche_neighbors/MESV_niche_neighbors_v3B_perslice.csv")

set.seed(20260612 + 58)
n_perm <- 999L
k_values <- c(6L, 12L)
top_fraction <- 0.10

make_hot <- function(x, frac) {
  x >= as.numeric(quantile(x, 1 - frac, na.rm = TRUE))
}

calc_crossK_old <- function(D, A_idx, B_idx, radius) {
  DB <- D[A_idx, B_idx, drop = FALSE]
  mean(rowSums(DB <= radius))
}

emp_p <- function(obs, nul) {
  (1 + sum(nul >= obs, na.rm = TRUE)) / (sum(is.finite(nul)) + 1)
}

weights <- as.data.table(qs2::qs_read(weights_path))
thresholds <- fread(threshold_path)[, .(distance_threshold = median(distance_threshold, na.rm = TRUE)), by = .(slice, k)]
locked <- fread(locked_path)
audit57 <- fread(audit57_path)

feat <- data.table(
  spot_id = weights$spot_id,
  slice = weights$slice,
  x = weights$x,
  y = weights$y,
  `MES-V` = weights$Subtype3,
  vascular = weights$Endothelial + weights[["Mural cells"]]
)

restore_rows <- list()
i <- 0L
for (sl in sort(unique(feat$slice))) {
  dt <- feat[slice == sl]
  D <- as.matrix(dist(as.matrix(dt[, .(x, y)])))
  diag(D) <- 0
  A <- make_hot(dt[["MES-V"]], top_fraction)
  B <- make_hot(dt[["vascular"]], top_fraction)
  n <- nrow(dt)
  nA <- sum(A)
  nB <- sum(B)
  for (kk in k_values) {
    radius <- thresholds[slice == sl & k == kk, distance_threshold][1]
    obs <- calc_crossK_old(D, which(A), which(B), radius)
    null <- numeric(n_perm)
    count_ok <- logical(n_perm)
    for (p in seq_len(n_perm)) {
      Ai <- sample.int(n, nA)
      Bi <- sample.int(n, nB)
      count_ok[p] <- length(Ai) == nA && length(Bi) == nB
      null[p] <- calc_crossK_old(D, Ai, Bi, radius)
    }
    i <- i + 1L
    restore_rows[[i]] <- data.table(
      slice = sl,
      k = kk,
      response = "MES-V",
      predictor = "vascular",
      top_rule = "top10",
      top_fraction = top_fraction,
      n_spots = n,
      n_response_hot = nA,
      n_predictor_hot = nB,
      radius = radius,
      observed_crossK_count_restore = obs,
      null_crossK_median_restore = median(null),
      delta_crossK_restore = obs - median(null),
      p_emp_crossK_restore = emp_p(obs, null),
      emp_sig_restore = emp_p(obs, null) < 0.05,
      count_preserved_all_perm_restore = all(count_ok)
    )
  }
  message("restore done ", sl)
}

restore <- rbindlist(restore_rows)
locked_vasc <- locked[row_type == "slice" & response == "MES-V" & predictor == "vascular" & top_rule == "top10",
  .(slice, k, locked_delta_crossK = delta_crossK, locked_p_emp = p_emp_crossK_greater,
    locked_emp_sig = emp_sig, locked_n_response_hot = n_response_hot,
    locked_n_predictor_hot = n_predictor_hot, locked_count_preserved = count_preserved_all_perm)
]
audit57_vasc <- audit57[response == "MES-V" & predictor == "vascular" & top_rule == "top10",
  .(slice, k, audit57_delta_crossK = delta_crossK, audit57_p_emp = p_emp_crossK_greater,
    audit57_FDR = FDR_crossK_by_slice_response_k, audit57_crossK_pass_FDR = crossK_pass_FDR,
    audit57_delta_nn_closer = delta_nn_closer, audit57_p_nn = p_emp_nn_less,
    audit57_FDR_nn = FDR_nn_by_slice_response_k, audit57_nn_pass_FDR = nn_pass_FDR)
]

compare <- Reduce(function(x, y) merge(x, y, by = c("slice", "k"), all = TRUE),
  list(restore, locked_vasc, audit57_vasc)
)
compare[, restore_vs_locked_delta_absdiff := abs(delta_crossK_restore - locked_delta_crossK)]
compare[, restore_vs_57_delta_absdiff := abs(delta_crossK_restore - audit57_delta_crossK)]
compare[, locked_vs_57_delta_absdiff := abs(locked_delta_crossK - audit57_delta_crossK)]
compare[, locked_restore_same_emp_sig := emp_sig_restore == locked_emp_sig]

summary_restore <- compare[, .(
  n_slices = .N,
  restore_positive = sum(delta_crossK_restore > 0, na.rm = TRUE),
  restore_emp_sig = sum(emp_sig_restore, na.rm = TRUE),
  locked_positive = sum(locked_delta_crossK > 0, na.rm = TRUE),
  locked_emp_sig = sum(locked_emp_sig, na.rm = TRUE),
  audit57_positive = sum(audit57_delta_crossK > 0, na.rm = TRUE),
  audit57_emp_sig = sum(audit57_p_emp < 0.05, na.rm = TRUE),
  audit57_FDR_sig = sum(audit57_crossK_pass_FDR, na.rm = TRUE),
  median_restore_p = median(p_emp_crossK_restore, na.rm = TRUE),
  median_locked_p = median(locked_p_emp, na.rm = TRUE),
  median_audit57_p = median(audit57_p_emp, na.rm = TRUE),
  max_absdiff_restore_locked_delta = max(restore_vs_locked_delta_absdiff, na.rm = TRUE),
  median_absdiff_restore_locked_delta = median(restore_vs_locked_delta_absdiff, na.rm = TRUE),
  count_preserved_restore_all = all(count_preserved_all_perm_restore)
), by = k]

drop_diagnosis <- data.table(
  gate_component = c(
    "locked_source_empirical_crossK",
    "current_restore_old_conditions_empirical_crossK",
    "57_raw_empirical_crossK_with_overlap_exclusion",
    "57_multiclass_FDR_crossK",
    "57_nearest_distance_FDR",
    "57_combination_gate_crossK_and_NN"
  ),
  k6_pass_slices = c(
    summary_restore[k == 6, locked_emp_sig],
    summary_restore[k == 6, restore_emp_sig],
    summary_restore[k == 6, audit57_emp_sig],
    summary_restore[k == 6, audit57_FDR_sig],
    compare[k == 6, sum(audit57_nn_pass_FDR, na.rm = TRUE)],
    compare[k == 6, sum(audit57_crossK_pass_FDR & audit57_nn_pass_FDR, na.rm = TRUE)]
  ),
  k12_pass_slices = c(
    summary_restore[k == 12, locked_emp_sig],
    summary_restore[k == 12, restore_emp_sig],
    summary_restore[k == 12, audit57_emp_sig],
    summary_restore[k == 12, audit57_FDR_sig],
    compare[k == 12, sum(audit57_nn_pass_FDR, na.rm = TRUE)],
    compare[k == 12, sum(audit57_crossK_pass_FDR & audit57_nn_pass_FDR, na.rm = TRUE)]
  )
)

locked_safe <- all(summary_restore$restore_emp_sig >= 15) &&
  all(summary_restore$locked_emp_sig >= 15) &&
  all(summary_restore$max_absdiff_restore_locked_delta < 0.10)

interpretation <- if (locked_safe) {
  "possible_1_confirmed_57_gate_penalty_not_locked_A_problem"
} else {
  "red_flag_locked_anchor_requires_dedicated_v3B_audit"
}

fwrite(restore, file.path(out_dir, "MESV_vascular_restore_old_conditions_current_data.csv"))
fwrite(compare, file.path(out_dir, "MESV_vascular_anchor_locked_vs_57_compare_perslice.csv"))
fwrite(summary_restore, file.path(out_dir, "MESV_vascular_anchor_restore_summary.csv"))
fwrite(drop_diagnosis, file.path(out_dir, "MESV_vascular_anchor_gate_drop_diagnosis.csv"))

stop_lines <- c(
  "# R9 MES-V vascular anchor restore diagnostic STOP",
  "",
  "性质声明：本诊断只保护已锁 v3-B / 档 A 主证据，区分 57_ 更严组合闸导致的无害掉线，还是档 A 本身不稳。它不新增方法、不替换档 A。",
  "",
  "## Restored Original Conditions",
  "- 条件：MES-V response + vascular predictor only；top10；cross-K only；no multi-predictor FDR；no focal-overlap exclusion；no NN combination gate；random-labeling null 999；current RCTD weights。",
  paste0("- k=6 restore empirical pass：", summary_restore[k == 6, restore_emp_sig], "/18；locked source empirical pass：", summary_restore[k == 6, locked_emp_sig], "/18；57 raw empirical cross-K：", summary_restore[k == 6, audit57_emp_sig], "/18；57 FDR cross-K：", summary_restore[k == 6, audit57_FDR_sig], "/18；57 NN FDR：", drop_diagnosis[gate_component == "57_nearest_distance_FDR", k6_pass_slices], "/18。"),
  paste0("- k=12 restore empirical pass：", summary_restore[k == 12, restore_emp_sig], "/18；locked source empirical pass：", summary_restore[k == 12, locked_emp_sig], "/18；57 raw empirical cross-K：", summary_restore[k == 12, audit57_emp_sig], "/18；57 FDR cross-K：", summary_restore[k == 12, audit57_FDR_sig], "/18；57 NN FDR：", drop_diagnosis[gate_component == "57_nearest_distance_FDR", k12_pass_slices], "/18。"),
  paste0("- restore vs locked max abs delta difference：k=6 ", signif(summary_restore[k == 6, max_absdiff_restore_locked_delta], 4), "; k=12 ", signif(summary_restore[k == 12, max_absdiff_restore_locked_delta], 4), "。"),
  paste0("- restore count preservation：", all(restore$count_preserved_all_perm_restore), "。"),
  "",
  "## Gate Drop Diagnosis",
  paste(capture.output(print(drop_diagnosis)), collapse = "\n"),
  "",
  "## Interpretation",
  paste0("- call：", interpretation, "。"),
  if (locked_safe) {
    "- 结论：档 A 在原始 v3-B cross-K 条件下由当前数据复跑仍回到 15-17/18，且与已锁 source data 逐片 delta 基本一致。57_ 掉线来自新增组合闸，主要是 focal-overlap exclusion/多类 FDR 后 cross-K 降至 12/18，以及 NN 指标仅 3/18。正文档 A 不动；57_ 继续 internal-audit。"
  } else {
    "- 结论：还原条件未能回到锁定表现或与锁定 source 不一致；应暂停正文档 A 写法并专门审 v3-B。"
  },
  "",
  "## Source Data",
  paste0("- locked source：`", locked_path, "`"),
  paste0("- 57 source：`", audit57_path, "`"),
  paste0("- output dir：`", out_dir, "`")
)
writeLines(stop_lines, file.path(doc_dir, "R9_MESV_vascular_anchor_restore_diagnostic_STOP.md"), useBytes = TRUE)

message("Wrote anchor restore diagnostic outputs to: ", out_dir)
message("Wrote STOP to: ", file.path(doc_dir, "R9_MESV_vascular_anchor_restore_diagnostic_STOP.md"))
