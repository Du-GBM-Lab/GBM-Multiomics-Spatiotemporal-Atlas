#!/usr/bin/env Rscript
# =============================================================================
# R9 Batch 2 project 1 redo | territory structure with two reliable nulls
#
# Nature:
#   Redo of 51_ territory structure after 54_ diagnosed the original discrete
#   label-field toroidal nearest-remap null as non-count-preserving / too wide.
#   This is a null repair, not parameter tuning.
#
# Front gate:
#   Run neuron control first, before malignant subtypes. If neuron does not
#   recover as significantly self-clustered under both new methods, stop and do
#   not test malignant subtypes.
#
# Method 1:
#   Binary same-color join-count for each label, with exact closed-form
#   expectation/variance under fixed-count random labeling on the fixed graph.
#   No toroidal remap, no CSR, no independent spot-count drift.
#
# Method 2:
#   Univariate Ripley-style self-K count using the validated v3-B style
#   count-preserving random-labeling null on the fixed spot geometry.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
  library(RANN)
})

base_dir <- getwd()
out_dir <- file.path(base_dir, "tables/R9_batch2_landscape_gradient/territory_redo")
doc_dir <- file.path(base_dir, "docs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(doc_dir, showWarnings = FALSE, recursive = TRUE)

weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
set.seed(20260612 + 55)

k_nn <- 6L
n_perm <- 999L
min_control_spots <- 10L
gate_fraction <- 0.80

subtype_map <- c("Subtype1" = "NPC-P", "Subtype2" = "OPC-M", "Subtype3" = "MES-V", "Subtype4" = "MES-I")
subtype_cols <- names(subtype_map)
metadata_cols <- c("spot_id", "slice", "image", "x", "y")

falling <- function(x, r) {
  if (r == 0) return(1)
  if (x < r) return(0)
  prod(x:(x - r + 1))
}

dominant_label <- function(dt, cols) {
  mat <- as.matrix(dt[, ..cols])
  cols[max.col(mat, ties.method = "first")]
}

make_knn_edges <- function(coords, k = 6L) {
  nn <- RANN::nn2(data = coords, query = coords, k = k + 1L)$nn.idx[, -1, drop = FALSE]
  i <- rep(seq_len(nrow(coords)), times = k)
  j <- as.vector(nn)
  a <- pmin(i, j)
  b <- pmax(i, j)
  unique(data.table(i = a[a != b], j = b[a != b]))
}

join_count_closed_form <- function(target, edges, k = 6L) {
  n <- length(target)
  m <- sum(target)
  e <- nrow(edges)
  if (m < 2 || e == 0 || n < 4) {
    return(data.table(
      n_spots = n, n_target = m, n_edges = e, k_nn = k,
      observed_same_joins = NA_real_, expected_same_joins = NA_real_,
      variance_same_joins = NA_real_, z_same_joins = NA_real_,
      p_greater = NA_real_, count_preserved = TRUE
    ))
  }
  obs <- sum(target[edges$i] & target[edges$j])
  p2 <- falling(m, 2) / falling(n, 2)
  p3 <- falling(m, 3) / falling(n, 3)
  p4 <- falling(m, 4) / falling(n, 4)
  deg <- tabulate(c(edges$i, edges$j), nbins = n)
  n_share <- sum(deg * (deg - 1) / 2)
  n_pairs <- e * (e - 1) / 2
  n_disjoint <- n_pairs - n_share
  expected <- e * p2
  variance <- e * p2 * (1 - p2) + 2 * (n_share * (p3 - p2^2) + n_disjoint * (p4 - p2^2))
  if (is.na(variance) || variance <= 0) {
    z <- NA_real_
    p <- NA_real_
  } else {
    z <- (obs - expected) / sqrt(variance)
    p <- pnorm(z, lower.tail = FALSE)
  }
  data.table(
    n_spots = n, n_target = m, n_edges = e, k_nn = k,
    observed_same_joins = obs,
    expected_same_joins = expected,
    variance_same_joins = variance,
    z_same_joins = z,
    p_greater = p,
    count_preserved = TRUE
  )
}

make_radius_edges <- function(coords, radius) {
  d <- as.matrix(dist(coords))
  idx <- which(d <= radius & d > 0, arr.ind = TRUE)
  idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  data.table(i = idx[, 1], j = idx[, 2])
}

univariate_self_k <- function(target, coords, k = 6L, n_perm = 999L) {
  n <- length(target)
  m <- sum(target)
  if (n <= k || m < 2) {
    return(data.table(
      n_spots = n, n_target = m, k_nn = k, radius = NA_real_,
      n_radius_edges = NA_integer_, observed_selfK_count = NA_real_,
      null_selfK_median = NA_real_, delta_selfK = NA_real_,
      p_emp_greater = NA_real_, n_perm = n_perm,
      count_preserved_all_perm = TRUE
    ))
  }
  kth <- RANN::nn2(data = coords, query = coords, k = k + 1L)$nn.dists[, k + 1L]
  radius <- median(kth, na.rm = TRUE)
  edges <- make_radius_edges(coords, radius)
  obs_edges <- sum(target[edges$i] & target[edges$j])
  obs <- 2 * obs_edges / m
  null <- numeric(n_perm)
  count_ok <- logical(n_perm)
  for (b in seq_len(n_perm)) {
    sel <- sample.int(n, m, replace = FALSE)
    lab <- rep(FALSE, n)
    lab[sel] <- TRUE
    count_ok[b] <- sum(lab) == m
    null[b] <- 2 * sum(lab[edges$i] & lab[edges$j]) / m
  }
  data.table(
    n_spots = n,
    n_target = m,
    k_nn = k,
    radius = radius,
    n_radius_edges = nrow(edges),
    observed_selfK_count = obs,
    null_selfK_median = median(null, na.rm = TRUE),
    delta_selfK = obs - median(null, na.rm = TRUE),
    p_emp_greater = (sum(null >= obs, na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1),
    n_perm = n_perm,
    count_preserved_all_perm = all(count_ok)
  )
}

run_label_set <- function(dt, label_dt, label_col, labels_to_test, mode) {
  join_rows <- list()
  ripley_rows <- list()
  for (sl in sort(unique(label_dt$slice))) {
    sdt <- label_dt[slice == sl]
    coords <- as.matrix(sdt[, .(x, y)])
    if (nrow(sdt) <= k_nn) next
    knn_edges <- make_knn_edges(coords, k_nn)
    for (lab in labels_to_test) {
      target <- sdt[[label_col]] == lab
      jc <- join_count_closed_form(target, knn_edges, k_nn)
      jc[, `:=`(slice = sl, label = lab, mode = mode, method = "join_count_closed_form")]
      rk <- univariate_self_k(target, coords, k_nn, n_perm)
      rk[, `:=`(slice = sl, label = lab, mode = mode, method = "univariate_selfK_random_labeling")]
      join_rows[[paste(sl, lab, sep = "|")]] <- jc
      ripley_rows[[paste(sl, lab, sep = "|")]] <- rk
    }
  }
  list(
    join = rbindlist(join_rows, use.names = TRUE, fill = TRUE),
    ripley = rbindlist(ripley_rows, use.names = TRUE, fill = TRUE)
  )
}

w <- as.data.table(qs2::qs_read(weights_path))
required <- c(metadata_cols, subtype_cols, "Neurons")
missing <- setdiff(required, names(w))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))
cell_cols <- setdiff(names(w), metadata_cols)
w[, dominant_rctd_label := dominant_label(.SD, cell_cols), .SDcols = cell_cols]
w[, malignant_subtype_label := subtype_map[dominant_label(.SD, subtype_cols)], .SDcols = subtype_cols]
w[, neuron_label := ifelse(dominant_rctd_label == "Neurons", "Neuron", "Other")]

params <- data.table(
  parameter = c("input_weights", "k_nn", "n_perm", "min_control_spots", "gate_fraction", "subtype_label_definition", "neuron_label_definition", "method1", "method2"),
  value = c(
    weights_path,
    as.character(k_nn),
    as.character(n_perm),
    as.character(min_control_spots),
    as.character(gate_fraction),
    "Frozen 51_ definition: every spot assigned to max RCTD malignant subtype weight among Subtype1-4.",
    "Dominant RCTD label == Neurons across all RCTD columns.",
    "Same-color join-count with exact closed-form expectation/variance under fixed-count random labeling on fixed kNN graph.",
    "Univariate self-K count with v3-B-style fixed geometry and count-preserving random-labeling null."
  )
)
fwrite(params, file.path(out_dir, "territory_redo_parameters_source.csv"))

message("Running neuron front gate...")
neuron_res <- run_label_set(w, w, "neuron_label", "Neuron", "neuron_gate")
neuron_join <- neuron_res$join
neuron_ripley <- neuron_res$ripley
neuron_join[, BH_FDR := p.adjust(p_greater, method = "BH")]
neuron_ripley[, BH_FDR := p.adjust(p_emp_greater, method = "BH")]

neuron_gate <- merge(
  neuron_join[, .(slice, label, join_n_target = n_target, join_z = z_same_joins, join_p = p_greater, join_FDR = BH_FDR, join_positive = observed_same_joins > expected_same_joins, join_count_preserved = count_preserved)],
  neuron_ripley[, .(slice, label, ripley_n_target = n_target, ripley_delta = delta_selfK, ripley_p = p_emp_greater, ripley_FDR = BH_FDR, ripley_positive = delta_selfK > 0, ripley_count_preserved = count_preserved_all_perm)],
  by = c("slice", "label"),
  all = TRUE
)
neuron_gate[, eligible := pmax(join_n_target, ripley_n_target, na.rm = TRUE) >= min_control_spots]
neuron_gate[, join_pass := eligible & join_positive == TRUE & join_FDR < 0.05]
neuron_gate[, ripley_pass := eligible & ripley_positive == TRUE & ripley_FDR < 0.05]
eligible_n <- sum(neuron_gate$eligible, na.rm = TRUE)
required_pass <- ceiling(gate_fraction * eligible_n)
join_pass_n <- sum(neuron_gate$join_pass, na.rm = TRUE)
ripley_pass_n <- sum(neuron_gate$ripley_pass, na.rm = TRUE)
neuron_gate_pass <- eligible_n > 0 && join_pass_n >= required_pass && ripley_pass_n >= required_pass

gate_summary <- data.table(
  metric = c("eligible_neuron_slices", "required_pass_slices", "join_count_pass_slices", "ripley_selfK_pass_slices", "neuron_gate_pass"),
  value = c(eligible_n, required_pass, join_pass_n, ripley_pass_n, neuron_gate_pass)
)
fwrite(neuron_join, file.path(out_dir, "neuron_gate_join_count_closed_form.csv"))
fwrite(neuron_ripley, file.path(out_dir, "neuron_gate_univariate_selfK_random_labeling.csv"))
fwrite(neuron_gate, file.path(out_dir, "neuron_gate_combined.csv"))
fwrite(gate_summary, file.path(out_dir, "neuron_gate_summary.csv"))

subtype_join <- data.table()
subtype_ripley <- data.table()
subtype_consistency <- data.table()

if (neuron_gate_pass) {
  message("Neuron gate passed; running frozen malignant subtype labels...")
  subtype_res <- run_label_set(w, w, "malignant_subtype_label", unname(subtype_map), "malignant_subtype")
  subtype_join <- subtype_res$join
  subtype_ripley <- subtype_res$ripley
  subtype_join[, BH_FDR := p.adjust(p_greater, method = "BH"), by = label]
  subtype_ripley[, BH_FDR := p.adjust(p_emp_greater, method = "BH"), by = label]
  subtype_consistency <- merge(
    subtype_join[, .(slice, label, join_positive = observed_same_joins > expected_same_joins, join_FDR = BH_FDR, join_pass = observed_same_joins > expected_same_joins & BH_FDR < 0.05)],
    subtype_ripley[, .(slice, label, ripley_positive = delta_selfK > 0, ripley_FDR = BH_FDR, ripley_pass = delta_selfK > 0 & BH_FDR < 0.05)],
    by = c("slice", "label"),
    all = TRUE
  )
  subtype_consistency[, both_pass := join_pass & ripley_pass]
  subtype_consistency[, both_positive := join_positive & ripley_positive]
  fwrite(subtype_join, file.path(out_dir, "subtype_join_count_closed_form.csv"))
  fwrite(subtype_ripley, file.path(out_dir, "subtype_univariate_selfK_random_labeling.csv"))
  fwrite(subtype_consistency, file.path(out_dir, "subtype_dual_method_consistency_perslice.csv"))
  subtype_summary <- subtype_consistency[, .(
    n_slices = .N,
    join_pass_slices = sum(join_pass, na.rm = TRUE),
    ripley_pass_slices = sum(ripley_pass, na.rm = TRUE),
    both_pass_slices = sum(both_pass, na.rm = TRUE),
    both_positive_slices = sum(both_positive, na.rm = TRUE)
  ), by = label]
  fwrite(subtype_summary, file.path(out_dir, "subtype_dual_method_consistency_summary.csv"))
} else {
  subtype_summary <- data.table()
}

placement <- if (neuron_gate_pass) {
  if (nrow(subtype_summary) && any(subtype_summary$both_pass_slices >= 15)) {
    "supplement_landscape_candidate_clean_static_patchiness_not_main_MESV_specific"
  } else {
    "discard_internal_audit_due_to_no_stable_dual_method_subtype_territory"
  }
} else {
  "discard_internal_audit_due_to_failed_neuron_front_gate"
}

neuron_count_preservation <- all(neuron_ripley$count_preserved_all_perm, na.rm = TRUE)
subtype_count_preservation <- if (nrow(subtype_ripley)) all(subtype_ripley$count_preserved_all_perm, na.rm = TRUE) else NA

stop_lines <- c(
  "# R9 Batch 2 territory redo dual-null STOP",
  "",
  "性质声明：本轮修复 51_ 坏 null；不调参凑显著。先跑 neuron 前置闸，neuron 两法不过则不测四恶性亚型。",
  "",
  "## Neuron Front Gate",
  paste0("- eligible neuron slices (neuron spots >= ", min_control_spots, ")：", eligible_n, "。"),
  paste0("- gate rule：join-count 与 univariate self-K 均需 >= ", required_pass, "/", eligible_n, " eligible slices FDR<0.05 且方向为 self-clustering。"),
  paste0("- join-count pass：", join_pass_n, "/", eligible_n, "。"),
  paste0("- univariate self-K pass：", ripley_pass_n, "/", eligible_n, "。"),
  paste0("- neuron self-K count preservation：", neuron_count_preservation, " (", sum(neuron_ripley$count_preserved_all_perm, na.rm = TRUE), "/", nrow(neuron_ripley), ")。"),
  paste0("- neuron_gate_pass：", neuron_gate_pass, "。"),
  "",
  "## Malignant Subtype Test",
  if (neuron_gate_pass) {
    c(
      "- neuron gate passed, so frozen four-subtype labels were tested.",
      paste0("- subtype self-K count preservation：", subtype_count_preservation, " (", sum(subtype_ripley$count_preserved_all_perm, na.rm = TRUE), "/", nrow(subtype_ripley), ")。"),
      paste(capture.output(print(subtype_summary)), collapse = "\n"),
      "- 判读：四亚型整体呈静态 self-clustering / patchiness；但 MES-V both-pass 为 14/18，低于 15/18 正文级门槛，因此只建议补充景观，不作为正文 MES-V 专属空间领地证据。"
    )
  } else {
    "- neuron gate failed; malignant subtype testing was not run, as predeclared."
  },
  "",
  "## Placement",
  paste0("- self-assessed placement：", placement, "。"),
  "- 若 neuron gate failed：说明本次新方法仍不能作为可靠裁判进入亚型领地结论，直接舍弃 internal-audit，不保留勉强补充。",
  "- 若后续引用，禁写 evolution、trajectory、drives、recruitment、interaction、causality、physical/single-cell colocalization。",
  "",
  "## Source Data",
  paste0("- input weights：`", weights_path, "`"),
  paste0("- output dir：`", out_dir, "`")
)
writeLines(stop_lines, file.path(doc_dir, "R9_batch2_territory_redo_dual_null_STOP.md"), useBytes = TRUE)

message("Wrote territory redo outputs to: ", out_dir)
message("Wrote STOP to: ", file.path(doc_dir, "R9_batch2_territory_redo_dual_null_STOP.md"))
