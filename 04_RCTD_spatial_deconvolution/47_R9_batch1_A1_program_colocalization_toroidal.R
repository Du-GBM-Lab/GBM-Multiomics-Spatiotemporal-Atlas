#!/usr/bin/env Rscript
# =============================================================================
# R9 Batch1 | A1 program colocalization toroidal-shift primary test
#
# Nature of this test:
#   This is the only pre-declared A1 relation-level primary test. It tests
#   whether the continuous MES-V program score and vascular program score are
#   spatially co-enriched under a self-autocorrelation-preserving toroidal-shift
#   null. The statistic, null, neighborhood, and promotion gate are locked before
#   seeing results. If this fails the pre-specified gate, A1 remains a strong
#   supplementary/source-data layer and no fourth A1 method is pursued.
#
# Scores:
#   Uses existing per-spot continuous scores from the locked R9 score table:
#   MES-V, vascular, neuron_control. These fields were already used by R9 A2
#   and Stage 3 audits. No new gene set is invented here, and all three scores
#   enter the same continuous-score workflow.
#
# Statistic:
#   k=6 bivariate Moran-style spatial lag. For each fixed spot i, z(A_i) is
#   multiplied by the mean z(vascular) among the 6 nearest vascular score-points.
#   The slice statistic is the mean product across all spots. A is MES-V for the
#   primary test and neuron_control for the control test.
#
# Toroidal null:
#   For each slice, vascular score-coordinate pairs are shifted on the bounding
#   box torus by integer multiples of the median coordinate grid step. No score
#   interpolation, no label permutation, no CSR, no random score shuffling. The
#   shifted vascular point cloud may pass through non-tissue holes; holes are not
#   filled and no points are removed. Spot count and vascular score distribution
#   are therefore preserved exactly in every permutation.
#
# Gate:
#   Candidate main-text only if MES-V x vascular is positive and BH FDR<0.05 in
#   >=15/18 slices, neuron_control x vascular has 0/18 FDR<0.05 positive slices,
#   and cross-slice direction is consistent. Otherwise A1 falls back to strong
#   supplementary/source-data status; relation-level evidence remains v3-B.
#
# Forbidden interpretation:
#   No physical colocalization, single-cell colocalization, interaction,
#   recruitment, drives, causality, or TAM recruitment wording.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(FNN)
})

set.seed(20260612)

base_dir <- getwd()
out_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A1_program_colocalization_toroidal")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

score_path <- file.path(base_dir, "tables/R9_stage3_signature_scores/R9_stage3_IVY_hypoxia_scores_per_spot.csv")
a1_clean_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A1_SVG_cleaned")
a1_provenance_path <- file.path(a1_clean_dir, "A1_target_gene_set_provenance.csv")

if (!file.exists(score_path)) stop("Missing score table: ", score_path)
if (!file.exists(a1_provenance_path)) stop("Missing A1 gene-set provenance: ", a1_provenance_path)

n_perm <- 999L
k_nn <- 6L
required_cols <- c("spot_id", "slice", "x", "y", "MES-V", "vascular", "neuron_control")

params <- data.table(
  parameter = c(
    "test_nature",
    "score_source",
    "primary_pair",
    "control_pair",
    "statistic",
    "neighborhood",
    "null",
    "n_perm",
    "shift_grid_step",
    "hole_handling",
    "empirical_p",
    "fdr_scope",
    "promotion_gate",
    "random_seed"
  ),
  value = c(
    "only pre-declared A1 relation-level primary test; no metric/null switching after run",
    score_path,
    "MES-V x vascular",
    "neuron_control x vascular",
    "mean_i z(A_i) * mean_kNN z(vascular_neighbor), bivariate Moran-style spatial lag",
    "k=6 nearest vascular score-points; observed statistic excludes same coordinate when present by querying k+1 and removing distance-zero self",
    "toroidal shift of vascular score-coordinate pairs on slice bounding-box torus; no label permutation, no CSR, no score shuffling",
    as.character(n_perm),
    "integer multiples of the median positive coordinate step for x and y within each slice",
    "non-rectangular holes are not filled; shifted vascular points retain their scores and may pass through holes; no points are dropped",
    "(1 + number of null statistics >= observed) / (n_perm + 1)",
    "BH correction across the 18 per-slice tests within each pair",
    "MES-V positive and FDR<0.05 in >=15/18 slices; neuron positive and FDR<0.05 in 0/18 slices",
    "20260612"
  )
)
fwrite(params, file.path(out_dir, "A1_program_toroidal_parameters.csv"))

scores <- fread(score_path)
missing_cols <- setdiff(required_cols, names(scores))
if (length(missing_cols)) stop("Missing score columns: ", paste(missing_cols, collapse = ", "))
scores <- scores[, ..required_cols]
if (anyNA(scores)) stop("NA detected in required score/coordinate columns.")

score_provenance <- data.table(
  score = c("MES-V", "vascular", "neuron_control"),
  source_path = score_path,
  source_note = c(
    "locked R9 per-spot MES-V continuous score used in A2/Stage3",
    "locked R9 per-spot vascular continuous score used in A2/Stage3",
    "locked R9 per-spot neuron_control continuous score used in A2/Stage3"
  ),
  a1_marker_provenance_context = c(
    a1_provenance_path,
    a1_provenance_path,
    "neuron_control is not recompiled as a marker list in this script; existing locked continuous score is used to avoid improvising a neuron gene set"
  )
)
fwrite(score_provenance, file.path(out_dir, "A1_program_toroidal_score_provenance.csv"))

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

positive_step <- function(x) {
  ux <- sort(unique(x))
  dx <- diff(ux)
  dx <- dx[dx > 0]
  if (!length(dx)) return(1)
  as.numeric(stats::median(dx))
}

wrap_shift <- function(coords, dx, dy, xmin, ymin, width, height) {
  data.table(
    x = ((coords$x - xmin + dx) %% width) + xmin,
    y = ((coords$y - ymin + dy) %% height) + ymin
  )
}

calc_stat <- function(query_coords, shifted_vascular_coords, z_A, z_vascular, k = 6L) {
  nn <- FNN::get.knnx(
    data = as.matrix(shifted_vascular_coords[, .(x, y)]),
    query = as.matrix(query_coords[, .(x, y)]),
    k = k
  )
  lag_v <- rowMeans(matrix(z_vascular[nn$nn.index], nrow = nrow(nn$nn.index)))
  mean(z_A * lag_v, na.rm = TRUE)
}

calc_stat_observed <- function(coords, z_A, z_vascular, k = 6L) {
  nn <- FNN::get.knnx(
    data = as.matrix(coords[, .(x, y)]),
    query = as.matrix(coords[, .(x, y)]),
    k = k + 1L
  )
  idx <- nn$nn.index
  dist <- nn$nn.dist
  keep_idx <- matrix(NA_integer_, nrow = nrow(idx), ncol = k)
  for (i in seq_len(nrow(idx))) {
    candidates <- idx[i, dist[i, ] > 1e-9]
    if (length(candidates) < k) candidates <- idx[i, seq_len(k)]
    keep_idx[i, ] <- candidates[seq_len(k)]
  }
  lag_v <- rowMeans(matrix(z_vascular[keep_idx], nrow = nrow(keep_idx)))
  mean(z_A * lag_v, na.rm = TRUE)
}

emp_p_greater <- function(obs, null_vec) {
  null_vec <- null_vec[is.finite(null_vec)]
  (1 + sum(null_vec >= obs)) / (length(null_vec) + 1)
}

slice_levels <- sort(unique(scores$slice))
pair_map <- data.table(
  pair = c("MESV_vascular", "neuron_vascular"),
  response_score = c("MES-V", "neuron_control"),
  predictor_score = c("vascular", "vascular"),
  role = c("primary", "neuron_control")
)

per_slice_rows <- list()
null_audit_rows <- list()
null_store_rows <- list()
ri <- 0L
ai <- 0L
ni <- 0L

for (sl in slice_levels) {
  sd <- scores[slice == sl]
  coords <- sd[, .(x = as.numeric(x), y = as.numeric(y))]
  n <- nrow(sd)
  x_step <- positive_step(coords$x)
  y_step <- positive_step(coords$y)
  xmin <- min(coords$x)
  ymin <- min(coords$y)
  width <- (max(coords$x) - xmin) + x_step
  height <- (max(coords$y) - ymin) + y_step
  if (!is.finite(width) || width <= 0 || !is.finite(height) || height <= 0) {
    stop("Invalid torus bounds for slice: ", sl)
  }

  z_vascular <- zscore(sd[["vascular"]])
  z_list <- list(
    `MES-V` = zscore(sd[["MES-V"]]),
    neuron_control = zscore(sd[["neuron_control"]])
  )

  obs_by_pair <- pair_map[, .(
    pair,
    response_score,
    predictor_score,
    role,
    obs_stat = vapply(response_score, function(resp) {
      calc_stat_observed(coords, z_list[[resp]], z_vascular, k = k_nn)
    }, numeric(1))
  )]

  null_mat <- matrix(NA_real_, nrow = n_perm, ncol = nrow(pair_map))
  colnames(null_mat) <- pair_map$pair
  count_ok <- logical(n_perm)
  dist_ok <- logical(n_perm)
  shift_rows <- integer(n_perm)
  shift_cols <- integer(n_perm)
  null_score_sum <- numeric(n_perm)
  null_score_sumsq <- numeric(n_perm)

  for (p in seq_len(n_perm)) {
    repeat {
      sx <- sample.int(max(2L, ceiling(width / x_step)), 1L) - 1L
      sy <- sample.int(max(2L, ceiling(height / y_step)), 1L) - 1L
      if (!(sx == 0L && sy == 0L)) break
    }
    shift_rows[p] <- sy
    shift_cols[p] <- sx
    shifted_coords <- wrap_shift(coords, dx = sx * x_step, dy = sy * y_step,
                                 xmin = xmin, ymin = ymin, width = width, height = height)
    count_ok[p] <- nrow(shifted_coords) == n && length(z_vascular) == n
    dist_ok[p] <- isTRUE(all.equal(sort(z_vascular), sort(z_vascular), tolerance = 0))
    null_score_sum[p] <- sum(z_vascular)
    null_score_sumsq[p] <- sum(z_vascular^2)

    for (j in seq_len(nrow(pair_map))) {
      resp <- pair_map$response_score[j]
      null_mat[p, j] <- calc_stat(coords, shifted_coords, z_list[[resp]], z_vascular, k = k_nn)
    }
  }

  for (j in seq_len(nrow(pair_map))) {
    ri <- ri + 1L
    nul <- null_mat[, j]
    obs <- obs_by_pair[j, obs_stat]
    per_slice_rows[[ri]] <- data.table(
      slice = sl,
      pair = pair_map$pair[j],
      role = pair_map$role[j],
      response_score = pair_map$response_score[j],
      predictor_score = pair_map$predictor_score[j],
      n_spots = n,
      k = k_nn,
      obs_stat = obs,
      null_median = median(nul, na.rm = TRUE),
      null_q25 = as.numeric(quantile(nul, 0.25, na.rm = TRUE)),
      null_q75 = as.numeric(quantile(nul, 0.75, na.rm = TRUE)),
      delta_vs_null_median = obs - median(nul, na.rm = TRUE),
      emp_p_greater = emp_p_greater(obs, nul),
      n_perm = n_perm,
      torus_width = width,
      torus_height = height,
      x_step = x_step,
      y_step = y_step,
      count_preserved_all_perm = all(count_ok),
      score_distribution_preserved_all_perm = all(abs(null_score_sum - sum(z_vascular)) < 1e-8) &&
        all(abs(null_score_sumsq - sum(z_vascular^2)) < 1e-8)
    )
  }

  ai <- ai + 1L
  null_audit_rows[[ai]] <- data.table(
    slice = sl,
    n_spots = n,
    n_perm = n_perm,
    k = k_nn,
    x_step = x_step,
    y_step = y_step,
    torus_width = width,
    torus_height = height,
    count_preserved_all_perm = all(count_ok),
    score_distribution_preserved_all_perm = all(abs(null_score_sum - sum(z_vascular)) < 1e-8) &&
      all(abs(null_score_sumsq - sum(z_vascular^2)) < 1e-8),
    zero_shift_used = any(shift_rows == 0L & shift_cols == 0L),
    unique_shifts_used = uniqueN(paste(shift_cols, shift_rows, sep = "_"))
  )

  for (j in seq_len(nrow(pair_map))) {
    ni <- ni + 1L
    null_store_rows[[ni]] <- data.table(
      slice = sl,
      pair = pair_map$pair[j],
      null_stat_mean = mean(null_mat[, j], na.rm = TRUE),
      null_stat_sd = sd(null_mat[, j], na.rm = TRUE),
      null_stat_min = min(null_mat[, j], na.rm = TRUE),
      null_stat_max = max(null_mat[, j], na.rm = TRUE)
    )
  }
}

per_slice <- rbindlist(per_slice_rows, use.names = TRUE)
per_slice[, fdr_bh_by_pair := p.adjust(emp_p_greater, method = "BH"), by = pair]
per_slice[, positive := delta_vs_null_median > 0]
per_slice[, pass_positive_fdr := positive & fdr_bh_by_pair < 0.05]
fwrite(per_slice, file.path(out_dir, "A1_program_toroidal_perslice.csv"))

null_audit <- rbindlist(null_audit_rows, use.names = TRUE)
fwrite(null_audit, file.path(out_dir, "A1_program_toroidal_null_audit.csv"))
fwrite(rbindlist(null_store_rows, use.names = TRUE), file.path(out_dir, "A1_program_toroidal_null_distribution_summary.csv"))

summary_dt <- per_slice[, .(
  n_slices = .N,
  n_positive = sum(positive, na.rm = TRUE),
  n_pass_positive_fdr = sum(pass_positive_fdr, na.rm = TRUE),
  median_obs_stat = median(obs_stat, na.rm = TRUE),
  median_null = median(null_median, na.rm = TRUE),
  median_delta_vs_null = median(delta_vs_null_median, na.rm = TRUE),
  q25_delta_vs_null = as.numeric(quantile(delta_vs_null_median, 0.25, na.rm = TRUE)),
  q75_delta_vs_null = as.numeric(quantile(delta_vs_null_median, 0.75, na.rm = TRUE)),
  median_emp_p = median(emp_p_greater, na.rm = TRUE),
  median_fdr = median(fdr_bh_by_pair, na.rm = TRUE),
  all_count_preserved = all(count_preserved_all_perm),
  all_score_distribution_preserved = all(score_distribution_preserved_all_perm)
), by = .(pair, role, response_score, predictor_score)]

mesv_pass <- summary_dt[pair == "MESV_vascular", n_pass_positive_fdr] >= 15L &&
  summary_dt[pair == "MESV_vascular", n_positive] >= 15L
neuron_clean <- summary_dt[pair == "neuron_vascular", n_pass_positive_fdr] == 0L
count_ok_all <- all(null_audit$count_preserved_all_perm) && all(null_audit$score_distribution_preserved_all_perm)

gate <- data.table(
  gate_item = c(
    "MES-V x vascular positive FDR slices >=15",
    "neuron_control x vascular positive FDR slices ==0",
    "count and score distribution preservation",
    "overall_A1_primary_gate"
  ),
  pass = c(
    mesv_pass,
    neuron_clean,
    count_ok_all,
    mesv_pass && neuron_clean && count_ok_all
  ),
  value = c(
    paste0(summary_dt[pair == "MESV_vascular", n_pass_positive_fdr], "/18 pass; ",
           summary_dt[pair == "MESV_vascular", n_positive], "/18 positive"),
    paste0(summary_dt[pair == "neuron_vascular", n_pass_positive_fdr], "/18 pass; ",
           summary_dt[pair == "neuron_vascular", n_positive], "/18 positive"),
    paste0(sum(null_audit$count_preserved_all_perm & null_audit$score_distribution_preserved_all_perm), "/18 slices preserved"),
    ifelse(mesv_pass && neuron_clean && count_ok_all, "candidate_main_text_pending_user_judgment", "fallback_strong_supplement")
  )
)

fwrite(summary_dt, file.path(out_dir, "A1_program_toroidal_summary.csv"))
fwrite(gate, file.path(out_dir, "A1_program_toroidal_gate.csv"))

decision <- if (gate[gate_item == "overall_A1_primary_gate", pass]) {
  "candidate for main-text overview via A1 primary toroidal relation test; await user/guidance judgment"
} else {
  "fallback to strong supplement/source-data; A1 relation-level claim not promoted, v3-B remains relation-level evidence"
}

stop_lines <- c(
  "# R9 Batch1 A1 program colocalization toroidal-shift STOP",
  "",
  paste0("- Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "- Nature: this is the only pre-declared A1 relation-level primary test. Statistic/null/gate were locked before running; no metric/null switching or slice-picking is allowed after results.",
  "- Scores: existing locked continuous per-spot `MES-V`, `vascular`, and `neuron_control` scores from `R9_stage3_IVY_hypoxia_scores_per_spot.csv`; no new neuron gene set was improvised.",
  "- Statistic: k=6 bivariate Moran-style spatial lag, `mean_i z(A_i) * mean z(vascular_kNN_i)`, with A=`MES-V` or `neuron_control`.",
  "- Null: 999 toroidal shifts of vascular score-coordinate pairs on each slice bounding-box torus; no label permutation, no CSR, no random score shuffling.",
  "- Non-rectangular/hole handling: holes are not filled and no points are dropped; shifted vascular points keep their scores and may pass through holes; count and score distribution preservation are audited.",
  "",
  "## Gate Summary",
  paste0("- MES-V x vascular: ", summary_dt[pair == "MESV_vascular", n_pass_positive_fdr], "/18 positive and BH FDR<0.05; ",
         summary_dt[pair == "MESV_vascular", n_positive], "/18 positive; median delta vs null = ",
         signif(summary_dt[pair == "MESV_vascular", median_delta_vs_null], 4),
         "; median empirical p = ", signif(summary_dt[pair == "MESV_vascular", median_emp_p], 4), "."),
  paste0("- neuron_control x vascular: ", summary_dt[pair == "neuron_vascular", n_pass_positive_fdr], "/18 positive and BH FDR<0.05; ",
         summary_dt[pair == "neuron_vascular", n_positive], "/18 positive; median delta vs null = ",
         signif(summary_dt[pair == "neuron_vascular", median_delta_vs_null], 4),
         "; median empirical p = ", signif(summary_dt[pair == "neuron_vascular", median_emp_p], 4), "."),
  paste0("- Count/score preservation: ", sum(null_audit$count_preserved_all_perm & null_audit$score_distribution_preserved_all_perm), "/18 slices TRUE."),
  paste0("- STOP placement: ", decision, "."),
  "",
  "## Output Files",
  paste0("- ", normalizePath(file.path(out_dir, "A1_program_toroidal_perslice.csv"), winslash = "/")),
  paste0("- ", normalizePath(file.path(out_dir, "A1_program_toroidal_summary.csv"), winslash = "/")),
  paste0("- ", normalizePath(file.path(out_dir, "A1_program_toroidal_gate.csv"), winslash = "/")),
  paste0("- ", normalizePath(file.path(out_dir, "A1_program_toroidal_null_audit.csv"), winslash = "/")),
  paste0("- ", normalizePath(file.path(out_dir, "A1_program_toroidal_parameters.csv"), winslash = "/")),
  paste0("- ", normalizePath(file.path(out_dir, "A1_program_toroidal_score_provenance.csv"), winslash = "/")),
  "",
  "## Wording Boundary",
  "- Allowed only if promoted after judgment: spatially co-enriched / spatial preference / consistent with v3-B.",
  "- Forbidden: physical colocalization, single-cell colocalization, interaction, recruitment, drives, causality, TAM recruitment."
)
writeLines(stop_lines, file.path(base_dir, "docs/R9_batch1_A1_program_colocalization_toroidal_STOP.md"), useBytes = TRUE)

message("Done. Outputs written to: ", out_dir)
