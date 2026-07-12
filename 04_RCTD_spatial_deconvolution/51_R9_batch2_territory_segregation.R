#!/usr/bin/env Rscript
# =============================================================================
# R9 Batch 2 project 1 | malignant subtype territory / spatial segregation
#
# Nature:
#   Descriptive spatial landscape. This script asks whether NPC-P/OPC-M/MES-V/
#   MES-I subtype fields are patch-like or intermixed at the Visium spot level.
#   It does not infer evolution, transition, recruitment, interaction, or
#   causality. Relation-level manuscript claims must return to locked v3-B.
#
# Locked inputs:
#   - RCTD weights from C3_4 local niche
#   - Subtype1=NPC-P, Subtype2=OPC-M, Subtype3=MES-V, Subtype4=MES-I
#
# Null note:
#   A toroidal coordinate-shift label-field null is used as an autocorrelation-
#   preserving audit for same-neighbor enrichment. It shifts the existing label
#   field on the slice bounding box and remaps to nearest observed spots; no CSR
#   or independent label permutation is used.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
  library(RANN)
})

base_dir <- getwd()
out_dir <- file.path(base_dir, "tables/R9_batch2_landscape_gradient/territory_segregation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
set.seed(20260612)
k_nn <- 6L
n_perm <- 99L

subtype_map <- c(
  "Subtype1" = "NPC-P",
  "Subtype2" = "OPC-M",
  "Subtype3" = "MES-V",
  "Subtype4" = "MES-I"
)
subtype_cols <- names(subtype_map)

read_weights <- function(path) {
  if (!file.exists(path)) stop("Missing RCTD weights: ", path)
  as.data.table(qs2::qs_read(path))
}

dominant_label <- function(dt, cols) {
  mat <- as.matrix(dt[, ..cols])
  cols[max.col(mat, ties.method = "first")]
}

shannon_entropy <- function(x) {
  p <- table(x) / length(x)
  -sum(as.numeric(p) * log2(as.numeric(p)))
}

same_stats <- function(labels, nn_idx) {
  n <- length(labels)
  neighbor_labels <- matrix(labels[as.vector(nn_idx)], nrow = nrow(nn_idx), ncol = ncol(nn_idx))
  same_mat <- neighbor_labels == labels
  same_by_spot <- rowMeans(same_mat, na.rm = TRUE)
  entropy_by_spot <- vapply(seq_len(n), function(i) {
    shannon_entropy(labels[nn_idx[i, ]])
  }, numeric(1))
  lab_tab <- prop.table(table(labels))
  expected_same <- sum(as.numeric(lab_tab)^2)
  observed_same <- mean(same_by_spot, na.rm = TRUE)
  data.table(
    observed_same_neighbor_fraction = observed_same,
    expected_same_from_label_frequency = expected_same,
    assortativity_index = ifelse(expected_same < 1, (observed_same - expected_same) / (1 - expected_same), NA_real_),
    median_local_entropy = median(entropy_by_spot, na.rm = TRUE),
    mean_local_entropy = mean(entropy_by_spot, na.rm = TRUE)
  )
}

same_neighbor_fraction <- function(labels, nn_idx) {
  neighbor_labels <- matrix(labels[as.vector(nn_idx)], nrow = nrow(nn_idx), ncol = ncol(nn_idx))
  mean(neighbor_labels == labels, na.rm = TRUE)
}

label_stats_by_class <- function(labels, nn_idx) {
  rbindlist(lapply(sort(unique(labels)), function(lab) {
    idx <- which(labels == lab)
    data.table(
      label = lab,
      n_spots = length(idx),
      same_neighbor_fraction = mean(vapply(idx, function(i) {
        mean(labels[nn_idx[i, ]] == lab, na.rm = TRUE)
      }, numeric(1)), na.rm = TRUE)
    )
  }))
}

neuron_stats <- function(is_neuron, nn_idx) {
  idx <- which(is_neuron)
  if (!length(idx)) {
    return(data.table(neuron_spots = 0L, neuron_same_neighbor_fraction = NA_real_))
  }
  data.table(
    neuron_spots = length(idx),
    neuron_same_neighbor_fraction = mean(rowMeans(matrix(is_neuron[as.vector(nn_idx[idx, , drop = FALSE])], nrow = length(idx)), na.rm = TRUE), na.rm = TRUE)
  )
}

neuron_same_fraction <- function(is_neuron, nn_idx) {
  idx <- which(is_neuron)
  if (!length(idx)) return(NA_real_)
  mean(rowMeans(matrix(is_neuron[as.vector(nn_idx[idx, , drop = FALSE])], nrow = length(idx)), na.rm = TRUE), na.rm = TRUE)
}

toroidal_shift_labels <- function(coords, labels) {
  xr <- range(coords[, 1])
  yr <- range(coords[, 2])
  width <- diff(xr)
  height <- diff(yr)
  if (width == 0 || height == 0) return(labels)
  dx <- runif(1, 0, width)
  dy <- runif(1, 0, height)
  shifted <- cbind(
    ((coords[, 1] - xr[1] + dx) %% width) + xr[1],
    ((coords[, 2] - yr[1] + dy) %% height) + yr[1]
  )
  map <- RANN::nn2(data = shifted, query = coords, k = 1)$nn.idx[, 1]
  labels[map]
}

empirical_right_p <- function(obs, null) {
  (sum(null >= obs, na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1)
}

w <- read_weights(weights_path)
required <- c("spot_id", "slice", "image", "x", "y", subtype_cols, "Neurons")
missing <- setdiff(required, names(w))
if (length(missing)) stop("Missing columns in RCTD weights: ", paste(missing, collapse = ", "))

metadata_cols <- c("spot_id", "slice", "image", "x", "y")
cell_cols <- setdiff(names(w), metadata_cols)

w[, dominant_rctd_label := dominant_label(.SD, cell_cols), .SDcols = cell_cols]
w[, malignant_subtype_label := subtype_map[dominant_label(.SD, subtype_cols)], .SDcols = subtype_cols]
w[, neuron_dominant := dominant_rctd_label == "Neurons"]

slice_rows <- list()
class_rows <- list()
null_rows <- list()
neuron_rows <- list()

for (sl in sort(unique(w$slice))) {
  dt <- w[slice == sl]
  coords <- as.matrix(dt[, .(x, y)])
  if (nrow(dt) <= k_nn) next
  nn_idx <- RANN::nn2(data = coords, query = coords, k = k_nn + 1L)$nn.idx[, -1, drop = FALSE]

  labels <- dt$malignant_subtype_label
  obs <- same_stats(labels, nn_idx)
  obs[, `:=`(
    slice = sl,
    n_spots = nrow(dt),
    k_nn = k_nn,
    n_perm = n_perm
  )]

  class_dt <- label_stats_by_class(labels, nn_idx)
  class_dt[, `:=`(slice = sl, k_nn = k_nn)]

  neuron_obs <- neuron_stats(dt$neuron_dominant, nn_idx)
  neuron_obs[, `:=`(slice = sl, n_spots = nrow(dt), k_nn = k_nn, n_perm = n_perm)]

  null_same <- numeric(n_perm)
  null_neuron <- numeric(n_perm)
  subtype_count_max_abs_delta <- integer(n_perm)
  neuron_count_delta <- integer(n_perm)
  for (b in seq_len(n_perm)) {
    shifted_labels <- toroidal_shift_labels(coords, labels)
    shifted_neuron <- toroidal_shift_labels(coords, dt$neuron_dominant)
    null_same[b] <- same_neighbor_fraction(shifted_labels, nn_idx)
    null_neuron[b] <- neuron_same_fraction(shifted_neuron, nn_idx)
    subtype_count_max_abs_delta[b] <- max(abs(table(factor(shifted_labels, levels = unique(labels))) - table(factor(labels, levels = unique(labels)))))
    neuron_count_delta[b] <- abs(sum(shifted_neuron) - sum(dt$neuron_dominant))
  }

  obs[, `:=`(
    null_median_same_neighbor_fraction = median(null_same, na.rm = TRUE),
    null_q025_same_neighbor_fraction = quantile(null_same, 0.025, na.rm = TRUE),
    null_q975_same_neighbor_fraction = quantile(null_same, 0.975, na.rm = TRUE),
    emp_p_same_neighbor_enriched = empirical_right_p(observed_same_neighbor_fraction, null_same),
    toroidal_label_count_exact_all_perm = all(subtype_count_max_abs_delta == 0L),
    toroidal_label_count_max_abs_delta = max(subtype_count_max_abs_delta)
  )]
  neuron_obs[, `:=`(
    null_median_neuron_same_neighbor_fraction = median(null_neuron, na.rm = TRUE),
    emp_p_neuron_segregated = empirical_right_p(neuron_same_neighbor_fraction, null_neuron),
    toroidal_neuron_count_exact_all_perm = all(neuron_count_delta == 0L),
    toroidal_neuron_count_max_abs_delta = max(neuron_count_delta)
  )]

  null_rows[[sl]] <- data.table(
    slice = sl,
    null_iteration = seq_len(n_perm),
    subtype_same_neighbor_fraction = null_same,
    neuron_same_neighbor_fraction = null_neuron,
    subtype_count_max_abs_delta = subtype_count_max_abs_delta,
    neuron_count_delta = neuron_count_delta
  )
  slice_rows[[sl]] <- obs
  class_rows[[sl]] <- class_dt
  neuron_rows[[sl]] <- neuron_obs
}

slice_summary <- rbindlist(slice_rows, use.names = TRUE, fill = TRUE)
class_summary <- rbindlist(class_rows, use.names = TRUE, fill = TRUE)
null_audit <- rbindlist(null_rows, use.names = TRUE, fill = TRUE)
neuron_summary <- rbindlist(neuron_rows, use.names = TRUE, fill = TRUE)

slice_summary[, BH_FDR_same_neighbor := p.adjust(emp_p_same_neighbor_enriched, method = "BH")]
neuron_summary[, BH_FDR_neuron_segregated := p.adjust(emp_p_neuron_segregated, method = "BH")]

global_summary <- data.table(
  metric = c(
    "slices",
    "subtype_same_neighbor_FDR_pass",
    "subtype_same_neighbor_positive_vs_null_median",
    "neuron_control_FDR_pass",
    "neuron_control_positive_vs_null_median",
    "median_assortativity_index",
    "median_local_entropy"
  ),
  value = c(
    nrow(slice_summary),
    sum(slice_summary$BH_FDR_same_neighbor < 0.05, na.rm = TRUE),
    sum(slice_summary$observed_same_neighbor_fraction > slice_summary$null_median_same_neighbor_fraction, na.rm = TRUE),
    sum(neuron_summary$BH_FDR_neuron_segregated < 0.05, na.rm = TRUE),
    sum(neuron_summary$neuron_same_neighbor_fraction > neuron_summary$null_median_neuron_same_neighbor_fraction, na.rm = TRUE),
    median(slice_summary$assortativity_index, na.rm = TRUE),
    median(slice_summary$median_local_entropy, na.rm = TRUE)
  )
)

params <- data.table(
  parameter = c("input_weights", "k_nn", "n_perm", "subtype_mapping", "dominant_neuron_definition", "null"),
  value = c(
    weights_path,
    as.character(k_nn),
    as.character(n_perm),
    paste(names(subtype_map), subtype_map, sep = "=", collapse = "; "),
    "dominant RCTD label == Neurons across all RCTD columns",
    "toroidal coordinate shift of label field, nearest observed spot remap; no CSR/independent label permutation"
  )
)

fwrite(slice_summary, file.path(out_dir, "territory_segregation_perslice.csv"))
fwrite(class_summary, file.path(out_dir, "territory_segregation_by_subtype.csv"))
fwrite(neuron_summary, file.path(out_dir, "territory_segregation_neuron_control.csv"))
fwrite(null_audit, file.path(out_dir, "territory_segregation_toroidal_null_audit.csv"))
fwrite(global_summary, file.path(out_dir, "territory_segregation_global_summary.csv"))
fwrite(params, file.path(out_dir, "territory_segregation_parameters_source.csv"))

qs2::qs_save(
  list(
    slice_summary = slice_summary,
    class_summary = class_summary,
    neuron_summary = neuron_summary,
    global_summary = global_summary,
    parameters = params
  ),
  file.path(out_dir, "territory_segregation_results.qs2")
)

message("Wrote project 1 territory segregation outputs to: ", out_dir)
