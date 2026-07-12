#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-1: build spatial neighborhoods and RCTD-derived niche scores
# Purpose: construct k=6 and k=12 neighborhoods with distance constraints, then
#          compute neighbor-mean niche predictors excluding focal spots.
# STOP  : no association tests, no spatial null model, no enrichment conclusion.
# =============================================================================

suppressPackageStartupMessages({
  library(qs2)
  library(data.table)
})

base_dir <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
weights_path <- file.path(base_dir, "outputs/R9_A2_RCTD_weights_allslices_long.qs2")
out_dir <- file.path(base_dir, "outputs/R9_C3_neighborhood_niche")
tab_dir <- file.path(base_dir, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

k_main <- 6L
k_sens <- 12L

big <- as.data.table(qs2::qs_read(weights_path))
required <- c("spot_id", "slice", "x", "y",
              "Subtype3", "Subtype4",
              "Endothelial", "Mural cells",
              "Macrophages", "Microglial", "Monocytes",
              "Neurons")
miss <- setdiff(required, names(big))
if (length(miss)) stop("missing required columns: ", paste(miss, collapse = ", "))

big[, `MES-V` := Subtype3]
big[, `MES-I` := Subtype4]
big[, `MES-lineage` := Subtype3 + Subtype4]
big[, vascular_self := Endothelial + `Mural cells`]
big[, myeloid_self := Macrophages + Microglial + Monocytes]
big[, microglia_self := Microglial]
big[, macrophage_monocyte_self := Macrophages + Monocytes]
big[, neuron_self := Neurons]

response_cols <- c("MES-lineage", "MES-V", "MES-I")
niche_self_cols <- c(
  vascular = "vascular_self",
  myeloid = "myeloid_self",
  microglia = "microglia_self",
  macrophage_monocyte = "macrophage_monocyte_self",
  neuron_control = "neuron_self"
)

slice_dist_and_neighbors <- function(dt, k, threshold) {
  xy <- as.matrix(dt[, .(x, y)])
  n <- nrow(xy)
  dmat <- as.matrix(dist(xy))
  diag(dmat) <- Inf
  ord <- t(apply(dmat, 1, order))
  kth <- dmat[cbind(seq_len(n), ord[, k])]
  neigh <- vector("list", n)
  neigh_dist <- vector("list", n)
  for (i in seq_len(n)) {
    cand <- ord[i, seq_len(k)]
    keep <- dmat[i, cand] <= threshold
    neigh[[i]] <- cand[keep]
    neigh_dist[[i]] <- dmat[i, cand[keep]]
  }
  list(neigh = neigh, neigh_dist = neigh_dist, kth = kth)
}

make_neighbor_table <- function(dt, neigh_obj, k, threshold, threshold_mode) {
  rows <- lapply(seq_along(neigh_obj$neigh), function(i) {
    idx <- neigh_obj$neigh[[i]]
    if (!length(idx)) {
      return(data.table(
        slice = dt$slice[i],
        spot_id = dt$spot_id[i],
        neighbor_id = character(0),
        neighbor_rank = integer(0),
        distance = numeric(0)
      ))
    }
    data.table(
      slice = dt$slice[i],
      spot_id = dt$spot_id[i],
      neighbor_id = dt$spot_id[idx],
      neighbor_rank = seq_along(idx),
      distance = neigh_obj$neigh_dist[[i]]
    )
  })
  rbindlist(rows)[, `:=`(k = k, distance_threshold = threshold, threshold_mode = threshold_mode)]
}

make_niche_scores <- function(dt, neigh_obj, k, threshold, threshold_mode) {
  score <- data.table(
    spot_id = dt$spot_id,
    slice = dt$slice,
    image = if ("image" %in% names(dt)) dt$image else NA_character_,
    x = dt$x,
    y = dt$y
  )
  for (nm in response_cols) score[, (nm) := dt[[nm]]]
  score[, `:=`(k = k, distance_threshold = threshold, threshold_mode = threshold_mode)]
  score[, realized_neighbors := lengths(neigh_obj$neigh)]
  for (nm in names(niche_self_cols)) {
    col <- niche_self_cols[[nm]]
    vals <- vapply(neigh_obj$neigh, function(idx) {
      if (!length(idx)) return(NA_real_)
      mean(dt[[col]][idx], na.rm = TRUE)
    }, numeric(1))
    score[, paste0(nm, "_niche") := vals]
  }
  score
}

all_scores <- list()
all_neighbors <- list()
qc_rows <- list()
threshold_rows <- list()

for (sl in unique(big$slice)) {
  dt <- big[slice == sl]
  message("== slice ", sl, " n=", nrow(dt))

  ## k=6: per-slice first-ring threshold. C3-1 v1 showed that the global
  ## threshold failed for high-resolution x/y slices, so thresholds must be
  ## calibrated within each slice.
  nb6_probe <- slice_dist_and_neighbors(dt, k_main, Inf)
  threshold_k6 <- 1.5 * median(nb6_probe$kth, na.rm = TRUE)
  nb6 <- slice_dist_and_neighbors(dt, k_main, threshold_k6)
  all_neighbors[[paste(sl, "k6", sep = "__")]] <-
    make_neighbor_table(dt, nb6, k_main, threshold_k6, "per_slice_6thNN_median_1p5x")
  all_scores[[paste(sl, "k6", sep = "__")]] <-
    make_niche_scores(dt, nb6, k_main, threshold_k6, "per_slice_6thNN_median_1p5x")
  qc_rows[[paste(sl, "k6", sep = "__")]] <- data.table(
    slice = sl, k = k_main,
    threshold = threshold_k6,
    threshold_mode = "per_slice_6thNN_median_1p5x",
    n_spots = nrow(dt),
    kth_distance_median = median(nb6_probe$kth, na.rm = TRUE),
    kth_distance_q05 = as.numeric(quantile(nb6_probe$kth, 0.05, na.rm = TRUE)),
    kth_distance_q95 = as.numeric(quantile(nb6_probe$kth, 0.95, na.rm = TRUE)),
    realized_min = min(lengths(nb6$neigh)),
    realized_q05 = as.numeric(quantile(lengths(nb6$neigh), 0.05)),
    realized_median = median(lengths(nb6$neigh)),
    realized_q95 = as.numeric(quantile(lengths(nb6$neigh), 0.95)),
    realized_max = max(lengths(nb6$neigh)),
    frac_with_zero_neighbors = mean(lengths(nb6$neigh) == 0),
    frac_with_full_k = mean(lengths(nb6$neigh) == k_main)
  )

  ## k=12: per-slice second-ring threshold from 12th-neighbor distance.
  nb12_probe <- slice_dist_and_neighbors(dt, k_sens, Inf)
  threshold_k12 <- 1.5 * median(nb12_probe$kth, na.rm = TRUE)
  nb12 <- slice_dist_and_neighbors(dt, k_sens, threshold_k12)
  all_neighbors[[paste(sl, "k12", sep = "__")]] <-
    make_neighbor_table(dt, nb12, k_sens, threshold_k12, "per_slice_12thNN_median_1p5x")
  all_scores[[paste(sl, "k12", sep = "__")]] <-
    make_niche_scores(dt, nb12, k_sens, threshold_k12, "per_slice_12thNN_median_1p5x")
  qc_rows[[paste(sl, "k12", sep = "__")]] <- data.table(
    slice = sl, k = k_sens,
    threshold = threshold_k12,
    threshold_mode = "per_slice_12thNN_median_1p5x",
    n_spots = nrow(dt),
    kth_distance_median = median(nb12_probe$kth, na.rm = TRUE),
    kth_distance_q05 = as.numeric(quantile(nb12_probe$kth, 0.05, na.rm = TRUE)),
    kth_distance_q95 = as.numeric(quantile(nb12_probe$kth, 0.95, na.rm = TRUE)),
    realized_min = min(lengths(nb12$neigh)),
    realized_q05 = as.numeric(quantile(lengths(nb12$neigh), 0.05)),
    realized_median = median(lengths(nb12$neigh)),
    realized_q95 = as.numeric(quantile(lengths(nb12$neigh), 0.95)),
    realized_max = max(lengths(nb12$neigh)),
    frac_with_zero_neighbors = mean(lengths(nb12$neigh) == 0),
    frac_with_full_k = mean(lengths(nb12$neigh) == k_sens)
  )
}

score_dt <- rbindlist(all_scores)
neighbor_dt <- rbindlist(all_neighbors)
qc_dt <- rbindlist(qc_rows)

summary_dt <- qc_dt[, .(
  n_slices = .N,
  median_threshold = median(threshold),
  threshold_q05 = quantile(threshold, 0.05),
  threshold_q95 = quantile(threshold, 0.95),
  median_realized_neighbors = median(realized_median),
  q05_realized_neighbors = quantile(realized_median, 0.05),
  q95_realized_neighbors = quantile(realized_median, 0.95),
  max_frac_zero_neighbors = max(frac_with_zero_neighbors),
  min_frac_full_k = min(frac_with_full_k)
), by = .(k, threshold_mode)]

qs2::qs_save(
  list(
    niche_scores = score_dt,
    neighbor_edges = neighbor_dt,
    neighbor_qc = qc_dt,
    neighbor_summary = summary_dt,
    response_cols = response_cols,
    niche_self_cols = niche_self_cols,
    note = "C3-1 v2 only: per-slice calibrated neighborhoods and RCTD-derived niche score construction; no association test"
  ),
  file.path(out_dir, "R9_C3_1_neighborhood_niche_scores.qs2")
)

fwrite(score_dt, file.path(tab_dir, "R9_C3_1_neighborhood_niche_scores.csv"))
fwrite(qc_dt, file.path(tab_dir, "R9_C3_1_neighbor_count_qc_by_slice.csv"))
fwrite(summary_dt, file.path(tab_dir, "R9_C3_1_neighbor_count_summary.csv"))

cat("\n== neighbor summary:\n")
print(summary_dt)
cat("\n== first rows of per-slice QC:\n")
print(head(qc_dt, 12))
cat("\n[STOP C3-1] Confirm realized neighbor counts before association testing / spatial null model.\n")
