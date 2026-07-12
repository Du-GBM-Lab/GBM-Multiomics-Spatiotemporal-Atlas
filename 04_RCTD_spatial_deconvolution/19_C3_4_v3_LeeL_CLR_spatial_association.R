#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-4 v3 method A: Lee's L spatial association on CLR RCTD weights
# Purpose:
#   Use a standard bivariate spatial association statistic instead of hand-made
#   subset/Jaccard tests. This script performs:
#     1) per-slice W audit for k=6/k=12 distance-constrained neighbor graphs;
#     2) Lee's L for focal CLR MES vs focal CLR niche components;
#     3) conditional permutation null by shuffling the niche labels within slice.
#
# STOP:
#   Run Lee's L only. No Ripley/cross-K, no Jaccard, no figures.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
  library(FNN)
  library(Matrix)
})

set.seed(1)

base_dir <- getwd()
score_path <- file.path(base_dir, "tables/C3_niche_preflight/R9_C3_1_neighborhood_niche_scores.csv")
weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
out_dir <- file.path(base_dir, "tables/C3_4_local_niche")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

responses <- c("MES-lineage", "MES-V", "MES-I")
predictors <- c("vascular_niche", "myeloid_niche", "neuron_control_niche")
k_values <- c(6L, 12L)
n_perm <- 999L
eps <- 1e-6

scores <- fread(score_path)
weights <- as.data.table(qs2::qs_read(weights_path))
stopifnot(all(c("spot_id", "slice", "x", "y", "k", "distance_threshold") %in% names(scores)))

meta_cols <- c("spot_id", "slice", "image", "x", "y")
celltype_cols <- setdiff(names(weights), meta_cols)
stopifnot(all(c("Subtype3", "Subtype4", "Endothelial", "Mural cells",
                "Macrophages", "Microglial", "Monocytes", "Neurons") %in% celltype_cols))

make_clr <- function(wdt) {
  mat <- as.matrix(wdt[, ..celltype_cols])
  gm <- rowMeans(log(mat + eps))
  data.table(
    spot_id = wdt$spot_id,
    slice = wdt$slice,
    x = wdt$x,
    y = wdt$y,
    `MES-lineage` = log(wdt$Subtype3 + wdt$Subtype4 + eps) - gm,
    `MES-V` = log(wdt$Subtype3 + eps) - gm,
    `MES-I` = log(wdt$Subtype4 + eps) - gm,
    vascular_niche = log(wdt$Endothelial + wdt[["Mural cells"]] + eps) - gm,
    myeloid_niche = log(wdt$Macrophages + wdt$Microglial + wdt$Monocytes + eps) - gm,
    neuron_control_niche = log(wdt$Neurons + eps) - gm
  )
}

make_W <- function(dt) {
  xy <- as.matrix(dt[, .(x, y)])
  threshold <- median(dt$distance_threshold, na.rm = TRUE)
  nn <- FNN::get.knnx(xy, xy, k = min(32L, nrow(dt)))
  i_vec <- integer(0)
  j_vec <- integer(0)
  for (ii in seq_len(nrow(dt))) {
    idx <- nn$nn.index[ii, ]
    dst <- nn$nn.dist[ii, ]
    keep <- idx != ii & dst <= threshold
    if (any(keep)) {
      i_vec <- c(i_vec, rep(ii, sum(keep)))
      j_vec <- c(j_vec, idx[keep])
    }
  }
  if (!length(i_vec)) stop("No neighbors recovered for slice=", dt$slice[1], " k=", dt$k[1])
  A <- sparseMatrix(i = i_vec, j = j_vec, x = 1, dims = c(nrow(dt), nrow(dt)))
  rs <- Matrix::rowSums(A)
  W <- Diagonal(x = ifelse(rs > 0, 1 / rs, 0)) %*% A
  list(W = W, threshold = threshold, neighbor_count = as.numeric(rs))
}

lee_L <- function(x, y, W) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 50) return(NA_real_)
  x <- x[ok]
  y <- y[ok]
  W2 <- W[ok, ok, drop = FALSE]
  rs <- Matrix::rowSums(W2)
  keep <- rs > 0
  if (sum(keep) < 50) return(NA_real_)
  x <- x[keep]
  y <- y[keep]
  W2 <- W2[keep, keep, drop = FALSE]
  rs <- Matrix::rowSums(W2)
  W2 <- Diagonal(x = 1 / rs) %*% W2

  xc <- x - mean(x)
  yc <- y - mean(y)
  den <- sqrt(sum(xc^2)) * sqrt(sum(yc^2))
  if (!is.finite(den) || den == 0) return(NA_real_)
  S2 <- sum(Matrix::rowSums(W2)^2)
  as.numeric((length(x) / S2) * sum((W2 %*% xc) * (W2 %*% yc)) / den)
}

lee_with_perm <- function(x, y, W, nperm = n_perm) {
  obs <- lee_L(x, y, W)
  nul <- replicate(nperm, lee_L(x, sample(y), W))
  nul <- nul[is.finite(nul)]
  data.table(
    observed_L = obs,
    n_perm = length(nul),
    null_median = median(nul, na.rm = TRUE),
    null_q25 = as.numeric(quantile(nul, 0.25, na.rm = TRUE)),
    null_q75 = as.numeric(quantile(nul, 0.75, na.rm = TRUE)),
    p_emp_two_sided = if (is.finite(obs) && length(nul)) (1 + sum(abs(nul) >= abs(obs))) / (length(nul) + 1) else NA_real_,
    p_emp_greater = if (is.finite(obs) && length(nul)) (1 + sum(nul >= obs)) / (length(nul) + 1) else NA_real_,
    p_emp_less = if (is.finite(obs) && length(nul)) (1 + sum(nul <= obs)) / (length(nul) + 1) else NA_real_
  )
}

sign_p <- function(x, alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  x <- x[is.finite(x) & x != 0]
  if (!length(x)) return(NA_real_)
  binom.test(sum(x > 0), length(x), p = 0.5, alternative = alternative)$p.value
}

wilcox_p <- function(x, alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, alternative = alternative, exact = FALSE)$p.value)
}

summarize_lee <- function(dt) {
  dt[, .(
    n_slices = sum(is.finite(observed_L)),
    median_L = median(observed_L, na.rm = TRUE),
    q25_L = as.numeric(quantile(observed_L, 0.25, na.rm = TRUE)),
    q75_L = as.numeric(quantile(observed_L, 0.75, na.rm = TRUE)),
    n_positive = sum(observed_L > 0, na.rm = TRUE),
    n_negative = sum(observed_L < 0, na.rm = TRUE),
    sign_p_greater = sign_p(observed_L, "greater"),
    wilcox_p_greater = wilcox_p(observed_L, "greater"),
    median_emp_p_two_sided = median(p_emp_two_sided, na.rm = TRUE),
    n_slice_emp_p_lt_0p05 = sum(p_emp_two_sided < 0.05, na.rm = TRUE),
    median_perm = median(n_perm, na.rm = TRUE),
    min_perm = min(n_perm, na.rm = TRUE)
  ), by = .(k, response, predictor)]
}

message("== Building CLR features from full 21-class RCTD composition.")
clr <- make_clr(weights)
perslice <- list()
audit <- list()
ii <- aa <- 0L

for (sl in sort(unique(scores$slice))) {
  for (kk in k_values) {
    coord <- scores[slice == sl & k == kk,
                    .(spot_id, slice, x, y, k, distance_threshold)]
    dt <- merge(coord, clr[slice == sl], by = c("spot_id", "slice", "x", "y"),
                all.x = TRUE, sort = FALSE)
    wg <- make_W(dt)
    aa <- aa + 1L
    audit[[aa]] <- data.table(
      slice = sl,
      k = kk,
      n_spots = nrow(dt),
      distance_threshold = wg$threshold,
      min_neighbors = min(wg$neighbor_count),
      median_neighbors = median(wg$neighbor_count),
      max_neighbors = max(wg$neighbor_count),
      zero_neighbor_fraction = mean(wg$neighbor_count == 0),
      row_sum_median = median(Matrix::rowSums(wg$W)),
      row_sum_min = min(Matrix::rowSums(wg$W))
    )

    for (resp in responses) {
      for (pred in predictors) {
        ii <- ii + 1L
        ans <- lee_with_perm(dt[[resp]], dt[[pred]], wg$W)
        ans[, `:=`(slice = sl, k = kk, response = resp, predictor = pred)]
        perslice[[ii]] <- ans
      }
    }
    message("  done Lee's L slice=", sl, " k=", kk)
  }
}

perslice <- rbindlist(perslice, use.names = TRUE, fill = TRUE)
setcolorder(perslice, c("slice", "k", "response", "predictor",
                        setdiff(names(perslice), c("slice", "k", "response", "predictor"))))
audit <- rbindlist(audit)
summary <- summarize_lee(perslice)
summary_no265 <- summarize_lee(perslice[slice != "#UKF265_T_ST"])
summary[, sensitivity := "all_slices"]
summary_no265[, sensitivity := "exclude_UKF265"]
summary_all <- rbindlist(list(summary, summary_no265), use.names = TRUE)

fwrite(audit, file.path(out_dir, "R9_C3_4_v3_LeeL_W_audit.csv"))
fwrite(perslice, file.path(out_dir, "R9_C3_4_v3_LeeL_perslice.csv"))
fwrite(summary_all, file.path(out_dir, "R9_C3_4_v3_LeeL_summary.csv"))

cat("\n== W audit summary:\n")
print(audit[, .(
  n_slice_k = .N,
  min_median_neighbors = min(median_neighbors),
  median_median_neighbors = median(median_neighbors),
  max_zero_neighbor_fraction = max(zero_neighbor_fraction),
  min_row_sum_min = min(row_sum_min)
)])

cat("\n== Lee's L summary, k=6 focus pairs:\n")
print(summary_all[k == 6 & sensitivity == "all_slices" &
                    response %in% responses &
                    predictor %in% c("vascular_niche", "myeloid_niche", "neuron_control_niche")]
      [order(response, predictor)])

cat("\n== Lee's L summary, k=6 focus pairs, excluding UKF265:\n")
print(summary_all[k == 6 & sensitivity == "exclude_UKF265" &
                    response %in% responses &
                    predictor %in% c("vascular_niche", "myeloid_niche", "neuron_control_niche")]
      [order(response, predictor)])

cat("\n[STOP C3-4 v3 Lee's L] Review W audit, vascular effect, and neuron_control.",
    " Do not proceed to Ripley's K until this STOP is checked.\n")
