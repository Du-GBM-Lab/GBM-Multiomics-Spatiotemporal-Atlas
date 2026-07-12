#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-2a v2: physical vector-shift audit from nearest-neighbor geometry
# Purpose: x/y axes are not aligned to array rows/cols. Infer hex-like basis
#          vectors from local nearest-neighbor vectors and test whether spatial
#          field shifts can be matched back to real spots with sufficient coverage.
# STOP  : no association tests, no p-values, no biological conclusion.
# =============================================================================

suppressPackageStartupMessages({
  library(qs2)
  library(data.table)
})

set.seed(1)

base_dir <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
weights_path <- file.path(base_dir, "outputs/R9_A2_RCTD_weights_allslices_long.qs2")
out_dir <- file.path(base_dir, "outputs/R9_C3_lattice_audit")
tab_dir <- file.path(base_dir, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

big <- as.data.table(qs2::qs_read(weights_path))
stopifnot(all(c("spot_id", "slice", "x", "y") %in% names(big)))

nearest_match <- function(query, ref, max_dist) {
  if (requireNamespace("FNN", quietly = TRUE)) {
    nn <- FNN::get.knnx(ref, query, k = 1)
    return(list(index = nn$nn.index[, 1], dist = nn$nn.dist[, 1],
                valid = nn$nn.dist[, 1] <= max_dist))
  }
  ## Fallback: chunked brute force.
  idx <- integer(nrow(query))
  dst <- numeric(nrow(query))
  chunk <- 500L
  for (i0 in seq(1L, nrow(query), by = chunk)) {
    i1 <- min(nrow(query), i0 + chunk - 1L)
    q <- query[i0:i1, , drop = FALSE]
    d <- sapply(seq_len(nrow(q)), function(i) {
      (ref[, 1] - q[i, 1])^2 + (ref[, 2] - q[i, 2])^2
    })
    if (is.null(dim(d))) d <- matrix(d, ncol = 1)
    mi <- apply(d, 2, which.min)
    md <- sqrt(d[cbind(mi, seq_along(mi))])
    idx[i0:i1] <- mi
    dst[i0:i1] <- md
  }
  list(index = idx, dist = dst, valid = dst <= max_dist)
}

infer_axes <- function(xy) {
  dmat <- as.matrix(dist(xy))
  diag(dmat) <- Inf
  nn1 <- apply(dmat, 1, min)
  pitch <- median(nn1, na.rm = TRUE)

  ord <- t(apply(dmat, 1, order))
  k <- min(6L, nrow(xy) - 1L)
  pairs <- rbindlist(lapply(seq_len(nrow(xy)), function(i) {
    j <- ord[i, seq_len(k)]
    dd <- dmat[i, j]
    keep <- dd <= 1.25 * pitch
    if (!any(keep)) return(NULL)
    data.table(dx = xy[j[keep], 1] - xy[i, 1],
               dy = xy[j[keep], 2] - xy[i, 2],
               dist = dd[keep])
  }))
  if (nrow(pairs) < 30) stop("Too few nearest-neighbor vectors for axis inference.")

  ang <- atan2(pairs$dy, pairs$dx)
  ang[ang < 0] <- ang[ang < 0] + pi
  feat <- cbind(cos(2 * ang), sin(2 * ang))
  km <- kmeans(feat, centers = 3, nstart = 25)

  axes <- rbindlist(lapply(sort(unique(km$cluster)), function(cl) {
    v <- pairs[km$cluster == cl]
    ctr <- colMeans(feat[km$cluster == cl, , drop = FALSE])
    theta <- atan2(ctr[2], ctr[1]) / 2
    if (theta < 0) theta <- theta + pi
    u <- c(cos(theta), sin(theta))
    vv <- as.matrix(v[, .(dx, dy)])
    flip <- (vv[, 1] * u[1] + vv[, 2] * u[2]) < 0
    vv[flip, ] <- -vv[flip, ]
    data.table(
      axis = cl,
      theta = theta,
      dx = median(vv[, 1]),
      dy = median(vv[, 2]),
      length = median(sqrt(rowSums(vv^2))),
      n_vectors = nrow(v)
    )
  }))
  axes[, theta_deg := theta * 180 / pi]
  axes[order(theta)]
}

audit_slice <- function(dt) {
  xy <- as.matrix(dt[, .(x, y)])
  axes <- infer_axes(xy)
  axes[, slice := dt$slice[1]]

  ## Choose two axes closest to 60 degrees apart as basis.
  comb <- CJ(a = seq_len(nrow(axes)), b = seq_len(nrow(axes)))[a < b]
  comb[, angle_diff := abs(axes$theta[a] - axes$theta[b])]
  comb[, angle_diff := pmin(angle_diff, pi - angle_diff)]
  comb[, target_diff := abs(angle_diff - pi / 3)]
  best <- comb[order(target_diff)][1]
  b1 <- as.numeric(axes[best$a, .(dx, dy)])
  b2 <- as.numeric(axes[best$b, .(dx, dy)])

  shifts <- CJ(delta_a = -3:3, delta_b = -3:3)
  shifts <- shifts[!(delta_a == 0 & delta_b == 0)]
  shifts[, `:=`(
    shift_x = delta_a * b1[1] + delta_b * b2[1],
    shift_y = delta_a * b1[2] + delta_b * b2[2]
  )]
  shifts[, shift_dist := sqrt(shift_x^2 + shift_y^2)]
  ## Exclude tiny shifts and extremely large shifts for null feasibility audit.
  shifts <- shifts[shift_dist >= median(axes$length) & shift_dist <= 4 * median(axes$length)]

  max_match_dist <- 0.35 * median(axes$length)
  cov <- shifts[, {
    query <- cbind(xy[, 1] + shift_x, xy[, 2] + shift_y)
    nn <- nearest_match(query, xy, max_match_dist)
    .(
      coverage = mean(nn$valid),
      median_match_dist = median(nn$dist[nn$valid], na.rm = TRUE),
      max_match_dist = max_match_dist
    )
  }, by = .(delta_a, delta_b, shift_x, shift_y, shift_dist)]
  cov[, slice := dt$slice[1]]

  summary <- data.table(
    slice = dt$slice[1],
    n_spots = nrow(dt),
    pitch = median(axes$length),
    basis_axis_a = axes$axis[best$a],
    basis_axis_b = axes$axis[best$b],
    basis_angle_diff_deg = best$angle_diff * 180 / pi,
    max_match_dist = max_match_dist,
    n_candidate_shifts = nrow(cov),
    median_shift_coverage = median(cov$coverage),
    q25_shift_coverage = quantile(cov$coverage, 0.25),
    q75_shift_coverage = quantile(cov$coverage, 0.75),
    max_shift_coverage = max(cov$coverage),
    n_shift_ge_0p50 = sum(cov$coverage >= 0.50),
    n_shift_ge_0p70 = sum(cov$coverage >= 0.70),
    n_shift_ge_0p80 = sum(cov$coverage >= 0.80)
  )

  list(axes = axes, shift_coverage = cov, summary = summary)
}

res <- lapply(split(big, big$slice), audit_slice)
axes_dt <- rbindlist(lapply(res, `[[`, "axes"))
shift_dt <- rbindlist(lapply(res, `[[`, "shift_coverage"))
summary_dt <- rbindlist(lapply(res, `[[`, "summary"))
global_summary <- summary_dt[, .(
  n_slices = .N,
  median_pitch = median(pitch),
  median_shift_coverage = median(median_shift_coverage),
  min_median_shift_coverage = min(median_shift_coverage),
  median_max_shift_coverage = median(max_shift_coverage),
  min_n_shift_ge_0p50 = min(n_shift_ge_0p50),
  min_n_shift_ge_0p70 = min(n_shift_ge_0p70),
  min_n_shift_ge_0p80 = min(n_shift_ge_0p80)
)]

qs2::qs_save(
  list(
    axes = axes_dt,
    shift_coverage = shift_dt,
    slice_summary = summary_dt,
    global_summary = global_summary,
    note = "C3-2a v2 vector-shift feasibility audit only; no p-values"
  ),
  file.path(out_dir, "R9_C3_2a2_vector_shift_audit.qs2")
)

fwrite(axes_dt, file.path(tab_dir, "R9_C3_2a2_inferred_hex_axes.csv"))
fwrite(shift_dt, file.path(tab_dir, "R9_C3_2a2_vector_shift_coverage.csv"))
fwrite(summary_dt, file.path(tab_dir, "R9_C3_2a2_vector_shift_slice_summary.csv"))
fwrite(global_summary, file.path(tab_dir, "R9_C3_2a2_vector_shift_global_summary.csv"))

cat("\n== vector-shift global summary:\n")
print(global_summary)
cat("\n== vector-shift slice summary:\n")
print(summary_dt)
cat("\n[STOP C3-2a2] Decide whether vector-shift coverage is sufficient or switch to spatial block permutation.\n")
