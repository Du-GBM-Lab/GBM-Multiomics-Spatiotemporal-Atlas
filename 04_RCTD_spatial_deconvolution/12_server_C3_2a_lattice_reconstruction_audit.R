#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-2a: lattice reconstruction audit from x/y coordinates
# Purpose: reconstruct slice-wise lattice indices from physical x/y coordinates
#          and test whether field-shift/toroidal-style spatial null is feasible.
# STOP  : no association tests, no permutation p-values, no biological conclusion.
# =============================================================================

suppressPackageStartupMessages({
  library(qs2)
  library(data.table)
})

base_dir <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
weights_path <- file.path(base_dir, "outputs/R9_A2_RCTD_weights_allslices_long.qs2")
out_dir <- file.path(base_dir, "outputs/R9_C3_lattice_audit")
tab_dir <- file.path(base_dir, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

big <- as.data.table(qs2::qs_read(weights_path))
stopifnot(all(c("spot_id", "slice", "x", "y") %in% names(big)))

## Cluster one coordinate axis into lattice rows/columns by grouping adjacent
## sorted unique values separated by <= tolerance.
cluster_axis <- function(v) {
  u <- sort(unique(as.numeric(v)))
  if (length(u) < 3) return(list(idx = match(v, u), centers = u, tol = NA_real_))
  dif <- diff(u)
  pos_dif <- dif[dif > 0]
  ## For exact/integer image coordinates, repeated row/column values produce
  ## large gaps between true grid lines. For jittered coordinates, use a small
  ## tolerance relative to the smallest positive step.
  tol <- max(1e-6, quantile(pos_dif, 0.10, na.rm = TRUE) * 0.25)
  group <- integer(length(u))
  group[1] <- 1L
  for (i in 2:length(u)) {
    group[i] <- group[i - 1] + as.integer((u[i] - u[i - 1]) > tol)
  }
  centers <- as.numeric(tapply(u, group, median))
  idx <- group[match(as.numeric(v), u)]
  list(idx = idx, centers = centers, tol = tol)
}

nearest_distance <- function(dt) {
  xy <- as.matrix(dt[, .(x, y)])
  dmat <- as.matrix(dist(xy))
  diag(dmat) <- Inf
  apply(dmat, 1, min)
}

audit_one_slice <- function(dt) {
  ## Try both orientations and keep the one where row grouping is more compact:
  ## orientation A: y -> row, x -> col; orientation B: x -> row, y -> col.
  ca_y <- cluster_axis(dt$y)
  ca_x <- cluster_axis(dt$x)

  cand <- list(
    y_as_row = list(row = ca_y$idx, col = ca_x$idx,
                    n_row = length(ca_y$centers), n_col = length(ca_x$centers),
                    row_tol = ca_y$tol, col_tol = ca_x$tol),
    x_as_row = list(row = ca_x$idx, col = ca_y$idx,
                    n_row = length(ca_x$centers), n_col = length(ca_y$centers),
                    row_tol = ca_x$tol, col_tol = ca_y$tol)
  )

  score_cand <- rbindlist(lapply(names(cand), function(nm) {
    z <- cand[[nm]]
    key <- paste(z$row, z$col, sep = "_")
    data.table(
      orientation = nm,
      n_row = z$n_row,
      n_col = z$n_col,
      n_duplicate_cells = length(key) - uniqueN(key),
      occupancy = uniqueN(key) / (z$n_row * z$n_col),
      median_spots_per_row = median(tabulate(z$row, nbins = z$n_row)),
      median_spots_per_col = median(tabulate(z$col, nbins = z$n_col)),
      row_tol = z$row_tol,
      col_tol = z$col_tol
    )
  }))

  ## Prefer no duplicate lattice cells and the orientation with higher occupancy;
  ## if tied, use y_as_row as conventional image-row orientation.
  score_cand[, duplicate_penalty := n_duplicate_cells > 0]
  best_name <- score_cand[order(duplicate_penalty, -occupancy, orientation)][1, orientation]
  best <- cand[[best_name]]

  lat <- copy(dt[, .(spot_id, slice, x, y)])
  lat[, `:=`(lattice_row = best$row, lattice_col = best$col,
             lattice_orientation = best_name)]

  nd <- nearest_distance(dt)

  ## Coverage for candidate field shifts. Shift each occupied lattice coordinate
  ## by delta row/col with wrap over the reconstructed bounding box; retain only
  ## target coordinates corresponding to real spots.
  occupied <- unique(lat[, .(lattice_row, lattice_col)])
  key0 <- paste(occupied$lattice_row, occupied$lattice_col, sep = "_")
  key_set <- new.env(hash = TRUE, parent = emptyenv())
  for (k in key0) assign(k, TRUE, envir = key_set)

  shifts <- CJ(delta_row = -3:3, delta_col = -3:3)
  shifts <- shifts[!(delta_row == 0 & delta_col == 0)]
  nrow_box <- best$n_row
  ncol_box <- best$n_col
  cov <- shifts[, {
    tr <- ((occupied$lattice_row + delta_row - 1L) %% nrow_box) + 1L
    tc <- ((occupied$lattice_col + delta_col - 1L) %% ncol_box) + 1L
    target_keys <- paste(tr, tc, sep = "_")
    .(coverage = mean(vapply(target_keys, exists, logical(1), envir = key_set, inherits = FALSE)))
  }, by = .(delta_row, delta_col)]

  list(
    lattice = lat,
    orientation_audit = score_cand,
    summary = data.table(
      slice = dt$slice[1],
      n_spots = nrow(dt),
      chosen_orientation = best_name,
      n_lattice_rows = best$n_row,
      n_lattice_cols = best$n_col,
      n_duplicate_lattice_cells = score_cand[orientation == best_name, n_duplicate_cells],
      lattice_occupancy = score_cand[orientation == best_name, occupancy],
      nn_dist_median = median(nd),
      nn_dist_q05 = as.numeric(quantile(nd, 0.05)),
      nn_dist_q95 = as.numeric(quantile(nd, 0.95)),
      shift_coverage_median = median(cov$coverage),
      shift_coverage_q05 = as.numeric(quantile(cov$coverage, 0.05)),
      shift_coverage_min = min(cov$coverage),
      n_shift_ge_0p70 = sum(cov$coverage >= 0.70),
      n_shift_ge_0p80 = sum(cov$coverage >= 0.80),
      n_shift_ge_0p90 = sum(cov$coverage >= 0.90)
    ),
    shift_coverage = cov[, slice := dt$slice[1]][]
  )
}

res <- lapply(split(big, big$slice), audit_one_slice)
lattice_dt <- rbindlist(lapply(res, `[[`, "lattice"))
orientation_audit <- rbindlist(lapply(res, `[[`, "orientation_audit"), idcol = "slice_name")
summary_dt <- rbindlist(lapply(res, `[[`, "summary"))
shift_dt <- rbindlist(lapply(res, `[[`, "shift_coverage"))

global_summary <- summary_dt[, .(
  n_slices = .N,
  n_slices_with_duplicates = sum(n_duplicate_lattice_cells > 0),
  median_occupancy = median(lattice_occupancy),
  q05_occupancy = quantile(lattice_occupancy, 0.05),
  median_shift_coverage = median(shift_coverage_median),
  min_shift_coverage_median = min(shift_coverage_median),
  min_n_shift_ge_0p70 = min(n_shift_ge_0p70),
  min_n_shift_ge_0p80 = min(n_shift_ge_0p80),
  min_n_shift_ge_0p90 = min(n_shift_ge_0p90)
)]

qs2::qs_save(
  list(
    lattice = lattice_dt,
    slice_summary = summary_dt,
    shift_coverage = shift_dt,
    orientation_audit = orientation_audit,
    global_summary = global_summary,
    note = "C3-2a lattice reconstruction audit only; no spatial null p-values"
  ),
  file.path(out_dir, "R9_C3_2a_lattice_reconstruction_audit.qs2")
)

fwrite(lattice_dt, file.path(tab_dir, "R9_C3_2a_lattice_index_by_spot.csv"))
fwrite(summary_dt, file.path(tab_dir, "R9_C3_2a_lattice_slice_summary.csv"))
fwrite(shift_dt, file.path(tab_dir, "R9_C3_2a_shift_coverage_by_slice.csv"))
fwrite(orientation_audit, file.path(tab_dir, "R9_C3_2a_lattice_orientation_audit.csv"))
fwrite(global_summary, file.path(tab_dir, "R9_C3_2a_lattice_global_summary.csv"))

cat("\n== lattice global summary:\n")
print(global_summary)
cat("\n== lattice slice summary:\n")
print(summary_dt)
cat("\n[STOP C3-2a] Confirm lattice reconstruction and shift coverage before C3-2 association/null model.\n")
