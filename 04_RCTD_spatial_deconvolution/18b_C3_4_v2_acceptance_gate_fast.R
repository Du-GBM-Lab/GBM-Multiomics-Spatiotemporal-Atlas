#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-4 v2 fast acceptance gate
# Minimal gate check after the full v2 script proved too heavy:
#   1) CLR observed slice-level associations: does neuron_control return near 0?
#   2) Lattice-toroidal hotspot-count audit for k=6 MES-lineage focus pairs:
#      can shifted/original hotspot count pass median >=0.95 and no slice <0.9?
#
# This is an audit-only STOP script: no Jaccard, no adjacency, no figures.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
  library(FNN)
})

set.seed(1)

base_dir <- getwd()
score_path <- file.path(base_dir, "tables/C3_niche_preflight/R9_C3_1_neighborhood_niche_scores.csv")
weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
shift_path <- file.path(base_dir, "tables/C3_niche_preflight/R9_C3_2a2_vector_shift_coverage.csv")
out_dir <- file.path(base_dir, "tables/C3_4_local_niche")

responses <- c("MES-lineage", "MES-V", "MES-I")
predictors <- c("vascular_niche", "myeloid_niche", "neuron_control_niche")
focus_predictors <- c("vascular_niche", "myeloid_niche", "neuron_control_niche")
hotspot_cut <- 0.10
eps <- 1e-6

scores <- fread(score_path)
scores6 <- scores[k == 6]
weights <- as.data.table(qs2::qs_read(weights_path))
old_shifts <- fread(shift_path)

meta_cols <- c("spot_id", "slice", "image", "x", "y")
celltype_cols <- setdiff(names(weights), meta_cols)

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 50) return(NA_real_)
  x <- x[ok]; y <- y[ok]
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

make_clr_features <- function(wdt) {
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

neighbor_mean <- function(dt, feature) {
  xy <- as.matrix(dt[, .(x, y)])
  threshold <- median(dt$distance_threshold, na.rm = TRUE)
  nn <- FNN::get.knnx(xy, xy, k = min(nrow(dt), 16L))
  out <- rep(NA_real_, nrow(dt))
  for (i in seq_len(nrow(dt))) {
    idx <- nn$nn.index[i, ]
    dst <- nn$nn.dist[i, ]
    keep <- idx != i & dst <= threshold
    if (any(keep)) out[i] <- mean(feature[idx[keep]], na.rm = TRUE)
  }
  out
}

fit_basis_shifts <- function(sl, max_delta = 18L) {
  sh0 <- old_shifts[slice == sl]
  fitx <- lm(shift_x ~ 0 + delta_a + delta_b, data = sh0)
  fity <- lm(shift_y ~ 0 + delta_a + delta_b, data = sh0)
  b1 <- c(coef(fitx)[["delta_a"]], coef(fity)[["delta_a"]])
  b2 <- c(coef(fitx)[["delta_b"]], coef(fity)[["delta_b"]])
  grid <- CJ(delta_a = -max_delta:max_delta, delta_b = -max_delta:max_delta)
  grid <- grid[!(delta_a == 0 & delta_b == 0)]
  grid[, `:=`(
    shift_x = delta_a * b1[1] + delta_b * b2[1],
    shift_y = delta_a * b1[2] + delta_b * b2[2],
    b1x = b1[1], b1y = b1[2], b2x = b2[1], b2y = b2[2]
  )]
  grid
}

infer_lattice <- function(dt, shifts) {
  b1 <- c(shifts$b1x[1], shifts$b1y[1])
  b2 <- c(shifts$b2x[1], shifts$b2y[1])
  B <- cbind(b1, b2)
  xy <- as.matrix(dt[, .(x, y)])
  origin <- xy[1, ]
  uv <- t(solve(B, t(xy - matrix(origin, nrow(xy), 2, byrow = TRUE))))
  ij <- round(uv)
  recon <- t(B %*% t(ij)) + matrix(origin, nrow(ij), 2, byrow = TRUE)
  data.table(row_id = seq_len(nrow(dt)),
             i = ij[, 1], j = ij[, 2],
             recon_error = sqrt(rowSums((xy - recon)^2)))
}

toroidal_ratio <- function(dt, shifts, predictor) {
  lat <- infer_lattice(dt, shifts)
  dup <- duplicated(lat[, .(i, j)]) | duplicated(lat[, .(i, j)], fromLast = TRUE)
  lat_u <- lat[!dup]
  map_row <- setNames(lat_u$row_id, paste(lat_u$i, lat_u$j))
  ir <- range(lat_u$i); jr <- range(lat_u$j)
  ni <- diff(ir) + 1L; nj <- diff(jr) + 1L
  x <- dt[[predictor]]
  hot <- is.finite(x) & x >= quantile(x, 1 - hotspot_cut, na.rm = TRUE)
  orig <- sum(hot)
  ratios <- numeric(nrow(shifts))
  for (s in seq_len(nrow(shifts))) {
    ti <- ((lat_u$i + shifts$delta_a[s] - ir[1]) %% ni) + ir[1]
    tj <- ((lat_u$j + shifts$delta_b[s] - jr[1]) %% nj) + jr[1]
    target <- map_row[paste(ti, tj)]
    ok <- !is.na(target)
    shifted <- rep(FALSE, nrow(dt))
    shifted[lat_u$row_id[ok]] <- hot[as.integer(target[ok])]
    ratios[s] <- sum(shifted) / orig
  }
  data.table(
    predictor = predictor,
    n_spots = nrow(dt),
    observed_hotspots = orig,
    n_lattice_unique = nrow(lat_u),
    n_lattice_duplicate_spots = sum(dup),
    median_recon_error = median(lat$recon_error),
    q95_recon_error = as.numeric(quantile(lat$recon_error, 0.95)),
    n_shifts = length(ratios),
    count_ratio_median = median(ratios, na.rm = TRUE),
    count_ratio_min = min(ratios, na.rm = TRUE),
    count_ratio_q25 = as.numeric(quantile(ratios, 0.25, na.rm = TRUE)),
    count_ratio_q75 = as.numeric(quantile(ratios, 0.75, na.rm = TRUE))
  )
}

clr <- make_clr_features(weights)
clr_perslice <- list()
tor_perslice <- list()
i1 <- i2 <- 0L

for (sl in sort(unique(scores6$slice))) {
  dt0 <- scores6[slice == sl]
  dtc <- merge(dt0[, .(spot_id, slice, x, y, distance_threshold)],
               clr[slice == sl],
               by = c("spot_id", "slice", "x", "y"),
               all.x = TRUE, sort = FALSE)
  for (pred in predictors) {
    dtc[[pred]] <- neighbor_mean(dtc, dtc[[pred]])
  }
  for (resp in responses) {
    for (pred in predictors) {
      i1 <- i1 + 1L
      clr_perslice[[i1]] <- data.table(
        slice = sl, k = 6L, response = resp, predictor = pred,
        observed_rho = safe_cor(dtc[[resp]], dtc[[pred]])
      )
    }
  }

  shifts <- fit_basis_shifts(sl)
  for (pred in focus_predictors) {
    i2 <- i2 + 1L
    tmp <- toroidal_ratio(dt0, shifts, pred)
    tmp[, `:=`(slice = sl, k = 6L, response = "MES-lineage")]
    tor_perslice[[i2]] <- tmp
  }
  message("done ", sl)
}

clr_perslice <- rbindlist(clr_perslice)
clr_summary <- clr_perslice[, .(
  n_slices = sum(is.finite(observed_rho)),
  median_rho = median(observed_rho, na.rm = TRUE),
  q25_rho = as.numeric(quantile(observed_rho, 0.25, na.rm = TRUE)),
  q75_rho = as.numeric(quantile(observed_rho, 0.75, na.rm = TRUE)),
  n_positive = sum(observed_rho > 0, na.rm = TRUE),
  n_negative = sum(observed_rho < 0, na.rm = TRUE)
), by = .(k, response, predictor)]
clr_summary[, method := "clr_full_composition_observed_only"]

tor_perslice <- rbindlist(tor_perslice)
tor_perslice[, pass_slice := count_ratio_median >= 0.95 & count_ratio_min >= 0.90]
tor_summary <- tor_perslice[, .(
  n_slices = .N,
  global_median_count_ratio = median(count_ratio_median, na.rm = TRUE),
  min_slice_count_ratio_median = min(count_ratio_median, na.rm = TRUE),
  min_slice_count_ratio_min = min(count_ratio_min, na.rm = TRUE),
  n_slice_pass = sum(pass_slice),
  max_duplicate_spots = max(n_lattice_duplicate_spots),
  median_q95_recon_error = median(q95_recon_error)
), by = .(k, response, predictor)]

fwrite(clr_perslice, file.path(out_dir, "R9_C3_4_v2_fast_CLR_observed_perslice.csv"))
fwrite(clr_summary, file.path(out_dir, "R9_C3_4_v2_fast_CLR_observed_summary.csv"))
fwrite(tor_perslice, file.path(out_dir, "R9_C3_4_v2_fast_toroidal_hotspot_count_perslice.csv"))
fwrite(tor_summary, file.path(out_dir, "R9_C3_4_v2_fast_toroidal_hotspot_count_summary.csv"))

cat("\n== CLR observed summary, k=6:\n")
print(clr_summary[order(response, predictor)])
cat("\n== Toroidal hotspot count summary, k=6 MES-lineage:\n")
print(tor_summary[order(predictor)])
cat("\n[STOP C3-4 v2 fast gate] Review whether CLR neuron_control is near zero and whether toroidal count gate passes.\n")
