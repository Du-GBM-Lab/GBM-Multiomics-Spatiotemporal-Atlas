#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-4 v2 audit: method-1 replacements + hotspot toroidal null audit
# Purpose:
#   v2 supersedes the old subset-based method 1 and the cropped hotspot null.
#   This script only checks the two acceptance gates:
#     1) raw no-subset and full-composition CLR associations must bring
#        neuron_control back near zero.
#     2) hotspot null must preserve hotspot counts:
#        shifted/original median >= 0.95 and no slice < 0.90.
#
# STOP:
#   No Jaccard, no adjacency statistics, no figures, no biological claim.
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
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

responses <- c("MES-lineage", "MES-V", "MES-I")
predictors <- c("vascular_niche", "myeloid_niche", "neuron_control_niche")
k_values <- c(6L, 12L)
hotspot_cut <- 0.10
max_delta <- 18L
n_shift_target <- 999L
min_mapped_spots <- 50L
eps <- 1e-6

scores <- fread(score_path)
old_shifts <- fread(shift_path)
weights <- as.data.table(qs2::qs_read(weights_path))

stopifnot(all(c("spot_id", "slice", "x", "y", "k", responses, predictors) %in% names(scores)))
stopifnot(all(c("spot_id", "slice", "x", "y") %in% names(weights)))

meta_cols <- c("spot_id", "slice", "image", "x", "y")
celltype_cols <- setdiff(names(weights), meta_cols)
stopifnot(all(c("Subtype3", "Subtype4", "Endothelial", "Mural cells",
                "Macrophages", "Microglial", "Monocytes", "Neurons") %in% celltype_cols))

safe_spearman <- function(x, y, min_pairs = 50L) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < min_pairs) return(NA_real_)
  x <- x[ok]; y <- y[ok]
  if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

emp_p <- function(obs, nul, side = c("two.sided", "greater", "less")) {
  side <- match.arg(side)
  nul <- nul[is.finite(nul)]
  if (!is.finite(obs) || !length(nul)) return(NA_real_)
  if (side == "two.sided") return((1 + sum(abs(nul) >= abs(obs))) / (length(nul) + 1))
  if (side == "greater") return((1 + sum(nul >= obs)) / (length(nul) + 1))
  (1 + sum(nul <= obs)) / (length(nul) + 1)
}

wilcox_p <- function(x, alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, alternative = alternative, exact = FALSE)$p.value)
}

sign_p <- function(x, alternative = c("two.sided", "greater", "less")) {
  alternative <- match.arg(alternative)
  x <- x[is.finite(x) & x != 0]
  if (!length(x)) return(NA_real_)
  binom.test(sum(x > 0), length(x), p = 0.5, alternative = alternative)$p.value
}

build_expanded_shifts <- function(dt, old_dt) {
  xy <- as.matrix(dt[, .(x, y)])
  fitx <- lm(shift_x ~ 0 + delta_a + delta_b, data = old_dt)
  fity <- lm(shift_y ~ 0 + delta_a + delta_b, data = old_dt)
  b1 <- c(coef(fitx)[["delta_a"]], coef(fity)[["delta_a"]])
  b2 <- c(coef(fitx)[["delta_b"]], coef(fity)[["delta_b"]])
  grid <- CJ(delta_a = -max_delta:max_delta, delta_b = -max_delta:max_delta)
  grid <- grid[!(delta_a == 0 & delta_b == 0)]
  grid[, `:=`(
    shift_x = delta_a * b1[1] + delta_b * b2[1],
    shift_y = delta_a * b1[2] + delta_b * b2[2]
  )]
  grid[, shift_dist := sqrt(shift_x^2 + shift_y^2)]
  max_match_dist <- median(old_dt$max_match_dist, na.rm = TRUE)
  coverage <- numeric(nrow(grid))
  mapped <- integer(nrow(grid))
  for (ii in seq_len(nrow(grid))) {
    nn <- FNN::get.knnx(xy, cbind(xy[, 1] + grid$shift_x[ii],
                                  xy[, 2] + grid$shift_y[ii]), k = 1)
    valid <- nn$nn.dist[, 1] <= max_match_dist
    coverage[ii] <- mean(valid)
    mapped[ii] <- sum(valid)
  }
  grid[, `:=`(coverage = coverage,
              mapped_spots = mapped,
              max_match_dist = max_match_dist,
              b1x = b1[1], b1y = b1[2], b2x = b2[1], b2y = b2[2])]
  grid <- grid[mapped_spots >= min_mapped_spots]
  setorder(grid, -coverage, shift_dist)
  grid[seq_len(min(.N, n_shift_target))]
}

build_shift_maps <- function(dt, shifts) {
  xy <- as.matrix(dt[, .(x, y)])
  lapply(seq_len(nrow(shifts)), function(ii) {
    nn <- FNN::get.knnx(xy, cbind(xy[, 1] + shifts$shift_x[ii],
                                  xy[, 2] + shifts$shift_y[ii]), k = 1)
    valid <- nn$nn.dist[, 1] <= shifts$max_match_dist[ii]
    list(index = nn$nn.index[, 1], valid = valid)
  })
}

make_neighbor_mean <- function(dt, feature, k_value) {
  xy <- as.matrix(dt[, .(x, y)])
  threshold <- median(dt$distance_threshold, na.rm = TRUE)
  nn <- FNN::get.knnx(xy, xy, k = min(nrow(dt), 32L))
  out <- rep(NA_real_, nrow(dt))
  for (ii in seq_len(nrow(dt))) {
    idx <- nn$nn.index[ii, ]
    dst <- nn$nn.dist[ii, ]
    keep <- idx != ii & dst <= threshold
    if (sum(keep) > 0) out[ii] <- mean(feature[idx[keep]], na.rm = TRUE)
  }
  out
}

make_clr_features <- function(wdt) {
  mat <- as.matrix(wdt[, ..celltype_cols])
  logm <- log(mat + eps)
  gm <- rowMeans(logm)
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

run_association <- function(dt, maps, response, predictor) {
  y <- dt[[response]]
  x <- dt[[predictor]]
  obs <- safe_spearman(y, x)
  nul <- vapply(maps, function(mp) {
    xs <- rep(NA_real_, length(x))
    xs[mp$valid] <- x[mp$index[mp$valid]]
    safe_spearman(y[mp$valid], xs[mp$valid])
  }, numeric(1))
  data.table(
    response = response,
    predictor = predictor,
    observed_rho = obs,
    n_null_shifts = sum(is.finite(nul)),
    null_rho_median = median(nul, na.rm = TRUE),
    p_emp_two_sided = emp_p(obs, nul, "two.sided"),
    p_emp_greater = emp_p(obs, nul, "greater"),
    p_emp_less = emp_p(obs, nul, "less")
  )
}

summarize_assoc <- function(dt) {
  dt[, .(
    n_slices = sum(is.finite(observed_rho)),
    median_rho = median(observed_rho, na.rm = TRUE),
    q25_rho = as.numeric(quantile(observed_rho, 0.25, na.rm = TRUE)),
    q75_rho = as.numeric(quantile(observed_rho, 0.75, na.rm = TRUE)),
    n_positive = sum(observed_rho > 0, na.rm = TRUE),
    n_negative = sum(observed_rho < 0, na.rm = TRUE),
    sign_p_greater = sign_p(observed_rho, "greater"),
    wilcox_p_greater = wilcox_p(observed_rho, "greater"),
    median_emp_p_two_sided = median(p_emp_two_sided, na.rm = TRUE),
    n_slice_emp_p_lt_0p05 = sum(p_emp_two_sided < 0.05, na.rm = TRUE),
    median_null_shifts = median(n_null_shifts, na.rm = TRUE),
    min_null_shifts = min(n_null_shifts, na.rm = TRUE)
  ), by = .(method, k, response, predictor)]
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
  data.table(
    spot_id = dt$spot_id,
    i = ij[, 1],
    j = ij[, 2],
    recon_error = sqrt(rowSums((xy - recon)^2))
  )
}

toroidal_count_audit <- function(dt, shifts, response, predictor) {
  lat <- infer_lattice(dt, shifts)
  ## Keep one spot per rounded lattice coordinate. Duplicates indicate imperfect
  ## lattice recovery and are excluded from toroidal count audit for honesty.
  lat[, row_id := seq_len(.N)]
  dup <- duplicated(lat[, .(i, j)]) | duplicated(lat[, .(i, j)], fromLast = TRUE)
  lat_u <- lat[!dup]
  key <- paste(lat_u$i, lat_u$j)
  map_row <- setNames(lat_u$row_id, key)
  irange <- range(lat_u$i)
  jrange <- range(lat_u$j)
  ni <- diff(irange) + 1L
  nj <- diff(jrange) + 1L

  x <- dt[[predictor]]
  thr <- as.numeric(quantile(x, 1 - hotspot_cut, na.rm = TRUE))
  hot <- is.finite(x) & x >= thr
  orig <- sum(hot)

  ratios <- numeric(nrow(shifts))
  valid_targets <- integer(nrow(shifts))
  for (ss in seq_len(nrow(shifts))) {
    ti <- ((lat_u$i + shifts$delta_a[ss] - irange[1]) %% ni) + irange[1]
    tj <- ((lat_u$j + shifts$delta_b[ss] - jrange[1]) %% nj) + jrange[1]
    target <- map_row[paste(ti, tj)]
    ok <- !is.na(target)
    shifted <- rep(FALSE, nrow(dt))
    shifted[lat_u$row_id[ok]] <- hot[as.integer(target[ok])]
    ratios[ss] <- sum(shifted) / orig
    valid_targets[ss] <- sum(ok)
  }
  data.table(
    response = response,
    predictor = predictor,
    hotspot_cut = hotspot_cut,
    n_spots = nrow(dt),
    n_lattice_unique = nrow(lat_u),
    n_lattice_duplicate_spots = sum(dup),
    median_recon_error = median(lat$recon_error),
    q95_recon_error = as.numeric(quantile(lat$recon_error, 0.95)),
    observed_hotspots = orig,
    n_null_shifts = length(ratios),
    count_ratio_median = median(ratios, na.rm = TRUE),
    count_ratio_min = min(ratios, na.rm = TRUE),
    count_ratio_q25 = as.numeric(quantile(ratios, 0.25, na.rm = TRUE)),
    count_ratio_q75 = as.numeric(quantile(ratios, 0.75, na.rm = TRUE)),
    valid_target_median = median(valid_targets),
    valid_target_min = min(valid_targets)
  )
}

message("== C3-4 v2 audit input: scores rows=", nrow(scores),
        "; weights rows=", nrow(weights), "; full composition columns=", length(celltype_cols))

clr_all <- make_clr_features(weights)
raw_assoc <- list()
clr_assoc <- list()
toroidal <- list()
idx_raw <- idx_clr <- idx_tor <- 0L

for (sl in sort(unique(scores$slice))) {
  for (kk in k_values) {
    dt_raw <- scores[slice == sl & k == kk]
    old_dt <- old_shifts[slice == sl]
    if (!nrow(dt_raw) || !nrow(old_dt)) next
    shifts <- build_expanded_shifts(dt_raw, old_dt)
    maps <- build_shift_maps(dt_raw, shifts)

    ## raw no-subset: response and predictor are the existing C3-1 fields.
    for (resp in responses) {
      for (pred in predictors) {
        idx_raw <- idx_raw + 1L
        rr <- run_association(dt_raw, maps, resp, pred)
        rr[, `:=`(method = "raw_no_subset", slice = sl, k = kk)]
        raw_assoc[[idx_raw]] <- rr
      }
    }

    ## CLR: focal response from CLR composition; predictor is neighbor mean of CLR field.
    dt_clr <- merge(
      dt_raw[, .(spot_id, slice, x, y, k, distance_threshold)],
      clr_all[slice == sl],
      by = c("spot_id", "slice", "x", "y"),
      all.x = TRUE,
      sort = FALSE
    )
    for (pred in predictors) {
      dt_clr[[pred]] <- make_neighbor_mean(dt_clr, dt_clr[[pred]], kk)
    }
    for (resp in responses) {
      for (pred in predictors) {
        idx_clr <- idx_clr + 1L
        cr <- run_association(dt_clr, maps, resp, pred)
        cr[, `:=`(method = "clr_full_composition", slice = sl, k = kk)]
        clr_assoc[[idx_clr]] <- cr
      }
    }

    ## Toroidal hotspot-count audit. This is only a count/null audit, not Jaccard.
    for (resp in responses) {
      for (pred in predictors) {
        idx_tor <- idx_tor + 1L
        ta <- toroidal_count_audit(dt_raw, shifts, resp, pred)
        ta[, `:=`(slice = sl, k = kk)]
        toroidal[[idx_tor]] <- ta
      }
    }
    message("  done slice=", sl, " k=", kk, " shifts=", nrow(shifts))
  }
}

assoc <- rbindlist(c(raw_assoc, clr_assoc), use.names = TRUE, fill = TRUE)
assoc[, sensitivity := "all_slices"]
assoc_no265 <- copy(assoc[slice != "#UKF265_T_ST"])
assoc_no265[, sensitivity := "exclude_UKF265"]
assoc2 <- rbindlist(list(assoc, assoc_no265), use.names = TRUE)
summary_assoc <- summarize_assoc(assoc2)

toroidal <- rbindlist(toroidal, use.names = TRUE, fill = TRUE)
toroidal[, pass_count_gate := count_ratio_median >= 0.95 & count_ratio_min >= 0.90]
toroidal_summary <- toroidal[, .(
  n_slices = .N,
  global_median_count_ratio = median(count_ratio_median, na.rm = TRUE),
  min_slice_count_ratio_median = min(count_ratio_median, na.rm = TRUE),
  min_slice_count_ratio_min = min(count_ratio_min, na.rm = TRUE),
  n_slice_pass = sum(pass_count_gate, na.rm = TRUE),
  max_duplicate_spots = max(n_lattice_duplicate_spots, na.rm = TRUE),
  median_q95_recon_error = median(q95_recon_error, na.rm = TRUE)
), by = .(k, response, predictor)]

fwrite(assoc, file.path(out_dir, "R9_C3_4_v2_method1_replacements_perslice.csv"))
fwrite(summary_assoc, file.path(out_dir, "R9_C3_4_v2_method1_replacements_summary.csv"))
fwrite(toroidal, file.path(out_dir, "R9_C3_4_v2_toroidal_hotspot_null_audit.csv"))
fwrite(toroidal_summary, file.path(out_dir, "R9_C3_4_v2_toroidal_hotspot_null_summary.csv"))

cat("\n== Method-1 replacement summary, k=6 focus pairs:\n")
print(summary_assoc[k == 6 & sensitivity == "all_slices" &
                      response %in% responses &
                      predictor %in% c("vascular_niche", "neuron_control_niche")]
      [order(method, response, predictor)])

cat("\n== Method-1 replacement summary, k=6 focus pairs, excluding UKF265:\n")
print(summary_assoc[k == 6 & sensitivity == "exclude_UKF265" &
                      response %in% responses &
                      predictor %in% c("vascular_niche", "neuron_control_niche")]
      [order(method, response, predictor)])

cat("\n== Toroidal hotspot null count audit, k=6 MES-lineage focus:\n")
print(toroidal_summary[k == 6 & response == "MES-lineage" &
                         predictor %in% c("vascular_niche", "myeloid_niche", "neuron_control_niche")]
      [order(predictor)])

cat("\n[STOP C3-4 v2 audit] Acceptance gates:\n",
    " 1) neuron_control should be near zero for method-1 replacements;\n",
    " 2) toroidal hotspot count ratio median >=0.95 and no slice <0.90.\n",
    "Do not run Jaccard/adjacency until both gates are reviewed.\n")
