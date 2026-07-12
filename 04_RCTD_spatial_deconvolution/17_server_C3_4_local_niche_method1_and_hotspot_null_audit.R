#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-4 step 1: local conditional association + hotspot-mask null audit
# Purpose:
#   1) Build an expanded vector-shift null from the audited C3-2a2 hex basis.
#   2) Run method 1: conditional/local Spearman association after removing
#      double-low background spots.
#   3) Audit method 2 null construction before running hotspot Jaccard:
#      shift whole hotspot masks and report null mask count/coverage behavior.
#
# Boundaries:
#   - Slice is the statistical unit. No pooled-spot inference.
#   - Vector-shift null is reported explicitly; bootstrap is not a null.
#   - Response is focal MES weight; predictor is neighborhood niche score.
#   - Claims are not made here. This is a STOP table-generation step.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(FNN)
})

set.seed(1)

## Local execution uses ASCII relative paths from the R9 working directory.
base_dir <- getwd()
score_path <- file.path(base_dir, "tables/C3_niche_preflight/R9_C3_1_neighborhood_niche_scores.csv")
shift_path <- file.path(base_dir, "tables/C3_niche_preflight/R9_C3_2a2_vector_shift_coverage.csv")
out_dir <- file.path(base_dir, "tables/C3_4_local_niche")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

responses <- c("MES-lineage", "MES-V", "MES-I")
predictors <- c("vascular_niche", "myeloid_niche", "neuron_control_niche")
subset_cuts <- c(tertile = 1/3, quartile = 1/4, decile = 1/10)
hotspot_cut <- 0.10
k_values <- c(6L, 12L)

## Expanded shift grid: enough candidates for 999-shift null where geometry allows.
max_delta <- 18L
n_shift_target <- 999L
min_mapped_spots <- 50L

scores <- fread(score_path)
old_shifts <- fread(shift_path)
need <- c("spot_id", "slice", "x", "y", "k", responses, predictors)
stopifnot(all(need %in% names(scores)))
stopifnot(all(c("slice", "delta_a", "delta_b", "shift_x", "shift_y",
                "coverage", "max_match_dist") %in% names(old_shifts)))

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

make_expanded_shifts <- function(dt, old_dt) {
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
              max_match_dist = max_match_dist)]
  grid <- grid[mapped_spots >= min_mapped_spots]
  setorder(grid, -coverage, shift_dist)
  grid[, shift_id := seq_len(.N)]
  grid[seq_len(min(.N, n_shift_target))]
}

make_shift_maps <- function(dt, shifts) {
  xy <- as.matrix(dt[, .(x, y)])
  lapply(seq_len(nrow(shifts)), function(ii) {
    nn <- FNN::get.knnx(xy, cbind(xy[, 1] + shifts$shift_x[ii],
                                  xy[, 2] + shifts$shift_y[ii]), k = 1)
    valid <- nn$nn.dist[, 1] <= shifts$max_match_dist[ii]
    list(index = nn$nn.index[, 1], valid = valid)
  })
}

observed_conditional <- function(y, x, cut) {
  y_thr <- as.numeric(quantile(y, 1 - cut, na.rm = TRUE))
  x_thr <- as.numeric(quantile(x, 1 - cut, na.rm = TRUE))
  use <- (y >= y_thr) | (x >= x_thr)
  list(rho = safe_spearman(y[use], x[use]),
       n_subset = sum(use, na.rm = TRUE),
       y_thr = y_thr, x_thr = x_thr)
}

null_conditional <- function(y, x, maps, cut) {
  vapply(maps, function(mp) {
    xs <- rep(NA_real_, length(x))
    xs[mp$valid] <- x[mp$index[mp$valid]]
    ok <- is.finite(xs) & is.finite(y)
    if (sum(ok) < min_mapped_spots) return(NA_real_)
    y_thr <- as.numeric(quantile(y[ok], 1 - cut, na.rm = TRUE))
    x_thr <- as.numeric(quantile(xs[ok], 1 - cut, na.rm = TRUE))
    use <- ok & ((y >= y_thr) | (xs >= x_thr))
    safe_spearman(y[use], xs[use])
  }, numeric(1))
}

hotspot_mask_audit <- function(dt, maps, response, predictor, cut = hotspot_cut) {
  y <- dt[[response]]
  x <- dt[[predictor]]
  x_thr <- as.numeric(quantile(x, 1 - cut, na.rm = TRUE))
  x_hot <- is.finite(x) & x >= x_thr
  null_counts <- vapply(maps, function(mp) {
    xs <- rep(FALSE, length(x_hot))
    xs[mp$valid] <- x_hot[mp$index[mp$valid]]
    sum(xs)
  }, integer(1))
  data.table(
    response = response,
    predictor = predictor,
    hotspot_cut = cut,
    observed_predictor_hotspots = sum(x_hot),
    null_hotspot_count_median = median(null_counts, na.rm = TRUE),
    null_hotspot_count_q25 = as.numeric(quantile(null_counts, 0.25, na.rm = TRUE)),
    null_hotspot_count_q75 = as.numeric(quantile(null_counts, 0.75, na.rm = TRUE)),
    null_hotspot_count_min = min(null_counts, na.rm = TRUE),
    null_hotspot_count_max = max(null_counts, na.rm = TRUE),
    null_count_ratio_median = median(null_counts) / sum(x_hot),
    null_count_ratio_min = min(null_counts) / sum(x_hot)
  )
}

message("== C3-4 input rows: ", nrow(scores),
        "; slices=", length(unique(scores$slice)),
        "; k=", paste(sort(unique(scores$k)), collapse = ","))

shift_audit <- list()
method1 <- list()
mask_audit <- list()
ii_shift <- ii_m1 <- ii_mask <- 0L

for (sl in sort(unique(scores$slice))) {
  for (kk in k_values) {
    dt <- scores[slice == sl & k == kk]
    old_dt <- old_shifts[slice == sl]
    if (!nrow(dt) || !nrow(old_dt)) next

    shifts <- make_expanded_shifts(dt, old_dt)
    maps <- make_shift_maps(dt, shifts)
    ii_shift <- ii_shift + 1L
    shift_audit[[ii_shift]] <- data.table(
      slice = sl,
      k = kk,
      n_spots = nrow(dt),
      n_candidate_shifts_retained = nrow(shifts),
      target_shifts = n_shift_target,
      min_mapped_spots_rule = min_mapped_spots,
      min_coverage_used = min(shifts$coverage),
      median_coverage_used = median(shifts$coverage),
      min_mapped_spots_used = min(shifts$mapped_spots),
      median_mapped_spots_used = median(shifts$mapped_spots)
    )

    for (resp in responses) {
      y <- dt[[resp]]
      for (pred in predictors) {
        x <- dt[[pred]]

        ## Method 2 pre-audit for the main top-decile hotspot mask null.
        ii_mask <- ii_mask + 1L
        ma <- hotspot_mask_audit(dt, maps, resp, pred, hotspot_cut)
        ma[, `:=`(slice = sl, k = kk, n_spots = nrow(dt),
                  n_null_shifts = length(maps))]
        mask_audit[[ii_mask]] <- ma

        ## Method 1 conditional association, main and sensitivity cuts.
        for (cut_name in names(subset_cuts)) {
          cut <- subset_cuts[[cut_name]]
          obs <- observed_conditional(y, x, cut)
          nul <- null_conditional(y, x, maps, cut)
          ii_m1 <- ii_m1 + 1L
          method1[[ii_m1]] <- data.table(
            slice = sl,
            k = kk,
            response = resp,
            predictor = pred,
            subset_rule = cut_name,
            top_fraction = cut,
            observed_rho = obs$rho,
            n_subset = obs$n_subset,
            response_threshold = obs$y_thr,
            predictor_threshold = obs$x_thr,
            n_null_shifts = sum(is.finite(nul)),
            null_rho_median = median(nul, na.rm = TRUE),
            null_rho_q25 = as.numeric(quantile(nul, 0.25, na.rm = TRUE)),
            null_rho_q75 = as.numeric(quantile(nul, 0.75, na.rm = TRUE)),
            p_emp_two_sided = emp_p(obs$rho, nul, "two.sided"),
            p_emp_greater = emp_p(obs$rho, nul, "greater"),
            p_emp_less = emp_p(obs$rho, nul, "less")
          )
        }
      }
    }
    message("  done slice=", sl, " k=", kk, " shifts=", nrow(shifts))
  }
}

shift_audit <- rbindlist(shift_audit, use.names = TRUE, fill = TRUE)
method1 <- rbindlist(method1, use.names = TRUE, fill = TRUE)
mask_audit <- rbindlist(mask_audit, use.names = TRUE, fill = TRUE)

summarize_method1 <- function(dt) {
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
    min_null_shifts = min(n_null_shifts, na.rm = TRUE),
    median_subset_n = median(n_subset, na.rm = TRUE)
  ), by = .(k, response, predictor, subset_rule, top_fraction)]
}

summary_method1 <- summarize_method1(method1)
summary_method1_no_ukf265 <- summarize_method1(method1[slice != "#UKF265_T_ST"])
summary_method1_no_ukf265[, sensitivity := "exclude_UKF265"]
summary_method1[, sensitivity := "all_slices"]
summary_method1_all <- rbindlist(list(summary_method1, summary_method1_no_ukf265), use.names = TRUE)

fwrite(shift_audit, file.path(out_dir, "R9_C3_4_step1_expanded_shift_audit.csv"))
fwrite(mask_audit, file.path(out_dir, "R9_C3_4_step1_hotspot_mask_null_audit.csv"))
fwrite(method1, file.path(out_dir, "R9_C3_4_method1_conditional_association_perslice.csv"))
fwrite(summary_method1_all, file.path(out_dir, "R9_C3_4_method1_conditional_association_summary.csv"))

cat("\n== Expanded vector-shift audit:\n")
print(shift_audit[, .(
  n_slice_k = .N,
  min_retained_shifts = min(n_candidate_shifts_retained),
  median_retained_shifts = median(n_candidate_shifts_retained),
  min_coverage = min(min_coverage_used),
  median_coverage = median(median_coverage_used),
  min_mapped = min(min_mapped_spots_used)
)])

cat("\n== Slices/k with fewer than target shifts or low minimum coverage:\n")
print(shift_audit[n_candidate_shifts_retained < n_shift_target | min_coverage_used < 0.2]
      [order(slice, k)])

cat("\n== Hotspot mask null audit, main pairs k=6:\n")
print(mask_audit[k == 6 & response == "MES-lineage" &
                   predictor %in% c("vascular_niche", "myeloid_niche", "neuron_control_niche")]
      [order(predictor, slice),
       .(slice, predictor, observed_predictor_hotspots, null_hotspot_count_median,
         null_count_ratio_median, null_count_ratio_min, n_null_shifts)])

cat("\n== Method 1 conditional association summary, main tertile k=6:\n")
print(summary_method1_all[k == 6 & subset_rule == "tertile" &
                            response %in% responses &
                            predictor %in% c("vascular_niche", "myeloid_niche", "neuron_control_niche")]
      [order(sensitivity, response, predictor)])

cat("\n[STOP C3-4 step 1] Report expanded-shift audit, hotspot-mask null audit,",
    "and method-1 conditional association summary. Do not proceed to Jaccard/",
    "adjacency figures until this STOP is reviewed.\n")
