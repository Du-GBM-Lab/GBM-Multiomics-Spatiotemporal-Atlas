#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-4 v3-B: Ripley-style distance statistics with random-labeling null
# Purpose:
#   Test whether high-MES spots are closer to high-vascular spots than expected
#   under a point-count preserving random-labeling null on the observed tissue
#   spot positions. CSR is explicitly not used.
#
# Statistics:
#   1) crossK_count: mean number of predictor-hot spots within the k=6/k=12
#      physical radius of each MES-hot spot. Higher = closer/enriched.
#   2) mean_nearest_distance: mean distance from each MES-hot spot to nearest
#      predictor-hot spot. Lower = closer.
#
# STOP:
#   Generate source tables only. No figures and no causal/spatial-contact claim.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
})

set.seed(1)

base_dir <- getwd()
weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
score_path <- file.path(base_dir, "tables/C3_niche_preflight/R9_C3_1_neighborhood_niche_scores.csv")
out_dir <- file.path(base_dir, "tables/C3_4_local_niche")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_perm <- 999L
k_values <- c(6L, 12L)
top_fracs <- c(top05 = 0.05, top10 = 0.10, top20 = 0.20)
responses <- c("MES-lineage", "MES-V", "MES-I")
predictors <- c("vascular", "myeloid", "neuron_control")

weights <- as.data.table(qs2::qs_read(weights_path))
scores <- fread(score_path)
stopifnot(all(c("spot_id", "slice", "x", "y", "Subtype3", "Subtype4",
                "Endothelial", "Mural cells", "Macrophages", "Microglial",
                "Monocytes", "Neurons") %in% names(weights)))

feat <- data.table(
  spot_id = weights$spot_id,
  slice = weights$slice,
  x = weights$x,
  y = weights$y,
  `MES-lineage` = weights$Subtype3 + weights$Subtype4,
  `MES-V` = weights$Subtype3,
  `MES-I` = weights$Subtype4,
  vascular = weights$Endothelial + weights[["Mural cells"]],
  myeloid = weights$Macrophages + weights$Microglial + weights$Monocytes,
  neuron_control = weights$Neurons
)

thresholds <- unique(scores[, .(slice, k, distance_threshold)])
thresholds <- thresholds[, .(distance_threshold = median(distance_threshold, na.rm = TRUE)),
                         by = .(slice, k)]

make_hot <- function(x, frac) {
  thr <- as.numeric(quantile(x, 1 - frac, na.rm = TRUE))
  is.finite(x) & x >= thr
}

calc_stats <- function(D, A, B, radius) {
  nA <- sum(A)
  nB <- sum(B)
  if (nA == 0 || nB == 0) {
    return(c(crossK_count = NA_real_, mean_nearest_distance = NA_real_))
  }
  DB <- D[A, B, drop = FALSE]
  c(crossK_count = mean(rowSums(DB <= radius)),
    mean_nearest_distance = mean(apply(DB, 1, min)))
}

emp_p <- function(obs, nul, side = c("greater", "less", "two.sided")) {
  side <- match.arg(side)
  nul <- nul[is.finite(nul)]
  if (!is.finite(obs) || !length(nul)) return(NA_real_)
  if (side == "greater") return((1 + sum(nul >= obs)) / (length(nul) + 1))
  if (side == "less") return((1 + sum(nul <= obs)) / (length(nul) + 1))
  (1 + sum(abs(nul) >= abs(obs))) / (length(nul) + 1)
}

run_one <- function(dt, D, radius, response, predictor, frac, frac_name) {
  A <- make_hot(dt[[response]], frac)
  B <- make_hot(dt[[predictor]], frac)
  n <- nrow(dt)
  nA <- sum(A)
  nB <- sum(B)
  obs <- calc_stats(D, A, B, radius)

  null_K <- numeric(n_perm)
  null_NN <- numeric(n_perm)
  count_ok <- logical(n_perm)
  for (pp in seq_len(n_perm)) {
    Ai <- rep(FALSE, n)
    Bi <- rep(FALSE, n)
    Ai[sample.int(n, nA)] <- TRUE
    Bi[sample.int(n, nB)] <- TRUE
    count_ok[pp] <- sum(Ai) == nA && sum(Bi) == nB
    st <- calc_stats(D, Ai, Bi, radius)
    null_K[pp] <- st[["crossK_count"]]
    null_NN[pp] <- st[["mean_nearest_distance"]]
  }

  data.table(
    response = response,
    predictor = predictor,
    top_rule = frac_name,
    top_fraction = frac,
    n_spots = n,
    radius = radius,
    n_response_hot = nA,
    n_predictor_hot = nB,
    count_preserved_all_perm = all(count_ok),
    observed_crossK_count = obs[["crossK_count"]],
    null_crossK_median = median(null_K, na.rm = TRUE),
    null_crossK_q25 = as.numeric(quantile(null_K, 0.25, na.rm = TRUE)),
    null_crossK_q75 = as.numeric(quantile(null_K, 0.75, na.rm = TRUE)),
    delta_crossK = obs[["crossK_count"]] - median(null_K, na.rm = TRUE),
    p_emp_crossK_greater = emp_p(obs[["crossK_count"]], null_K, "greater"),
    observed_mean_nn_dist = obs[["mean_nearest_distance"]],
    null_mean_nn_dist_median = median(null_NN, na.rm = TRUE),
    null_mean_nn_dist_q25 = as.numeric(quantile(null_NN, 0.25, na.rm = TRUE)),
    null_mean_nn_dist_q75 = as.numeric(quantile(null_NN, 0.75, na.rm = TRUE)),
    delta_nn_closer = median(null_NN, na.rm = TRUE) - obs[["mean_nearest_distance"]],
    p_emp_nn_less = emp_p(obs[["mean_nearest_distance"]], null_NN, "less"),
    n_perm = length(null_K)
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

summarize_one <- function(dt) {
  dt[, .(
    n_slices = .N,
    median_delta_crossK = median(delta_crossK, na.rm = TRUE),
    q25_delta_crossK = as.numeric(quantile(delta_crossK, 0.25, na.rm = TRUE)),
    q75_delta_crossK = as.numeric(quantile(delta_crossK, 0.75, na.rm = TRUE)),
    n_crossK_positive = sum(delta_crossK > 0, na.rm = TRUE),
    sign_p_crossK_greater = sign_p(delta_crossK, "greater"),
    wilcox_p_crossK_greater = wilcox_p(delta_crossK, "greater"),
    median_p_crossK = median(p_emp_crossK_greater, na.rm = TRUE),
    n_slice_crossK_p_lt_0p05 = sum(p_emp_crossK_greater < 0.05, na.rm = TRUE),
    median_delta_nn_closer = median(delta_nn_closer, na.rm = TRUE),
    q25_delta_nn_closer = as.numeric(quantile(delta_nn_closer, 0.25, na.rm = TRUE)),
    q75_delta_nn_closer = as.numeric(quantile(delta_nn_closer, 0.75, na.rm = TRUE)),
    n_nn_closer_positive = sum(delta_nn_closer > 0, na.rm = TRUE),
    sign_p_nn_greater = sign_p(delta_nn_closer, "greater"),
    wilcox_p_nn_greater = wilcox_p(delta_nn_closer, "greater"),
    median_p_nn = median(p_emp_nn_less, na.rm = TRUE),
    n_slice_nn_p_lt_0p05 = sum(p_emp_nn_less < 0.05, na.rm = TRUE),
    count_preserved_all = all(count_preserved_all_perm),
    median_response_hot = median(n_response_hot),
    median_predictor_hot = median(n_predictor_hot),
    median_perm = median(n_perm)
  ), by = .(k, response, predictor, top_rule, top_fraction)]
}

message("== Running v3-B random-labeling nearest/crossK audit; CSR is not used.")
res <- list()
idx <- 0L

for (sl in sort(unique(feat$slice))) {
  dt <- feat[slice == sl]
  xy <- as.matrix(dt[, .(x, y)])
  D <- as.matrix(dist(xy))
  diag(D) <- 0

  for (kk in k_values) {
    radius <- thresholds[slice == sl & k == kk, distance_threshold][1]
    for (frac_name in names(top_fracs)) {
      frac <- top_fracs[[frac_name]]
      for (resp in responses) {
        for (pred in predictors) {
          idx <- idx + 1L
          ans <- run_one(dt, D, radius, resp, pred, frac, frac_name)
          ans[, `:=`(slice = sl, k = kk)]
          res[[idx]] <- ans
        }
      }
    }
  }
  message("  done slice=", sl, " n=", nrow(dt))
}

perslice <- rbindlist(res, use.names = TRUE, fill = TRUE)
setcolorder(perslice, c("slice", "k", "response", "predictor", "top_rule",
                        setdiff(names(perslice), c("slice", "k", "response", "predictor", "top_rule"))))

summary <- summarize_one(perslice)
summary[, sensitivity := "all_slices"]
summary_no265 <- summarize_one(perslice[slice != "#UKF265_T_ST"])
summary_no265[, sensitivity := "exclude_UKF265"]
summary_all <- rbindlist(list(summary, summary_no265), use.names = TRUE)

fwrite(perslice, file.path(out_dir, "R9_C3_4_v3B_Ripley_random_labeling_perslice.csv"))
fwrite(summary_all, file.path(out_dir, "R9_C3_4_v3B_Ripley_random_labeling_summary.csv"))

cat("\n== v3-B focus summary, k=6 top10 all slices:\n")
print(summary_all[k == 6 & top_rule == "top10" & sensitivity == "all_slices" &
                    response %in% responses &
                    predictor %in% c("vascular", "myeloid", "neuron_control")]
      [order(response, predictor)])

cat("\n== v3-B focus summary, k=12 top10 all slices:\n")
print(summary_all[k == 12 & top_rule == "top10" & sensitivity == "all_slices" &
                    response %in% responses &
                    predictor %in% c("vascular", "myeloid", "neuron_control")]
      [order(response, predictor)])

cat("\n[STOP C3-4 v3-B Ripley/random-labeling] Report observed effect deltas,",
    " null p-values, count-preservation audit, and neuron_control gate.",
    " Do not make causal/contact claims.\n")
