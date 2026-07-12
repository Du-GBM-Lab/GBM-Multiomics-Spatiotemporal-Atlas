#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-4 v3-B fast main gate: MES-lineage vs vascular/neuron random-labeling
# Purpose:
#   Minimal Ripley-style STOP gate after the full sensitivity grid was too slow.
#   Runs only the primary response and gate predictors:
#     response = MES-lineage
#     predictors = vascular, neuron_control
#     top fraction = top10
#     k scale = 6 and 12
#   Null is random labeling on observed spot positions; CSR is not used.
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

n_perm <- 999L
k_values <- c(6L, 12L)
top_fraction <- 0.10
response <- "MES-lineage"
predictors <- c("vascular", "neuron_control")

weights <- as.data.table(qs2::qs_read(weights_path))
scores <- fread(score_path)
feat <- data.table(
  spot_id = weights$spot_id,
  slice = weights$slice,
  x = weights$x,
  y = weights$y,
  `MES-lineage` = weights$Subtype3 + weights$Subtype4,
  vascular = weights$Endothelial + weights[["Mural cells"]],
  neuron_control = weights$Neurons
)
thresholds <- unique(scores[, .(slice, k, distance_threshold)])
thresholds <- thresholds[, .(distance_threshold = median(distance_threshold, na.rm = TRUE)),
                         by = .(slice, k)]

make_hot <- function(x, frac) {
  x >= as.numeric(quantile(x, 1 - frac, na.rm = TRUE))
}

calc_stats <- function(D, A_idx, B_idx, radius) {
  DB <- D[A_idx, B_idx, drop = FALSE]
  c(
    crossK_count = mean(rowSums(DB <= radius)),
    mean_nearest_distance = mean(apply(DB, 1, min))
  )
}

emp_p <- function(obs, nul, side = c("greater", "less")) {
  side <- match.arg(side)
  if (side == "greater") return((1 + sum(nul >= obs, na.rm = TRUE)) / (sum(is.finite(nul)) + 1))
  (1 + sum(nul <= obs, na.rm = TRUE)) / (sum(is.finite(nul)) + 1)
}

run_one <- function(dt, D, radius, predictor) {
  A <- make_hot(dt[[response]], top_fraction)
  B <- make_hot(dt[[predictor]], top_fraction)
  n <- nrow(dt)
  nA <- sum(A)
  nB <- sum(B)
  obs <- calc_stats(D, which(A), which(B), radius)

  nullK <- numeric(n_perm)
  nullNN <- numeric(n_perm)
  count_ok <- logical(n_perm)
  for (p in seq_len(n_perm)) {
    Ai <- sample.int(n, nA)
    Bi <- sample.int(n, nB)
    count_ok[p] <- length(Ai) == nA && length(Bi) == nB
    st <- calc_stats(D, Ai, Bi, radius)
    nullK[p] <- st[["crossK_count"]]
    nullNN[p] <- st[["mean_nearest_distance"]]
  }

  data.table(
    response = response,
    predictor = predictor,
    top_rule = "top10",
    top_fraction = top_fraction,
    n_spots = n,
    n_response_hot = nA,
    n_predictor_hot = nB,
    radius = radius,
    count_preserved_all_perm = all(count_ok),
    observed_crossK_count = obs[["crossK_count"]],
    null_crossK_median = median(nullK),
    delta_crossK = obs[["crossK_count"]] - median(nullK),
    p_emp_crossK_greater = emp_p(obs[["crossK_count"]], nullK, "greater"),
    observed_mean_nn_dist = obs[["mean_nearest_distance"]],
    null_mean_nn_dist_median = median(nullNN),
    delta_nn_closer = median(nullNN) - obs[["mean_nearest_distance"]],
    p_emp_nn_less = emp_p(obs[["mean_nearest_distance"]], nullNN, "less"),
    n_perm = n_perm
  )
}

sign_p <- function(x, alternative = "greater") {
  x <- x[is.finite(x) & x != 0]
  if (!length(x)) return(NA_real_)
  binom.test(sum(x > 0), length(x), p = 0.5, alternative = alternative)$p.value
}

wilcox_p <- function(x, alternative = "greater") {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, alternative = alternative, exact = FALSE)$p.value)
}

res <- list()
i <- 0L
for (sl in sort(unique(feat$slice))) {
  dt <- feat[slice == sl]
  D <- as.matrix(dist(as.matrix(dt[, .(x, y)])))
  diag(D) <- 0
  for (kk in k_values) {
    radius <- thresholds[slice == sl & k == kk, distance_threshold][1]
    for (pred in predictors) {
      i <- i + 1L
      ans <- run_one(dt, D, radius, pred)
      ans[, `:=`(slice = sl, k = kk)]
      res[[i]] <- ans
    }
  }
  message("done ", sl)
}

perslice <- rbindlist(res, use.names = TRUE)
summary <- perslice[, .(
  n_slices = .N,
  median_delta_crossK = median(delta_crossK),
  q25_delta_crossK = as.numeric(quantile(delta_crossK, 0.25)),
  q75_delta_crossK = as.numeric(quantile(delta_crossK, 0.75)),
  n_crossK_positive = sum(delta_crossK > 0),
  sign_p_crossK_greater = sign_p(delta_crossK),
  wilcox_p_crossK_greater = wilcox_p(delta_crossK),
  median_p_crossK = median(p_emp_crossK_greater),
  n_slice_crossK_p_lt_0p05 = sum(p_emp_crossK_greater < 0.05),
  median_delta_nn_closer = median(delta_nn_closer),
  q25_delta_nn_closer = as.numeric(quantile(delta_nn_closer, 0.25)),
  q75_delta_nn_closer = as.numeric(quantile(delta_nn_closer, 0.75)),
  n_nn_closer_positive = sum(delta_nn_closer > 0),
  sign_p_nn_greater = sign_p(delta_nn_closer),
  wilcox_p_nn_greater = wilcox_p(delta_nn_closer),
  median_p_nn = median(p_emp_nn_less),
  n_slice_nn_p_lt_0p05 = sum(p_emp_nn_less < 0.05),
  count_preserved_all = all(count_preserved_all_perm),
  median_response_hot = median(n_response_hot),
  median_predictor_hot = median(n_predictor_hot)
), by = .(k, response, predictor, top_rule, top_fraction)]

summary[, sensitivity := "all_slices"]
summary_no265 <- perslice[slice != "#UKF265_T_ST", .(
  n_slices = .N,
  median_delta_crossK = median(delta_crossK),
  q25_delta_crossK = as.numeric(quantile(delta_crossK, 0.25)),
  q75_delta_crossK = as.numeric(quantile(delta_crossK, 0.75)),
  n_crossK_positive = sum(delta_crossK > 0),
  sign_p_crossK_greater = sign_p(delta_crossK),
  wilcox_p_crossK_greater = wilcox_p(delta_crossK),
  median_p_crossK = median(p_emp_crossK_greater),
  n_slice_crossK_p_lt_0p05 = sum(p_emp_crossK_greater < 0.05),
  median_delta_nn_closer = median(delta_nn_closer),
  q25_delta_nn_closer = as.numeric(quantile(delta_nn_closer, 0.25)),
  q75_delta_nn_closer = as.numeric(quantile(delta_nn_closer, 0.75)),
  n_nn_closer_positive = sum(delta_nn_closer > 0),
  sign_p_nn_greater = sign_p(delta_nn_closer),
  wilcox_p_nn_greater = wilcox_p(delta_nn_closer),
  median_p_nn = median(p_emp_nn_less),
  n_slice_nn_p_lt_0p05 = sum(p_emp_nn_less < 0.05),
  count_preserved_all = all(count_preserved_all_perm),
  median_response_hot = median(n_response_hot),
  median_predictor_hot = median(n_predictor_hot)
), by = .(k, response, predictor, top_rule, top_fraction)]
summary_no265[, sensitivity := "exclude_UKF265"]
summary_all <- rbindlist(list(summary, summary_no265), use.names = TRUE)

fwrite(perslice, file.path(out_dir, "R9_C3_4_v3B_fast_Ripley_random_labeling_perslice.csv"))
fwrite(summary_all, file.path(out_dir, "R9_C3_4_v3B_fast_Ripley_random_labeling_summary.csv"))

cat("\n== v3-B fast main gate summary:\n")
print(summary_all[order(k, sensitivity, predictor)])
cat("\n[STOP C3-4 v3-B fast] Evaluate vascular vs neuron_control gate under random-labeling null. CSR was not used.\n")
