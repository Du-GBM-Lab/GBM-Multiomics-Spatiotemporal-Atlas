#!/usr/bin/env Rscript
# =============================================================================
# R9-5 | PLAUR / PLAU positive-only top10 sensitivity
# Purpose:
#   Single-gene top10 by SCT can include many zero-expression spots when the
#   raw top10 cutoff is zero. This sensitivity defines gene-high spots only
#   among raw-detected spots. If detected spots are fewer than the nominal
#   top10 count, all detected spots are used.
#
# Scope:
#   responses  = PLAUR, PLAU
#   predictors = MES-V, vascular, neuron_control
#   top cut    = top10 only
#   distance   = k=6, k=12
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })
set.seed(1)

base_dir <- getwd()
score_path <- file.path(base_dir, "tables/R9_5_spatial_gene_coenrichment/R9_5_gene_spot_scores_SCT_and_raw.csv")
threshold_path <- file.path(base_dir, "tables/C3_niche_preflight/R9_C3_1_neighborhood_niche_scores.csv")
out_dir <- file.path(base_dir, "tables/R9_5_spatial_gene_coenrichment")

n_perm <- 999L
k_values <- c(6L, 12L)
frac <- 0.10
responses <- c("PLAUR", "PLAU")
predictors <- c("MES-V", "vascular", "neuron_control")

score_dt <- fread(score_path)
thresholds <- fread(threshold_path)
thresholds <- unique(thresholds[, .(slice, k, distance_threshold)])
thresholds <- thresholds[, .(distance_threshold = median(distance_threshold, na.rm = TRUE)), by = .(slice, k)]

make_hot <- function(x, frac) x >= as.numeric(quantile(x, 1 - frac, na.rm = TRUE))

make_gene_positive_hot <- function(dt, gene, frac) {
  score <- dt[[paste0(gene, "_SCT_data")]]
  detected <- dt[[paste0(gene, "_detected")]]
  pos <- which(detected)
  nominal_n <- ceiling(frac * nrow(dt))
  if (!length(pos)) return(integer(0))
  n_take <- min(nominal_n, length(pos))
  pos[order(score[pos], decreasing = TRUE)[seq_len(n_take)]]
}

calc_stats <- function(D, A_idx, B_idx, radius) {
  DB <- D[A_idx, B_idx, drop = FALSE]
  c(crossK_count = mean(rowSums(DB <= radius)),
    mean_nearest_distance = mean(apply(DB, 1, min)))
}

emp_p <- function(obs, nul, side = c("greater", "less")) {
  side <- match.arg(side)
  nul <- nul[is.finite(nul)]
  if (!is.finite(obs) || !length(nul)) return(NA_real_)
  if (side == "greater") return((1 + sum(nul >= obs)) / (length(nul) + 1))
  (1 + sum(nul <= obs)) / (length(nul) + 1)
}

make_null <- function(D, n, nA, nB, radius) {
  nullK <- numeric(n_perm); nullNN <- numeric(n_perm); count_ok <- logical(n_perm)
  for (p in seq_len(n_perm)) {
    Ai <- sample.int(n, nA); Bi <- sample.int(n, nB)
    count_ok[p] <- length(Ai) == nA && length(Bi) == nB
    st <- calc_stats(D, Ai, Bi, radius)
    nullK[p] <- st[["crossK_count"]]; nullNN[p] <- st[["mean_nearest_distance"]]
  }
  list(nullK = nullK, nullNN = nullNN, count_ok = count_ok)
}

sign_p <- function(x) {
  x <- x[is.finite(x) & x != 0]
  if (!length(x)) return(NA_real_)
  binom.test(sum(x > 0), length(x), p = 0.5, alternative = "greater")$p.value
}
wilcox_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, alternative = "greater", exact = FALSE)$p.value)
}

summarize_dt <- function(dt) {
  dt[, .(
    n_slices = .N,
    median_delta_crossK = median(delta_crossK, na.rm = TRUE),
    n_crossK_positive = sum(delta_crossK > 0, na.rm = TRUE),
    median_p_crossK = median(p_emp_crossK_greater, na.rm = TRUE),
    n_slice_crossK_p_lt_0p05 = sum(p_emp_crossK_greater < 0.05, na.rm = TRUE),
    sign_p_crossK_greater = sign_p(delta_crossK),
    wilcox_p_crossK_greater = wilcox_p(delta_crossK),
    median_delta_nn_closer = median(delta_nn_closer, na.rm = TRUE),
    n_nn_closer_positive = sum(delta_nn_closer > 0, na.rm = TRUE),
    median_p_nn = median(p_emp_nn_less, na.rm = TRUE),
    n_slice_nn_p_lt_0p05 = sum(p_emp_nn_less < 0.05, na.rm = TRUE),
    count_preserved_all = all(count_preserved_all_perm),
    median_response_hot = median(n_response_hot),
    median_response_hot_fraction = median(n_response_hot / n_spots),
    median_predictor_hot = median(n_predictor_hot),
    median_perm = median(n_perm)
  ), by = .(k, response, predictor, top_rule)]
}

res <- list(); null_audit <- list(); ri <- ai <- 0L
message("== R9-5 PLAUR/PLAU positive-only top10 sensitivity.")
for (sl in sort(unique(score_dt$slice))) {
  dt <- score_dt[slice == sl]
  n <- nrow(dt)
  D <- as.matrix(dist(as.matrix(dt[, .(x, y)]))); diag(D) <- 0
  for (kk in k_values) {
    radius <- thresholds[slice == sl & k == kk, distance_threshold][1]
    hot_pred <- lapply(predictors, function(pred) which(make_hot(dt[[pred]], frac)))
    names(hot_pred) <- predictors
    null_cache <- new.env(parent = emptyenv())
    for (resp in responses) {
      A_idx <- make_gene_positive_hot(dt, resp, frac)
      for (pred in predictors) {
        B_idx <- hot_pred[[pred]]
        nA <- length(A_idx); nB <- length(B_idx)
        cache_key <- paste(nA, nB, sep = "_")
        if (!exists(cache_key, envir = null_cache, inherits = FALSE)) {
          assign(cache_key, make_null(D, n, nA, nB, radius), envir = null_cache)
        }
        nul <- get(cache_key, envir = null_cache, inherits = FALSE)
        obs <- calc_stats(D, A_idx, B_idx, radius)
        ri <- ri + 1L
        res[[ri]] <- data.table(
          slice = sl, k = kk, response = resp, predictor = pred,
          top_rule = "positive_top10", n_spots = n, radius = radius,
          n_response_hot = nA, n_predictor_hot = nB,
          observed_crossK_count = obs[["crossK_count"]],
          null_crossK_median = median(nul$nullK),
          delta_crossK = obs[["crossK_count"]] - median(nul$nullK),
          p_emp_crossK_greater = emp_p(obs[["crossK_count"]], nul$nullK, "greater"),
          observed_mean_nn_dist = obs[["mean_nearest_distance"]],
          null_mean_nn_dist_median = median(nul$nullNN),
          delta_nn_closer = median(nul$nullNN) - obs[["mean_nearest_distance"]],
          p_emp_nn_less = emp_p(obs[["mean_nearest_distance"]], nul$nullNN, "less"),
          n_perm = n_perm,
          count_preserved_all_perm = all(nul$count_ok)
        )
      }
    }
    for (key in ls(null_cache)) {
      nul <- get(key, envir = null_cache)
      parts <- as.integer(strsplit(key, "_", fixed = TRUE)[[1]])
      ai <- ai + 1L
      null_audit[[ai]] <- data.table(
        slice = sl, k = kk, top_rule = "positive_top10",
        n_response_hot = parts[1], n_predictor_hot = parts[2],
        n_perm = n_perm, count_preserved_all_perm = all(nul$count_ok)
      )
    }
  }
  message("  done slice=", sl, " n=", n)
}

perslice <- rbindlist(res)
summary <- summarize_dt(perslice)
summary[, sensitivity := "positive_detected_top10_all_slices"]
summary_no265 <- summarize_dt(perslice[slice != "#UKF265_T_ST"])
summary_no265[, sensitivity := "positive_detected_top10_exclude_UKF265"]
summary_all <- rbindlist(list(summary, summary_no265), use.names = TRUE)
null_audit <- unique(rbindlist(null_audit))

fwrite(perslice, file.path(out_dir, "R9_5_PLAUR_PLAU_positive_only_top10_perslice.csv"))
fwrite(summary_all, file.path(out_dir, "R9_5_PLAUR_PLAU_positive_only_top10_summary.csv"))
fwrite(null_audit, file.path(out_dir, "R9_5_PLAUR_PLAU_positive_only_top10_null_audit.csv"))

cat("\n== Positive-only top10 summary, all slices:\n")
print(summary_all[sensitivity == "positive_detected_top10_all_slices"][order(response, k, predictor)])
cat("\n== Null audit:\n")
print(null_audit[, .(n_rows=.N, all_preserved=all(count_preserved_all_perm), min_perm=min(n_perm), max_perm=max(n_perm))])
