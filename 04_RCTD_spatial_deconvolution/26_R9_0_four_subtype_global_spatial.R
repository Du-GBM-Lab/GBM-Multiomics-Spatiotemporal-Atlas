#!/usr/bin/env Rscript
# =============================================================================
# R9-0 | Four-subtype global spatial analysis
# Stage 2 from the spatial spec:
#   1) Spatial separability / non-random structure of four subtype weights.
#      Metric: per-slice Moran's I with permutation null.
#   2) Four subtype x niche spatial preference matrix.
#      Metric: v3-B Ripley-style crossK + nearest-distance with random-labeling.
#
# No figures. No IVY/hypoxia. No region-restricted correlation.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
})

set.seed(1)

base_dir <- getwd()
weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
score_path <- file.path(base_dir, "tables/C3_niche_preflight/R9_C3_1_neighborhood_niche_scores.csv")
out_dir <- file.path(base_dir, "tables/R9_0_four_subtype_global")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_perm <- 999L
k_values <- c(6L, 12L)
top_frac <- 0.10
top_rule <- "top10"

subtypes <- c("NPC-P", "OPC-M", "MES-V", "MES-I")
predictors <- c("vascular", "myeloid", "microglia", "macrophage_monocyte", "neuron_control")

weights <- as.data.table(qs2::qs_read(weights_path))
scores <- fread(score_path)

feat <- data.table(
  spot_id = weights$spot_id,
  slice = weights$slice,
  x = weights$x,
  y = weights$y,
  `NPC-P` = weights$Subtype1,
  `OPC-M` = weights$Subtype2,
  `MES-V` = weights$Subtype3,
  `MES-I` = weights$Subtype4,
  vascular = weights$Endothelial + weights[["Mural cells"]],
  myeloid = weights$Macrophages + weights$Microglial + weights$Monocytes,
  microglia = weights$Microglial,
  macrophage_monocyte = weights$Macrophages + weights$Monocytes,
  neuron_control = weights$Neurons
)

thresholds <- unique(scores[, .(slice, k, distance_threshold)])
thresholds <- thresholds[, .(distance_threshold = median(distance_threshold, na.rm = TRUE)),
                         by = .(slice, k)]

make_hot <- function(x, frac) {
  x >= as.numeric(quantile(x, 1 - frac, na.rm = TRUE))
}

emp_p <- function(obs, nul, side = c("greater", "less")) {
  side <- match.arg(side)
  nul <- nul[is.finite(nul)]
  if (!is.finite(obs) || !length(nul)) return(NA_real_)
  if (side == "greater") return((1 + sum(nul >= obs)) / (length(nul) + 1))
  (1 + sum(nul <= obs)) / (length(nul) + 1)
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

## ---- Moran's I -------------------------------------------------------------
calc_moran <- function(x, ii, jj, n) {
  z <- x - mean(x, na.rm = TRUE)
  denom <- sum(z^2, na.rm = TRUE)
  if (!is.finite(denom) || denom <= 0 || !length(ii)) return(NA_real_)
  s0 <- length(ii)
  (n / s0) * sum(z[ii] * z[jj], na.rm = TRUE) / denom
}

message("== R9-0 stage 2A: Moran's I spatial structure.")
moran_res <- list()
mi <- 0L

for (sl in sort(unique(feat$slice))) {
  dt <- feat[slice == sl]
  n <- nrow(dt)
  D <- as.matrix(dist(as.matrix(dt[, .(x, y)])))
  diag(D) <- Inf

  for (kk in k_values) {
    radius <- thresholds[slice == sl & k == kk, distance_threshold][1]
    idx <- which(D <= radius, arr.ind = TRUE)
    ii <- idx[, 1]
    jj <- idx[, 2]
    n_edges <- length(ii)
    realized_k_mean <- n_edges / n

    for (st in subtypes) {
      x <- dt[[st]]
      obs <- calc_moran(x, ii, jj, n)
      nul <- numeric(n_perm)
      for (p in seq_len(n_perm)) {
        nul[p] <- calc_moran(sample(x), ii, jj, n)
      }
      mi <- mi + 1L
      moran_res[[mi]] <- data.table(
        slice = sl,
        k = kk,
        radius = radius,
        subtype = st,
        n_spots = n,
        n_edges = n_edges,
        realized_k_mean = realized_k_mean,
        moran_I = obs,
        p_emp_greater = emp_p(obs, nul, "greater"),
        null_median = median(nul, na.rm = TRUE),
        null_q95 = as.numeric(quantile(nul, 0.95, na.rm = TRUE)),
        n_perm = n_perm
      )
    }
  }
  message("  Moran done slice=", sl, " n=", n)
}

moran <- rbindlist(moran_res, use.names = TRUE)
moran_summary <- moran[, .(
  n_slices = .N,
  median_moran_I = median(moran_I, na.rm = TRUE),
  q25_moran_I = as.numeric(quantile(moran_I, 0.25, na.rm = TRUE)),
  q75_moran_I = as.numeric(quantile(moran_I, 0.75, na.rm = TRUE)),
  n_positive = sum(moran_I > 0, na.rm = TRUE),
  n_emp_p_lt_0p05 = sum(p_emp_greater < 0.05, na.rm = TRUE),
  median_emp_p = median(p_emp_greater, na.rm = TRUE),
  sign_p_positive = sign_p(moran_I),
  wilcox_p_greater = wilcox_p(moran_I)
), by = .(k, subtype)]

fwrite(moran, file.path(out_dir, "R9_0_four_subtype_MoranI_perslice.csv"))
fwrite(moran_summary, file.path(out_dir, "R9_0_four_subtype_MoranI_summary.csv"))

## ---- Ripley-style preference matrix ----------------------------------------
calc_stats <- function(D, A_idx, B_idx, radius) {
  DB <- D[A_idx, B_idx, drop = FALSE]
  c(
    crossK_count = mean(rowSums(DB <= radius)),
    mean_nearest_distance = mean(apply(DB, 1, min))
  )
}

make_null <- function(D, n, nA, nB, radius) {
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
  list(nullK = nullK, nullNN = nullNN, count_ok = count_ok)
}

summarize_pref <- function(dt) {
  dt[, .(
    n_slices = .N,
    median_delta_crossK = median(delta_crossK, na.rm = TRUE),
    q25_delta_crossK = as.numeric(quantile(delta_crossK, 0.25, na.rm = TRUE)),
    q75_delta_crossK = as.numeric(quantile(delta_crossK, 0.75, na.rm = TRUE)),
    n_crossK_positive = sum(delta_crossK > 0, na.rm = TRUE),
    sign_p_crossK_greater = sign_p(delta_crossK),
    wilcox_p_crossK_greater = wilcox_p(delta_crossK),
    median_p_crossK = median(p_emp_crossK_greater, na.rm = TRUE),
    n_slice_crossK_p_lt_0p05 = sum(p_emp_crossK_greater < 0.05, na.rm = TRUE),
    median_delta_nn_closer = median(delta_nn_closer, na.rm = TRUE),
    q25_delta_nn_closer = as.numeric(quantile(delta_nn_closer, 0.25, na.rm = TRUE)),
    q75_delta_nn_closer = as.numeric(quantile(delta_nn_closer, 0.75, na.rm = TRUE)),
    n_nn_closer_positive = sum(delta_nn_closer > 0, na.rm = TRUE),
    sign_p_nn_greater = sign_p(delta_nn_closer),
    wilcox_p_nn_greater = wilcox_p(delta_nn_closer),
    median_p_nn = median(p_emp_nn_less, na.rm = TRUE),
    n_slice_nn_p_lt_0p05 = sum(p_emp_nn_less < 0.05, na.rm = TRUE),
    count_preserved_all = all(count_preserved_all_perm),
    median_response_hot = median(n_response_hot),
    median_predictor_hot = median(n_predictor_hot),
    median_perm = median(n_perm)
  ), by = .(k, subtype, predictor, top_rule, top_fraction)]
}

message("== R9-0 stage 2B: four subtype x niche preference matrix.")
pref_res <- list()
null_audit <- list()
ri <- ai <- 0L

for (sl in sort(unique(feat$slice))) {
  dt <- feat[slice == sl]
  n <- nrow(dt)
  D <- as.matrix(dist(as.matrix(dt[, .(x, y)])))
  diag(D) <- 0

  for (kk in k_values) {
    radius <- thresholds[slice == sl & k == kk, distance_threshold][1]
    hot_idx <- list()
    for (nm in c(subtypes, predictors)) {
      hot_idx[[nm]] <- which(make_hot(dt[[nm]], top_frac))
    }

    null_cache <- new.env(parent = emptyenv())
    for (st in subtypes) {
      A_idx <- hot_idx[[st]]
      for (pred in predictors) {
        B_idx <- hot_idx[[pred]]
        nA <- length(A_idx)
        nB <- length(B_idx)
        cache_key <- paste(nA, nB, sep = "_")
        if (!exists(cache_key, envir = null_cache, inherits = FALSE)) {
          assign(cache_key, make_null(D, n, nA, nB, radius), envir = null_cache)
        }
        nul <- get(cache_key, envir = null_cache, inherits = FALSE)
        obs <- calc_stats(D, A_idx, B_idx, radius)

        ri <- ri + 1L
        pref_res[[ri]] <- data.table(
          slice = sl,
          k = kk,
          subtype = st,
          predictor = pred,
          top_rule = top_rule,
          top_fraction = top_frac,
          n_spots = n,
          radius = radius,
          n_response_hot = nA,
          n_predictor_hot = nB,
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
        slice = sl,
        k = kk,
        top_rule = top_rule,
        top_fraction = top_frac,
        n_response_hot = parts[1],
        n_predictor_hot = parts[2],
        n_perm = n_perm,
        count_preserved_all_perm = all(nul$count_ok)
      )
    }
  }
  message("  preference done slice=", sl, " n=", n)
}

pref <- rbindlist(pref_res, use.names = TRUE)
pref_null_audit <- unique(rbindlist(null_audit, use.names = TRUE))
pref_summary <- summarize_pref(pref)
pref_summary[, sensitivity := "all_slices"]
pref_summary_no265 <- summarize_pref(pref[slice != "#UKF265_T_ST"])
pref_summary_no265[, sensitivity := "exclude_UKF265"]
pref_summary_all <- rbindlist(list(pref_summary, pref_summary_no265), use.names = TRUE)

gate <- pref_summary_all[sensitivity == "all_slices", .(
  crossK_positive = median_delta_crossK > 0,
  crossK_emp_median_lt_0p05 = median_p_crossK < 0.05,
  crossK_emp_slices = n_slice_crossK_p_lt_0p05,
  nn_positive = median_delta_nn_closer > 0,
  nn_emp_median_lt_0p05 = median_p_nn < 0.05,
  nn_emp_slices = n_slice_nn_p_lt_0p05
), by = .(k, subtype, predictor, top_rule)]

fwrite(pref, file.path(out_dir, "R9_0_four_subtype_niche_preference_perslice.csv"))
fwrite(pref_summary_all, file.path(out_dir, "R9_0_four_subtype_niche_preference_summary.csv"))
fwrite(pref_null_audit, file.path(out_dir, "R9_0_four_subtype_niche_preference_null_audit.csv"))
fwrite(gate, file.path(out_dir, "R9_0_four_subtype_niche_preference_gate.csv"))

cat("\n== Moran's I summary:\n")
print(moran_summary[order(k, subtype)])

cat("\n== Four-subtype x niche preference summary (all slices):\n")
print(pref_summary_all[sensitivity == "all_slices"][order(k, subtype, predictor)])

cat("\n== Preference null count preservation audit:\n")
print(pref_null_audit[, .(
  n_rows = .N,
  all_preserved = all(count_preserved_all_perm),
  min_perm = min(n_perm),
  max_perm = max(n_perm)
)])

cat("\n[STOP R9-0 stage 2] Report subtype spatial structure and subtype x niche preference matrix before any IVY/hypoxia or region-restricted analysis.\n")
