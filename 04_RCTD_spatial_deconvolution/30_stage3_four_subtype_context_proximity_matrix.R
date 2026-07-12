#!/usr/bin/env Rscript
# =============================================================================
# R9 | Stage 3: four-subtype x context proximity matrix
# Scope:
#   responses  = NPC-P, OPC-M, MES-V, MES-I, MES-lineage
#   predictors = vascular, Buffa hypoxia, CTmvp, CTpan, neuron_control
#   top cut    = top5, top10, top20
#   distance   = k=6, k=12 physical radius
#
# Method:
#   Same v3-B random-labeling null as C3-4 and Stage 3 script 29.
#   Fixed observed spot positions and point counts; class labels are permuted.
#   CSR is not used. No subset correlation is performed.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
})

set.seed(1)

base_dir <- getwd()
weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
score_path <- file.path(base_dir, "tables/R9_stage3_signature_scores/R9_stage3_IVY_hypoxia_scores_per_spot.csv")
threshold_path <- file.path(base_dir, "tables/C3_niche_preflight/R9_C3_1_neighborhood_niche_scores.csv")
out_dir <- file.path(base_dir, "tables/R9_stage3_signature_scores")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_perm <- 999L
k_values <- c(6L, 12L)
top_fracs <- c(top05 = 0.05, top10 = 0.10, top20 = 0.20)
responses <- c("NPC-P", "OPC-M", "MES-V", "MES-I", "MES-lineage")
predictors <- c(
  vascular = "vascular",
  Hypoxia_Buffa = "Hypoxia_Buffa",
  CTmvp = "IVY_CTmvp",
  CTpan = "IVY_CTpan",
  neuron_control = "neuron_control"
)

weights <- as.data.table(qs2::qs_read(weights_path))
score_dt <- fread(score_path)
thresholds <- fread(threshold_path)
thresholds <- unique(thresholds[, .(slice, k, distance_threshold)])
thresholds <- thresholds[, .(distance_threshold = median(distance_threshold, na.rm = TRUE)),
                         by = .(slice, k)]

required_weights <- c("spot_id", "slice", "x", "y", "Subtype1", "Subtype2", "Subtype3", "Subtype4",
                      "Endothelial", "Mural cells", "Neurons")
missing_weights <- setdiff(required_weights, names(weights))
if (length(missing_weights)) {
  stop("missing columns in RCTD weights: ", paste(missing_weights, collapse = ", "))
}

required_scores <- c("spot_id", "slice", "Hypoxia_Buffa", "IVY_CTmvp", "IVY_CTpan")
missing_scores <- setdiff(required_scores, names(score_dt))
if (length(missing_scores)) {
  stop("missing columns in signature score table: ", paste(missing_scores, collapse = ", "))
}

feat <- data.table(
  spot_id = weights$spot_id,
  slice = weights$slice,
  x = weights$x,
  y = weights$y,
  `NPC-P` = weights$Subtype1,
  `OPC-M` = weights$Subtype2,
  `MES-V` = weights$Subtype3,
  `MES-I` = weights$Subtype4,
  `MES-lineage` = weights$Subtype3 + weights$Subtype4,
  vascular = weights$Endothelial + weights[["Mural cells"]],
  neuron_control = weights$Neurons
)

score_keep <- score_dt[, .(spot_id, slice, Hypoxia_Buffa, IVY_CTmvp, IVY_CTpan)]
feat <- merge(feat, score_keep, by = c("spot_id", "slice"), all.x = TRUE, sort = FALSE)

required <- c("spot_id", "slice", "x", "y", responses, unname(predictors))
missing <- setdiff(required, names(feat))
if (length(missing)) stop("missing columns in feature table: ", paste(missing, collapse = ", "))
if (anyNA(feat[, ..required])) stop("NA detected in required feature columns.")

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
  nul <- nul[is.finite(nul)]
  if (!is.finite(obs) || !length(nul)) return(NA_real_)
  if (side == "greater") return((1 + sum(nul >= obs)) / (length(nul) + 1))
  (1 + sum(nul <= obs)) / (length(nul) + 1)
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

summarize_dt <- function(dt) {
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
  ), by = .(k, response, predictor, predictor_col, top_rule, top_fraction)]
}

message("== Four-subtype x context proximity matrix: random-labeling null; CSR is not used.")
res <- list()
null_audit <- list()
ri <- ai <- 0L

for (sl in sort(unique(feat$slice))) {
  dt <- feat[slice == sl]
  n <- nrow(dt)
  D <- as.matrix(dist(as.matrix(dt[, .(x, y)])))
  diag(D) <- 0

  for (kk in k_values) {
    radius <- thresholds[slice == sl & k == kk, distance_threshold][1]

    for (frac_name in names(top_fracs)) {
      frac <- top_fracs[[frac_name]]
      hot_idx <- list()
      for (resp in responses) {
        hot_idx[[resp]] <- which(make_hot(dt[[resp]], frac))
      }
      for (pred in names(predictors)) {
        hot_idx[[pred]] <- which(make_hot(dt[[predictors[[pred]]]], frac))
      }

      null_cache <- new.env(parent = emptyenv())
      for (resp in responses) {
        A_idx <- hot_idx[[resp]]
        for (pred in names(predictors)) {
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
          res[[ri]] <- data.table(
            slice = sl,
            k = kk,
            response = resp,
            predictor = pred,
            predictor_col = predictors[[pred]],
            top_rule = frac_name,
            top_fraction = frac,
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
          top_rule = frac_name,
          top_fraction = frac,
          n_response_hot = parts[1],
          n_predictor_hot = parts[2],
          n_perm = n_perm,
          count_preserved_all_perm = all(nul$count_ok)
        )
      }
    }
  }
  message("  done slice=", sl, " n=", n)
}

perslice <- rbindlist(res, use.names = TRUE)
null_audit <- unique(rbindlist(null_audit, use.names = TRUE))
summary <- summarize_dt(perslice)
summary[, sensitivity := "all_slices"]
summary_no265 <- summarize_dt(perslice[slice != "#UKF265_T_ST"])
summary_no265[, sensitivity := "exclude_UKF265"]
summary_all <- rbindlist(list(summary, summary_no265), use.names = TRUE)

fwrite(perslice, file.path(out_dir, "R9_stage3_four_subtype_context_proximity_perslice.csv"))
fwrite(summary_all, file.path(out_dir, "R9_stage3_four_subtype_context_proximity_summary.csv"))
fwrite(null_audit, file.path(out_dir, "R9_stage3_four_subtype_context_proximity_null_audit.csv"))

cat("\n== Focus: vascular and Buffa columns, top10 all slices:\n")
print(summary_all[
  sensitivity == "all_slices" &
    top_rule == "top10" &
    predictor %in% c("vascular", "Hypoxia_Buffa", "neuron_control")
][order(predictor, k, response)])

cat("\n== Null count preservation audit:\n")
print(null_audit[, .(
  n_rows = .N,
  all_preserved = all(count_preserved_all_perm),
  min_perm = min(n_perm),
  max_perm = max(n_perm)
)])

cat("\n[STOP four-subtype context matrix] Report vascular and Buffa columns with neuron control before interpreting CTmvp/CTpan.\n")
