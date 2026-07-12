#!/usr/bin/env Rscript
# =============================================================================
# R9 batch 1 | A2 D02 gate review
#
# Goal:
#   Decide whether unbiased domain D02 has a real promotion path.
#
# Outputs:
#   1) Per-slice D02 composition and within-slice enrichment vs non-D02.
#   2) D02 as an independent binary region in v3-B-style crossK / nearest-distance
#      random-labeling null, with neuron_control as the control axis.
#
# Boundaries:
#   - Visium spots are not single cells.
#   - Per-slice statistics are the unit; no pooled spot claim.
#   - Random-labeling null fixes observed spot positions and point counts.
#   - CSR is not used.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

set.seed(1)

base_dir <- getwd()
label_path <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A2_unbiased_domain/A2_domain_labels_per_spot.csv")
score_path <- file.path(base_dir, "tables/R9_stage3_signature_scores/R9_stage3_IVY_hypoxia_scores_per_spot.csv")
threshold_path <- file.path(base_dir, "tables/C3_niche_preflight/R9_C3_1_neighborhood_niche_scores.csv")
out_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A2_D02_gate_review")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

domain_focus <- "D02"
feature_cols <- c("MES-lineage", "MES-V", "vascular", "neuron_control")
n_perm <- 999L
k_values <- c(6L, 12L)
top_fracs <- c(top05 = 0.05, top10 = 0.10, top20 = 0.20)
eps <- 1e-6

labels <- fread(label_path)
scores <- fread(score_path)
thresholds <- fread(threshold_path)
thresholds <- unique(thresholds[, .(slice, k, distance_threshold)])
thresholds <- thresholds[, .(distance_threshold = median(distance_threshold, na.rm = TRUE)), by = .(slice, k)]

required_labels <- c("spot_id", "slice", "image", "x", "y", "domain_id")
required_scores <- c("spot_id", "slice", feature_cols)
missing_labels <- setdiff(required_labels, names(labels))
missing_scores <- setdiff(required_scores, names(scores))
if (length(missing_labels)) stop("missing A2 label columns: ", paste(missing_labels, collapse = ", "))
if (length(missing_scores)) stop("missing score columns: ", paste(missing_scores, collapse = ", "))

dt <- merge(
  labels[, ..required_labels],
  scores[, ..required_scores],
  by = c("spot_id", "slice"),
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(dt[, ..feature_cols])) stop("NA detected in required feature columns after merge.")
dt[, D02_region := domain_id == domain_focus]

wilcox_greater <- function(x, y) {
  if (length(x) < 3 || length(y) < 3) return(NA_real_)
  suppressWarnings(wilcox.test(x, y, alternative = "greater", exact = FALSE)$p.value)
}

slice_rows <- list()
ri <- 0L
for (sl in sort(unique(dt$slice))) {
  sd <- dt[slice == sl]
  in_region <- sd[D02_region == TRUE]
  out_region <- sd[D02_region == FALSE]
  for (feat in feature_cols) {
    ri <- ri + 1L
    mean_region <- mean(in_region[[feat]], na.rm = TRUE)
    mean_bg <- mean(out_region[[feat]], na.rm = TRUE)
    slice_rows[[ri]] <- data.table(
      slice = sl,
      domain_id = domain_focus,
      feature = feat,
      n_region = nrow(in_region),
      n_background = nrow(out_region),
      slice_spots = nrow(sd),
      region_fraction = nrow(in_region) / nrow(sd),
      mean_region = mean_region,
      mean_background_nonD02 = mean_bg,
      delta_mean = mean_region - mean_bg,
      fold_enrichment = (mean_region + eps) / (mean_bg + eps),
      log2FC = log2((mean_region + eps) / (mean_bg + eps)),
      p_wilcox_greater = wilcox_greater(in_region[[feat]], out_region[[feat]])
    )
  }
}
enrich <- rbindlist(slice_rows, use.names = TRUE)
enrich[, fdr_bh_by_feature := p.adjust(p_wilcox_greater, method = "BH"), by = feature]
enrich[, fdr_bh_all_tests := p.adjust(p_wilcox_greater, method = "BH")]
enrich[, pass_positive_fdr_feature := log2FC > 0 & fdr_bh_by_feature < 0.05]

wide_gate <- dcast(
  enrich,
  slice + n_region + slice_spots + region_fraction ~ feature,
  value.var = c("mean_region", "mean_background_nonD02", "log2FC", "fdr_bh_by_feature", "pass_positive_fdr_feature")
)
slice_gate <- enrich[, .(
  pass_MESV_vascular = all(pass_positive_fdr_feature[feature %in% c("MES-V", "vascular")]),
  pass_MESlineage_MESV_vascular = all(pass_positive_fdr_feature[feature %in% c("MES-lineage", "MES-V", "vascular")])
), by = slice]
wide_gate <- merge(wide_gate, slice_gate, by = "slice", all.x = TRUE, sort = FALSE)

gate_summary <- enrich[, .(
  n_slices = .N,
  median_region_mean = median(mean_region, na.rm = TRUE),
  median_background = median(mean_background_nonD02, na.rm = TRUE),
  median_log2FC = median(log2FC, na.rm = TRUE),
  n_slice_log2FC_positive = sum(log2FC > 0, na.rm = TRUE),
  n_slice_fdr_lt_0p05 = sum(fdr_bh_by_feature < 0.05, na.rm = TRUE),
  n_slice_positive_and_fdr = sum(pass_positive_fdr_feature, na.rm = TRUE)
), by = feature]

gate_combined <- data.table(
  gate = c("D02_present", "D02_MESV_and_vascular_positive_FDR", "D02_MESlineage_MESV_vascular_positive_FDR"),
  n_slices = c(
    uniqueN(wide_gate[n_region > 0]$slice),
    sum(wide_gate$pass_MESV_vascular, na.rm = TRUE),
    sum(wide_gate$pass_MESlineage_MESV_vascular, na.rm = TRUE)
  ),
  denominator = uniqueN(dt$slice)
)

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

wilcox_one_sample <- function(x, alternative = "greater") {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, alternative = alternative, exact = FALSE)$p.value)
}

summarize_prox <- function(x) {
  x[, .(
    n_slices = .N,
    median_delta_crossK = median(delta_crossK, na.rm = TRUE),
    q25_delta_crossK = as.numeric(quantile(delta_crossK, 0.25, na.rm = TRUE)),
    q75_delta_crossK = as.numeric(quantile(delta_crossK, 0.75, na.rm = TRUE)),
    n_crossK_positive = sum(delta_crossK > 0, na.rm = TRUE),
    sign_p_crossK_greater = sign_p(delta_crossK),
    wilcox_p_crossK_greater = wilcox_one_sample(delta_crossK),
    median_p_crossK = median(p_emp_crossK_greater, na.rm = TRUE),
    n_slice_crossK_p_lt_0p05 = sum(p_emp_crossK_greater < 0.05, na.rm = TRUE),
    median_delta_nn_closer = median(delta_nn_closer, na.rm = TRUE),
    q25_delta_nn_closer = as.numeric(quantile(delta_nn_closer, 0.25, na.rm = TRUE)),
    q75_delta_nn_closer = as.numeric(quantile(delta_nn_closer, 0.75, na.rm = TRUE)),
    n_nn_closer_positive = sum(delta_nn_closer > 0, na.rm = TRUE),
    sign_p_nn_greater = sign_p(delta_nn_closer),
    wilcox_p_nn_greater = wilcox_one_sample(delta_nn_closer),
    median_p_nn = median(p_emp_nn_less, na.rm = TRUE),
    n_slice_nn_p_lt_0p05 = sum(p_emp_nn_less < 0.05, na.rm = TRUE),
    count_preserved_all = all(count_preserved_all_perm),
    median_n_response = median(n_response_points),
    median_n_predictor = median(n_predictor_points)
  ), by = .(k, top_rule, top_fraction, response, predictor)]
}

prox_rows <- list()
null_rows <- list()
pi <- ni <- 0L

message("== A2 D02 v3-B region proximity: random-labeling null; CSR is not used.")
for (sl in sort(unique(dt$slice))) {
  sd <- dt[slice == sl]
  n <- nrow(sd)
  D <- as.matrix(dist(as.matrix(sd[, .(x, y)])))
  diag(D) <- 0
  D02_idx <- which(sd$D02_region)

  for (kk in k_values) {
    radius <- thresholds[slice == sl & k == kk, distance_threshold][1]
    if (!is.finite(radius)) stop("missing radius for slice=", sl, " k=", kk)

    for (frac_name in names(top_fracs)) {
      frac <- top_fracs[[frac_name]]
      hot <- list(
        `MES-lineage` = which(make_hot(sd[["MES-lineage"]], frac)),
        `MES-V` = which(make_hot(sd[["MES-V"]], frac)),
        vascular = which(make_hot(sd[["vascular"]], frac)),
        neuron_control = which(make_hot(sd[["neuron_control"]], frac)),
        D02_region = D02_idx
      )

      contrasts <- rbindlist(list(
        data.table(response = "D02_region", predictor = c("MES-lineage", "MES-V", "vascular", "neuron_control")),
        data.table(response = c("MES-lineage", "MES-V", "vascular", "neuron_control"), predictor = "D02_region")
      ))

      null_cache <- new.env(parent = emptyenv())
      for (ii in seq_len(nrow(contrasts))) {
        resp <- contrasts$response[ii]
        pred <- contrasts$predictor[ii]
        A_idx <- hot[[resp]]
        B_idx <- hot[[pred]]
        nA <- length(A_idx)
        nB <- length(B_idx)
        cache_key <- paste(nA, nB, sep = "_")
        if (!exists(cache_key, envir = null_cache, inherits = FALSE)) {
          assign(cache_key, make_null(D, n, nA, nB, radius), envir = null_cache)
        }
        nul <- get(cache_key, envir = null_cache, inherits = FALSE)
        obs <- calc_stats(D, A_idx, B_idx, radius)

        pi <- pi + 1L
        prox_rows[[pi]] <- data.table(
          slice = sl,
          k = kk,
          top_rule = frac_name,
          top_fraction = frac,
          response = resp,
          predictor = pred,
          n_spots = n,
          radius = radius,
          n_response_points = nA,
          n_predictor_points = nB,
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

      for (key in ls(null_cache)) {
        nul <- get(key, envir = null_cache)
        parts <- as.integer(strsplit(key, "_", fixed = TRUE)[[1]])
        ni <- ni + 1L
        null_rows[[ni]] <- data.table(
          slice = sl,
          k = kk,
          top_rule = frac_name,
          top_fraction = frac,
          n_response_points = parts[1],
          n_predictor_points = parts[2],
          n_perm = n_perm,
          count_preserved_all_perm = all(nul$count_ok)
        )
      }
    }
  }
  message("  done slice=", sl, " n=", n, " D02=", length(D02_idx))
}

prox <- rbindlist(prox_rows, use.names = TRUE)
prox_summary <- summarize_prox(prox)
prox_summary[, sensitivity := "all_slices"]
prox_summary_no265 <- summarize_prox(prox[slice != "#UKF265_T_ST"])
prox_summary_no265[, sensitivity := "exclude_UKF265"]
prox_summary <- rbindlist(list(prox_summary, prox_summary_no265), use.names = TRUE)
prox_summary[, fdr_bh_crossK_by_topk := p.adjust(wilcox_p_crossK_greater, method = "BH"), by = .(k, top_rule, sensitivity)]
prox_summary[, fdr_bh_nn_by_topk := p.adjust(wilcox_p_nn_greater, method = "BH"), by = .(k, top_rule, sensitivity)]
null_audit <- unique(rbindlist(null_rows, use.names = TRUE))

fwrite(enrich, file.path(out_dir, "A2_D02_perslice_enrichment_vs_nonD02.csv"))
fwrite(wide_gate, file.path(out_dir, "A2_D02_perslice_gate_wide.csv"))
fwrite(gate_summary, file.path(out_dir, "A2_D02_enrichment_summary.csv"))
fwrite(gate_combined, file.path(out_dir, "A2_D02_gate_counts.csv"))
fwrite(prox, file.path(out_dir, "A2_D02_region_v3B_proximity_perslice.csv"))
fwrite(prox_summary, file.path(out_dir, "A2_D02_region_v3B_proximity_summary.csv"))
fwrite(null_audit, file.path(out_dir, "A2_D02_region_v3B_null_audit.csv"))

params <- data.table(
  parameter = c("domain_focus", "features", "background", "n_perm", "k_values", "top_fracs", "null", "boundary"),
  value = c(
    domain_focus,
    paste(feature_cols, collapse = ";"),
    "within-slice non-D02 spots",
    as.character(n_perm),
    paste(k_values, collapse = ";"),
    paste(names(top_fracs), top_fracs, sep = "=", collapse = ";"),
    "fixed observed spot positions and point counts; labels only are permuted; CSR not used",
    "A2 gate review only; no figure promoted automatically"
  )
)
fwrite(params, file.path(out_dir, "A2_D02_gate_review_parameters.csv"))

cat("\n== D02 gate counts:\n")
print(gate_combined)
cat("\n== D02 enrichment summary:\n")
print(gate_summary)
cat("\n== Focus proximity top10 all_slices:\n")
print(prox_summary[sensitivity == "all_slices" & top_rule == "top10"][order(k, response, predictor)])
cat("\n[STOP A2 D02 gate review] Review per-slice recurrence, enrichment/FDR, and neuron-controlled region proximity before any figure decision.\n")
