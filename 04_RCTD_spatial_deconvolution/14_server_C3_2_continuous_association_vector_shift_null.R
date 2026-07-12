#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-2: per-slice continuous association with vector-shift spatial null
# Input : C3-1 neighborhood niche score table + C3-2a2 vector-shift candidates.
# Method: for each slice/k/response/predictor, compute observed Spearman rho.
#         Then shift the predictor field by audited physical vectors and remap
#         to real spots to build an empirical spatial null.
# STOP  : report effect sizes and null p-values only. No pooled spot inference,
#         no causal interpretation, no figure decision.
# =============================================================================

suppressPackageStartupMessages({
  library(qs2)
  library(data.table)
})

set.seed(1)

base_dir <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
score_path <- file.path(base_dir, "tables/R9_C3_1_neighborhood_niche_scores.csv")
shift_path <- file.path(base_dir, "tables/R9_C3_2a2_vector_shift_coverage.csv")
out_dir <- file.path(base_dir, "outputs/R9_C3_continuous_association")
tab_dir <- file.path(base_dir, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

min_shift_coverage <- 0.80
min_pairs <- 50L

responses <- c("MES-lineage", "MES-V", "MES-I")
predictors <- c("vascular_niche", "myeloid_niche", "microglia_niche",
                "macrophage_monocyte_niche", "neuron_control_niche")

scores <- fread(score_path)
shifts <- fread(shift_path)
need_cols <- c("spot_id", "slice", "x", "y", "k", responses, predictors)
stopifnot(all(need_cols %in% names(scores)))
stopifnot(all(c("slice", "shift_x", "shift_y", "coverage", "max_match_dist") %in% names(shifts)))

nearest_match <- function(query, ref, max_dist) {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    stop("Package FNN is required for C3-2 vector-shift remapping on the server.")
  }
  nn <- FNN::get.knnx(ref, query, k = 1)
  list(index = nn$nn.index[, 1],
       dist = nn$nn.dist[, 1],
       valid = nn$nn.dist[, 1] <= max_dist)
}

safe_spearman <- function(x, y) {
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

one_combo <- function(dt, sh, response, predictor) {
  xy <- as.matrix(dt[, .(x, y)])
  y <- dt[[response]]
  x <- dt[[predictor]]
  obs <- safe_spearman(y, x)

  null_rho <- numeric(0)
  null_n <- integer(0)
  null_cov <- numeric(0)
  for (ii in seq_len(nrow(sh))) {
    query <- cbind(dt$x + sh$shift_x[ii], dt$y + sh$shift_y[ii])
    mt <- nearest_match(query, xy, sh$max_match_dist[ii])
    if (!any(mt$valid)) next
    yy <- y[mt$valid]
    xx <- x[mt$index[mt$valid]]
    rr <- safe_spearman(yy, xx)
    null_rho <- c(null_rho, rr)
    null_n <- c(null_n, sum(mt$valid))
    null_cov <- c(null_cov, mean(mt$valid))
  }
  null_rho <- null_rho[is.finite(null_rho)]

  data.table(
    response = response,
    predictor = predictor,
    observed_rho = obs,
    n_spots = nrow(dt),
    n_null_shifts = length(null_rho),
    null_rho_mean = if (length(null_rho)) mean(null_rho) else NA_real_,
    null_rho_sd = if (length(null_rho) > 1) sd(null_rho) else NA_real_,
    null_rho_median = if (length(null_rho)) median(null_rho) else NA_real_,
    p_emp_two_sided = emp_p(obs, null_rho, "two.sided"),
    p_emp_greater = emp_p(obs, null_rho, "greater"),
    p_emp_less = emp_p(obs, null_rho, "less"),
    min_shift_coverage_used = if (length(null_cov)) min(null_cov) else NA_real_,
    median_shift_coverage_used = if (length(null_cov)) median(null_cov) else NA_real_,
    median_mapped_spots = if (length(null_n)) median(null_n) else NA_real_
  )
}

message("== C3-2 input scores: ", nrow(scores), " rows; slices=", length(unique(scores$slice)),
        "; k=", paste(sort(unique(scores$k)), collapse = ","))
message("== vector shifts with coverage >= ", min_shift_coverage, " will be used.")

perslice <- list()
idx <- 0L
for (sl in sort(unique(scores$slice))) {
  for (kk in sort(unique(scores$k))) {
    dt <- scores[slice == sl & k == kk]
    sh <- shifts[slice == sl & coverage >= min_shift_coverage]
    if (!nrow(dt) || !nrow(sh)) next
    for (resp in responses) {
      for (pred in predictors) {
        idx <- idx + 1L
        ans <- one_combo(dt, sh, resp, pred)
        ans[, `:=`(slice = sl, k = kk)]
        perslice[[idx]] <- ans
      }
    }
  }
}
perslice <- rbindlist(perslice, use.names = TRUE, fill = TRUE)
setcolorder(perslice, c("slice", "k", "response", "predictor",
                        setdiff(names(perslice), c("slice", "k", "response", "predictor"))))

wilcox_p <- function(x, alternative = "two.sided") {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, alternative = alternative, exact = FALSE)$p.value)
}
sign_p <- function(x) {
  x <- x[is.finite(x) & x != 0]
  if (!length(x)) return(NA_real_)
  binom.test(sum(x > 0), length(x), p = 0.5, alternative = "two.sided")$p.value
}

summary <- perslice[, .(
  n_slices = sum(is.finite(observed_rho)),
  median_observed_rho = median(observed_rho, na.rm = TRUE),
  q25_observed_rho = quantile(observed_rho, 0.25, na.rm = TRUE),
  q75_observed_rho = quantile(observed_rho, 0.75, na.rm = TRUE),
  n_positive = sum(observed_rho > 0, na.rm = TRUE),
  n_negative = sum(observed_rho < 0, na.rm = TRUE),
  sign_test_p_two_sided = sign_p(observed_rho),
  wilcox_p_two_sided = wilcox_p(observed_rho, "two.sided"),
  wilcox_p_greater = wilcox_p(observed_rho, "greater"),
  median_emp_p_two_sided = median(p_emp_two_sided, na.rm = TRUE),
  n_slice_emp_p_lt_0p05 = sum(p_emp_two_sided < 0.05, na.rm = TRUE),
  median_null_shifts = median(n_null_shifts, na.rm = TRUE),
  min_null_shifts = min(n_null_shifts, na.rm = TRUE),
  min_shift_coverage_used = min(min_shift_coverage_used, na.rm = TRUE)
), by = .(k, response, predictor)]

fwrite(perslice, file.path(tab_dir, "R9_C3_2_perslice_vector_shift_association.csv"))
fwrite(summary, file.path(tab_dir, "R9_C3_2_slice_level_association_summary.csv"))
qs2::qs_save(list(perslice = perslice, summary = summary,
                  min_shift_coverage = min_shift_coverage,
                  responses = responses, predictors = predictors,
                  note = "per-slice Spearman plus vector-shift spatial null; no pooled spot inference"),
             file.path(out_dir, "R9_C3_2_continuous_association_vector_shift_null.qs2"))

cat("\n== C3-2 slice-level summary:\n")
print(summary[order(k, response, predictor)])
cat("\n== Main pairs, k=6:\n")
print(summary[k == 6 & response %in% c("MES-lineage", "MES-V", "MES-I") &
                predictor %in% c("vascular_niche", "myeloid_niche", "neuron_control_niche")]
      [order(response, predictor)])

cat("\n[STOP C3-2] Report main k=6 and k=12 sensitivity: median rho, positive/negative slices,",
    "slice-level Wilcoxon/sign tests, and vector-shift empirical p summaries.",
    "Do not write pooled spot significance or causal/physical-contact language.\n")
