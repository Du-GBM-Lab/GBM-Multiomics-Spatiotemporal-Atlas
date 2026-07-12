#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-3-stat: slice-level summary statistics for C3 continuous association
# Input : C3-2 per-slice observed rho and vector-shift empirical p values.
# Output: cross-slice summary with sign test, Wilcoxon, bootstrap median CI,
#         and Fisher-z mean effect CI. No new spatial statistic and no figures.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
})

set.seed(1)

base_dir <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
in_path <- file.path(base_dir, "tables/R9_C3_2_perslice_vector_shift_association.csv")
out_dir <- file.path(base_dir, "outputs/R9_C3_slice_summary")
tab_dir <- file.path(base_dir, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

n_boot <- 10000L
main_predictors <- c("vascular_niche", "myeloid_niche", "neuron_control_niche")
main_responses <- c("MES-lineage", "MES-V", "MES-I")

dt <- fread(in_path)
stopifnot(all(c("slice", "k", "response", "predictor", "observed_rho",
                "p_emp_two_sided", "n_null_shifts") %in% names(dt)))

boot_ci_median <- function(x, n = n_boot) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(c(lo = NA_real_, hi = NA_real_))
  b <- replicate(n, median(sample(x, length(x), replace = TRUE)))
  as.numeric(quantile(b, c(0.025, 0.975), na.rm = TRUE))
}

sign_test <- function(x, alternative = c("two.sided", "greater", "less")) {
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

fisher_z_summary <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) {
    return(data.table(fisher_z_mean_rho = NA_real_, fisher_z_ci_lo = NA_real_,
                      fisher_z_ci_hi = NA_real_, fisher_z_p_two_sided = NA_real_,
                      fisher_z_p_greater = NA_real_))
  }
  x <- pmin(pmax(x, -0.999999), 0.999999)
  z <- atanh(x)
  tt2 <- t.test(z, mu = 0, alternative = "two.sided")
  ttg <- t.test(z, mu = 0, alternative = "greater")
  data.table(
    fisher_z_mean_rho = tanh(mean(z)),
    fisher_z_ci_lo = tanh(tt2$conf.int[1]),
    fisher_z_ci_hi = tanh(tt2$conf.int[2]),
    fisher_z_p_two_sided = tt2$p.value,
    fisher_z_p_greater = ttg$p.value
  )
}

summarize_one <- function(x) {
  x <- x[is.finite(x)]
  ci <- boot_ci_median(x)
  zsum <- fisher_z_summary(x)
  cbind(data.table(
    n_slices = length(x),
    median_rho = median(x),
    median_rho_ci_lo = ci[1],
    median_rho_ci_hi = ci[2],
    q25_rho = as.numeric(quantile(x, 0.25)),
    q75_rho = as.numeric(quantile(x, 0.75)),
    mean_rho = mean(x),
    sd_rho = sd(x),
    n_positive = sum(x > 0),
    n_negative = sum(x < 0),
    sign_p_two_sided = sign_test(x, "two.sided"),
    sign_p_greater = sign_test(x, "greater"),
    wilcox_p_two_sided = wilcox_p(x, "two.sided"),
    wilcox_p_greater = wilcox_p(x, "greater")
  ), zsum)
}

sum_all <- dt[, summarize_one(observed_rho), by = .(k, response, predictor)]

emp_summary <- dt[, .(
  median_emp_p_two_sided = median(p_emp_two_sided, na.rm = TRUE),
  n_slice_emp_p_lt_0p05 = sum(p_emp_two_sided < 0.05, na.rm = TRUE),
  median_null_shifts = median(n_null_shifts, na.rm = TRUE),
  min_null_shifts = min(n_null_shifts, na.rm = TRUE)
), by = .(k, response, predictor)]

sum_all <- merge(sum_all, emp_summary, by = c("k", "response", "predictor"), all.x = TRUE)
sum_focus <- sum_all[k %in% c(6, 12) & response %in% main_responses & predictor %in% main_predictors]

fwrite(sum_all, file.path(tab_dir, "R9_C3_3stat_slice_level_summary_all_pairs.csv"))
fwrite(sum_focus, file.path(tab_dir, "R9_C3_3stat_slice_level_summary_focus_pairs.csv"))
qs2::qs_save(list(summary_all = sum_all, summary_focus = sum_focus,
                  n_boot = n_boot,
                  note = "Cross-slice summary of C3-2 per-slice rho; no pooled spot inference."),
             file.path(out_dir, "R9_C3_3stat_slice_level_summary.qs2"))

cat("\n== C3-3-stat focus pairs:\n")
print(sum_focus[order(k, response, predictor),
                .(k, response, predictor, n_slices, median_rho, median_rho_ci_lo,
                  median_rho_ci_hi, n_positive, n_negative, sign_p_greater,
                  wilcox_p_greater, fisher_z_mean_rho, fisher_z_ci_lo,
                  fisher_z_ci_hi, fisher_z_p_greater, median_emp_p_two_sided,
                  n_slice_emp_p_lt_0p05, min_null_shifts)])

cat("\n[STOP C3-3-stat] Use this table to annotate C3-3 forest plots.",
    "Main claim should remain effect size + cross-slice consistency,",
    "with neuron_control as negative control and empirical p values reported conservatively.\n")
