# 05_恶性细胞分亚群与Neftel对照/10c_permanova_implementation_verification.R
# Verify the custom fast two-group Euclidean PERMANOVA against vegan::adonis2
# on a 2000-cell S3/S4 subsample.

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(readr)
  library(dplyr)
  library(vegan)
})

set.seed(42)

params <- list(
  score_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10c_per_cell_metaprogram_scores.tsv"
  ),
  out_csv = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10c_permanova_implementation_verification.csv"
  ),
  out_session = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10c_permanova_implementation_verification_session_info.txt"
  ),
  score_cols = c("MP01", "MP02", "MP03", "MP04", "MP06"),
  n_subsample = 2000L,
  n_perm = 999L
)

fast_permanova_two_group <- function(score_mat, group, n_perm = 999, seed = 42) {
  group <- factor(group)
  if (nlevels(group) != 2) {
    stop("fast_permanova_two_group requires exactly two groups.", call. = FALSE)
  }
  X <- as.matrix(score_mat)
  n <- nrow(X)
  p <- ncol(X)
  g <- nlevels(group)
  overall <- colMeans(X)
  total_ss <- sum(rowSums((X - matrix(overall, n, p, byrow = TRUE))^2))

  between_ss_for_group <- function(grp) {
    sum(vapply(levels(grp), function(lv) {
      idx <- grp == lv
      n_g <- sum(idx)
      mu_g <- colMeans(X[idx, , drop = FALSE])
      n_g * sum((mu_g - overall)^2)
    }, numeric(1)))
  }

  ss_between <- between_ss_for_group(group)
  ss_within <- total_ss - ss_between
  f_obs <- (ss_between / (g - 1)) / (ss_within / (n - g))
  r2 <- ss_between / total_ss

  set.seed(seed)
  f_perm <- replicate(n_perm, {
    grp_perm <- factor(sample(group), levels = levels(group))
    ss_b <- between_ss_for_group(grp_perm)
    ss_w <- total_ss - ss_b
    (ss_b / (g - 1)) / (ss_w / (n - g))
  })
  p_value <- (sum(f_perm >= f_obs) + 1) / (n_perm + 1)

  list(F = f_obs, R2 = r2, p_value = p_value)
}

scores <- read_tsv(params$score_file, show_col_types = FALSE) |>
  filter(subtype_k4 %in% c("Subtype3", "Subtype4"))

stopifnot(all(params$score_cols %in% colnames(scores)))
stopifnot(nrow(scores) >= params$n_subsample)

sub_idx <- sample(seq_len(nrow(scores)), params$n_subsample)
sub_scores <- scores[sub_idx, , drop = FALSE]
score_mat <- as.matrix(sub_scores[, params$score_cols])
sub_label <- factor(sub_scores$subtype_k4)

set.seed(42)
vegan_result <- vegan::adonis2(
  score_mat ~ sub_label,
  permutations = params$n_perm,
  method = "euclidean"
)

custom_result <- fast_permanova_two_group(
  score_mat,
  sub_label,
  n_perm = params$n_perm,
  seed = 42
)

f_diff <- abs(vegan_result$F[1] - custom_result$F)
r2_diff <- abs(vegan_result$R2[1] - custom_result$R2)
p_diff <- abs(vegan_result$`Pr(>F)`[1] - custom_result$p_value)

verification <- tibble(
  n_cells = params$n_subsample,
  n_features = length(params$score_cols),
  n_perm = params$n_perm,
  vegan_F = vegan_result$F[1],
  custom_F = custom_result$F,
  F_abs_diff = f_diff,
  vegan_R2 = vegan_result$R2[1],
  custom_R2 = custom_result$R2,
  R2_abs_diff = r2_diff,
  vegan_p = vegan_result$`Pr(>F)`[1],
  custom_p = custom_result$p_value,
  p_abs_diff = p_diff,
  F_R2_identical = f_diff < 1e-6 & r2_diff < 1e-6
)

stopifnot(verification$F_R2_identical)
write_csv(verification, params$out_csv)
writeLines(capture.output(sessionInfo()), params$out_session)

print(verification)
