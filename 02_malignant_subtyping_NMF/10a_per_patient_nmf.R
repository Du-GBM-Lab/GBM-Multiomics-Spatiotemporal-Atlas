# 05_恶性细胞分亚群与Neftel对照/10a_per_patient_nmf.R
# Per-patient NMF validation, phase 10a.
#
# Architecture:
# - Master loads the v2 Seurat object once, then prepares patient-specific
#   sparse RNA data matrices with top 3000 HVGs.
# - Workers receive only lightweight sparse matrices and run RcppML::nmf().
# - Workers save only W, H, HVGs, cell IDs, and run metadata, not full fit objects.
# - Consensus/cophenetic/chosen-k metrics are computed serially after all NMF
#   runs complete.
#
# Patient-level k selection rule:
# 1. Compute cophenetic[k] for k in 4:7.
# 2. Compute drop[k] = cophenetic[k] - cophenetic[k-1] for k in 5:7.
# 3. If min(cophenetic[4:7]) >= 0.95, choose k=4 by parsimony.
# 4. Else choose the largest k such that all drops up to that k are >= -0.05,
#    i.e. the last k before the first significant drop.
# 5. If cophenetic[4] < 0.85, still report chosen_k but flag low_stability_patient.

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(qs2)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(Matrix)
  library(RcppML)
  library(future)
  library(furrr)
  library(R.utils)
  library(cluster)
})

set.seed(42)

params <- list(
  project_dir = "05_恶性细胞分亚群与Neftel对照",
  input_object = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.v2.qs2"
  ),
  nmf_dir = file.path("05_恶性细胞分亚群与Neftel对照", "outputs", "nmf", "per_patient"),
  tables_dir = file.path("05_恶性细胞分亚群与Neftel对照", "tables"),
  excluded_min_cells = 50L,
  hvg_n = 3000L,
  k_values = 4:7,
  seed_values = 1:5,
  workers = 4L,
  timeout_sec = 120L
)

dir.create(params$nmf_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(params$tables_dir, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

write_session_info <- function(path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  sink(con)
  print(sessionInfo())
  sink()
}

make_connectivity <- function(H) {
  prog <- apply(H, 2, which.max)
  outer(prog, prog, FUN = "==") * 1
}

calc_rss <- function(mat, W, H) {
  recon <- W %*% H
  sum((as.matrix(mat) - recon)^2)
}

dispersion <- function(consensus) {
  mean(4 * (consensus - 0.5)^2)
}

choose_k <- function(df) {
  df <- df |> arrange(k)
  k_vec <- df$k
  coph <- df$cophenetic
  names(coph) <- k_vec

  low_stability <- isTRUE(coph[["4"]] < 0.85)

  if (all(coph >= 0.95, na.rm = TRUE)) {
    return(tibble(
      chosen_k = 4L,
      cophenetic_at_chosen = unname(coph[["4"]]),
      max_drop_at_or_after_chosen = min(diff(coph), na.rm = TRUE),
      low_stability_flag = low_stability,
      justification = "all_k_cophenetic_ge_0.95_parsimony_k4"
    ))
  }

  drops <- diff(coph)
  significant_drop_pos <- which(drops < -0.05)
  if (length(significant_drop_pos) > 0) {
    chosen <- as.integer(names(coph)[significant_drop_pos[1]])
    justification <- paste0("first_significant_drop_after_k", chosen)
  } else {
    chosen <- max(k_vec)
    justification <- "no_significant_drop_choose_largest_k"
  }

  idx <- match(chosen, k_vec)
  drops_after <- if (idx < length(coph)) drops[idx:length(drops)] else numeric(0)
  max_drop_at_or_after <- if (length(drops_after) > 0) min(drops_after, na.rm = TRUE) else NA_real_

  tibble(
    chosen_k = chosen,
    cophenetic_at_chosen = unname(coph[as.character(chosen)]),
    max_drop_at_or_after_chosen = max_drop_at_or_after,
    low_stability_flag = low_stability,
    justification = justification
  )
}

run_one_nmf <- function(task_i, all_tasks, prep_list, out_dir, timeout_sec) {
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  Sys.setenv(OMP_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS = "1")
  suppressPackageStartupMessages({
    library(RcppML)
    library(R.utils)
  })

  task <- all_tasks[task_i, ]
  pt <- as.character(task$patient)
  k <- as.integer(task$k)
  seed_idx <- as.integer(task$seed_idx)
  nmf_seed <- as.integer(task$nmf_seed)
  pt_data <- prep_list[[pt]]

  t0 <- Sys.time()
  fit <- R.utils::withTimeout(
    RcppML::nmf(
      pt_data$mat,
      k = k,
      seed = nmf_seed,
      L1 = c(0, 0),
      maxit = 1000,
      tol = 1e-6,
      verbose = FALSE
    ),
    timeout = timeout_sec,
    onTimeout = "error"
  )
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (any(fit$w < 0) || any(fit$h < 0) || anyNA(fit$w) || anyNA(fit$h)) {
    stop("Invalid NMF fit for ", pt, " k=", k, " seed=", seed_idx, call. = FALSE)
  }

  rss <- calc_rss(pt_data$mat, fit$w, fit$h)
  residual_density <- rss / (nrow(pt_data$mat) * ncol(pt_data$mat))
  out_path <- file.path(out_dir, sprintf("%s_k%d_seed%d_nmf.rds", pt, k, seed_idx))

  out <- list(
    patient = pt,
    k = k,
    seed = seed_idx,
    nmf_seed = nmf_seed,
    W = fit$w,
    H = fit$h,
    rss = rss,
    residual_density = residual_density,
    elapsed_sec = elapsed,
    n_iter = fit$iter,
    tol = fit$tol,
    hvg = pt_data$hvg,
    cell_ids = pt_data$cell_ids
  )
  saveRDS(out, out_path, compress = "xz")

  data.frame(
    patient = pt,
    k = k,
    seed = seed_idx,
    nmf_seed = nmf_seed,
    n_cells = pt_data$n_cells,
    n_hvg = length(pt_data$hvg),
    rss = rss,
    residual_density = residual_density,
    elapsed_sec = round(elapsed, 3),
    n_iter = fit$iter,
    converged = isTRUE(fit$tol < 1e-6),
    rds_path = out_path,
    stringsAsFactors = FALSE
  )
}

start_time <- Sys.time()
msg("Loading v2 object:", params$input_object)
obj <- qs2::qs_read(params$input_object)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}
md <- obj@meta.data
stopifnot("Pt_number" %in% colnames(md))

pt_counts <- md |>
  count(Pt_number, name = "n_cells") |>
  arrange(Pt_number)

excluded <- pt_counts |>
  filter(n_cells < params$excluded_min_cells) |>
  mutate(reason = paste0("n_cells_lt_", params$excluded_min_cells))
eligible <- pt_counts |>
  filter(n_cells >= params$excluded_min_cells) |>
  arrange(Pt_number)
eligible_patients <- as.character(eligible$Pt_number)

write_csv(excluded, file.path(params$tables_dir, "10a_excluded_patients.csv"))

msg("Preparing patient-specific sparse matrices in master.")
prep_list <- list()
hvg_long <- list()
prep_summary <- list()
for (pt in eligible_patients) {
  cells_use <- rownames(md)[as.character(md$Pt_number) == pt]
  sub <- subset(obj, cells = cells_use)
  DefaultAssay(sub) <- "RNA"
  if (inherits(sub[["RNA"]], "Assay5")) {
    sub[["RNA"]] <- JoinLayers(sub[["RNA"]])
  }
  sub <- FindVariableFeatures(
    sub,
    assay = "RNA",
    selection.method = "vst",
    nfeatures = params$hvg_n,
    verbose = FALSE
  )
  hvg <- VariableFeatures(sub)
  mat <- GetAssayData(sub, assay = "RNA", layer = "data")[hvg, , drop = FALSE]
  stopifnot(inherits(mat, "dgCMatrix"))
  stopifnot(!anyNA(mat@x))
  stopifnot(min(mat@x) >= 0)

  prep_list[[pt]] <- list(
    patient = pt,
    mat = mat,
    hvg = hvg,
    cell_ids = colnames(mat),
    n_cells = ncol(mat),
    n_hvg = length(hvg),
    density_pct = length(mat@x) / (nrow(mat) * ncol(mat)) * 100
  )
  hvg_long[[pt]] <- tibble(patient = pt, hvg_gene = hvg, hvg_rank = seq_along(hvg))
  prep_summary[[pt]] <- tibble(
    patient = pt,
    n_cells = ncol(mat),
    n_hvg = length(hvg),
    density_pct = prep_list[[pt]]$density_pct
  )
  rm(sub)
  gc()
  msg("Prepared ", pt, ": ", ncol(mat), " cells, ", length(hvg), " HVGs.")
}

write_csv(bind_rows(hvg_long), file.path(params$tables_dir, "10a_per_patient_hvg.csv"))
write_csv(bind_rows(prep_summary), file.path(params$tables_dir, "10a_patient_matrix_summary.csv"))

rm(obj)
gc()

all_tasks <- expand.grid(
  patient = eligible_patients,
  k = params$k_values,
  seed_idx = params$seed_values,
  stringsAsFactors = FALSE
) |>
  arrange(patient, k, seed_idx) |>
  mutate(
    patient_index = match(patient, eligible_patients),
    nmf_seed = patient_index * 100L + seed_idx
  )

expected_n_runs <- length(eligible_patients) * length(params$k_values) * length(params$seed_values)
stopifnot(nrow(all_tasks) == expected_n_runs)

msg("Starting NMF sweep with ", params$workers, " workers and ", expected_n_runs, " runs.")
options(future.globals.maxSize = 8 * 1024^3)
future::plan(future::multisession, workers = params$workers)

results <- furrr::future_map(
  seq_len(nrow(all_tasks)),
  run_one_nmf,
  all_tasks = all_tasks,
  prep_list = prep_list,
  out_dir = params$nmf_dir,
  timeout_sec = params$timeout_sec,
  .options = furrr::furrr_options(seed = 42),
  .progress = TRUE
)

future::plan(future::sequential)

metrics_df <- bind_rows(results)
write_csv(metrics_df, file.path(params$tables_dir, "10a_per_patient_k_sweep_metrics.csv"))

expected_rds <- file.path(
  params$nmf_dir,
  sprintf("%s_k%d_seed%d_nmf.rds", all_tasks$patient, all_tasks$k, all_tasks$seed_idx)
)
missing_rds <- expected_rds[!file.exists(expected_rds)]
if (length(missing_rds) > 0) {
  stop("Missing NMF rds files: ", paste(head(missing_rds, 5), collapse = ", "), call. = FALSE)
}
stopifnot(nrow(metrics_df) == expected_n_runs)
stopifnot(all(is.finite(metrics_df$rss)))
stopifnot(all(metrics_df$elapsed_sec < params$timeout_sec))

msg("Computing consensus/cophenetic metrics serially.")
consensus_rows <- list()
for (pt in eligible_patients) {
  for (k in params$k_values) {
    files <- file.path(
      params$nmf_dir,
      sprintf("%s_k%d_seed%d_nmf.rds", pt, k, params$seed_values)
    )
    H_list <- lapply(files, function(path) readRDS(path)$H)
    consensus <- matrix(0, nrow = ncol(H_list[[1]]), ncol = ncol(H_list[[1]]))
    for (H in H_list) {
      consensus <- consensus + make_connectivity(H)
    }
    consensus <- consensus / length(H_list)

    d <- as.dist(1 - consensus)
    hc <- hclust(d, method = "average")
    coph <- suppressWarnings(cor(d, cophenetic(hc)))
    disp <- dispersion(consensus)
    prog <- apply(H_list[[1]], 2, which.max)
    sil <- tryCatch(
      mean(cluster::silhouette(prog, d)[, 3]),
      error = function(e) NA_real_
    )

    consensus_rows[[paste(pt, k, sep = "_")]] <- tibble(
      patient = pt,
      k = k,
      cophenetic = coph,
      dispersion = disp,
      mean_silhouette = sil
    )
    rm(H_list, consensus, d, hc)
    gc()
  }
}
consensus_df <- bind_rows(consensus_rows)
write_csv(consensus_df, file.path(params$tables_dir, "10a_per_patient_consensus_metrics.csv"))

chosen_df <- consensus_df |>
  group_by(patient) |>
  group_modify(~ choose_k(.x)) |>
  ungroup()
write_csv(chosen_df, file.path(params$tables_dir, "10a_per_patient_chosen_k.csv"))

runtime_summary <- tibble(
  started_at = format(start_time, "%Y-%m-%d %H:%M:%S"),
  completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  total_wall_min = as.numeric(difftime(Sys.time(), start_time, units = "mins")),
  n_workers = params$workers,
  n_runs = nrow(metrics_df),
  mean_per_run_sec = mean(metrics_df$elapsed_sec),
  median_per_run_sec = median(metrics_df$elapsed_sec),
  max_per_run_sec = max(metrics_df$elapsed_sec),
  n_patients_eligible = length(eligible_patients),
  n_patients_excluded = nrow(excluded)
)
write_csv(runtime_summary, file.path(params$tables_dir, "10a_runtime_summary.csv"))
write_session_info(file.path(params$tables_dir, "10a_session_info.txt"))

cat("\nRuntime summary:\n")
print(runtime_summary)
cat("\nChosen k:\n")
print(chosen_df)
cat("\nLow stability patients:\n")
print(chosen_df |> filter(low_stability_flag))
msg("Done.")
