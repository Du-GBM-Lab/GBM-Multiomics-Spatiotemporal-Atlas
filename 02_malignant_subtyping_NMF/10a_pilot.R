# 05_恶性细胞分亚群与Neftel对照/10a_pilot.R
# Single-run pilot for per-patient NMF timing and memory planning.
# This is a dry-run: it does not save NMF fit objects.

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(qs2)
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(RcppML)
})

set.seed(42)

params <- list(
  project_dir = "05_恶性细胞分亚群与Neftel对照",
  input_object = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.v2.qs2"
  ),
  tables_dir = file.path("05_恶性细胞分亚群与Neftel对照", "tables"),
  excluded_patient = "Pt1",
  hvg_n = 3000L,
  k = 5L,
  seed = 1L,
  total_runs_full_sweep = 23L * 4L * 5L
)

dir.create(params$tables_dir, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

gc_used_mb <- function(x) sum(as.numeric(x[, 2]))
gc_max_used_mb <- function(x) sum(as.numeric(x[, 6]))

write_session_info <- function(path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  sink(con)
  print(sessionInfo())
  sink()
}

msg("Input object:", params$input_object)
stopifnot(file.exists(params$input_object))

obj <- qs2::qs_read(params$input_object)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

md <- obj@meta.data
stopifnot("Pt_number" %in% colnames(md))

pt_counts <- md |>
  filter(Pt_number != params$excluded_patient) |>
  count(Pt_number, name = "n_cells") |>
  arrange(n_cells)

cat("\nEligible patient counts sorted by malignant cell count:\n")
print(as.data.frame(pt_counts))

pilot_idx <- ceiling(nrow(pt_counts) / 2)
pilot_patient <- as.character(pt_counts$Pt_number[pilot_idx])
n_cells <- pt_counts$n_cells[pilot_idx]

cat("\nPilot patient:", pilot_patient, "with", n_cells, "cells\n")

cells_use <- rownames(md)[md$Pt_number == pilot_patient]
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
if (length(hvg) < params$hvg_n) {
  warning("Fewer HVGs than requested: ", length(hvg), " < ", params$hvg_n)
}
stopifnot(length(hvg) > 0)

mat <- GetAssayData(sub, assay = "RNA", layer = "data")[hvg, , drop = FALSE]
stopifnot(inherits(mat, "dgCMatrix"))
stopifnot(!any(is.na(mat@x)))
stopifnot(min(mat@x) >= 0)

sparse_density_pct <- length(mat@x) / (nrow(mat) * ncol(mat)) * 100
cat(
  "\nInput matrix:", nrow(mat), "genes x", ncol(mat), "cells, sparse density",
  round(sparse_density_pct, 2), "%\n"
)

gc_before <- gc(reset = TRUE)
mem_before_mb <- gc_used_mb(gc_before)
cat("\nMem before NMF:", round(mem_before_mb, 1), "MB\n")

set.seed(42)
t0 <- Sys.time()
fit <- RcppML::nmf(
  mat,
  k = params$k,
  seed = params$seed,
  L1 = c(0, 0),
  maxit = 1000,
  tol = 1e-6,
  verbose = FALSE
)
t1 <- Sys.time()
elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))

gc_after <- gc()
mem_after_mb <- gc_max_used_mb(gc_after)
mem_delta_mb <- mem_after_mb - mem_before_mb

cat("\nNMF elapsed:", round(elapsed_sec, 2), "sec\n")
cat("Peak mem after NMF:", round(mem_after_mb, 1), "MB (delta", round(mem_delta_mb, 1), "MB)\n")

stopifnot(all(fit$w >= 0))
stopifnot(all(fit$h >= 0))
stopifnot(!any(is.na(fit$w)))
stopifnot(!any(is.na(fit$h)))
stopifnot(ncol(fit$w) == params$k)
stopifnot(nrow(fit$h) == params$k)
stopifnot(ncol(fit$h) == ncol(mat))

total_cells_eligible <- sum(pt_counts$n_cells)
runs_per_patient <- 4L * 5L
sec_per_cell_per_run <- elapsed_sec / n_cells
est_total_sec_singlecore <- sec_per_cell_per_run * total_cells_eligible * runs_per_patient

if (mem_delta_mb > 4000) {
  rec_workers <- 2L
} else if (mem_delta_mb > 2000) {
  rec_workers <- 4L
} else {
  rec_workers <- 6L
}

cat(
  "\nEstimated single-core total time:",
  round(est_total_sec_singlecore / 60, 1),
  "min\n"
)
cat(
  "Recommended workers:", rec_workers,
  "-> estimated wall time:",
  round(est_total_sec_singlecore / 60 / rec_workers, 1),
  "min\n"
)

pilot_row <- data.frame(
  patient = pilot_patient,
  n_cells = n_cells,
  n_hvg = length(hvg),
  k = params$k,
  seed = params$seed,
  elapsed_sec = round(elapsed_sec, 2),
  mem_used_before_MB = round(mem_before_mb, 1),
  mem_used_after_MB = round(mem_after_mb, 1),
  mem_delta_MB = round(mem_delta_mb, 1),
  sparse_density_pct = round(sparse_density_pct, 3),
  recommended_workers = rec_workers,
  est_total_min_singlecore = round(est_total_sec_singlecore / 60, 1),
  est_wall_min_parallel = round(est_total_sec_singlecore / 60 / rec_workers, 1),
  stringsAsFactors = FALSE
)

write.csv(
  pilot_row,
  file.path(params$tables_dir, "10a_pilot_timing.csv"),
  row.names = FALSE
)
write_session_info(file.path(params$tables_dir, "10a_pilot_session_info.txt"))

cat("\nPilot output:\n")
print(pilot_row)

rm(list = ls())
gc()
