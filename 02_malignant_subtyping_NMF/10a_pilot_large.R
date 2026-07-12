# 05_恶性细胞分亚群与Neftel对照/10a_pilot_large.R
# Worst-case pilot for per-patient NMF: largest patient x max k.
# This is a dry-run and does not save the NMF fit object.

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
  library(peakRAM)
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
  pilot_patient = "Pt15",
  hvg_n = 3000L,
  k = 7L,
  seed = 1L,
  free_memory_mib = 47000,
  main_process_and_os_mib = 8000,
  safety_margin_mib = 5000
)

dir.create(params$tables_dir, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

write_session_info <- function(path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  sink(con)
  print(sessionInfo())
  sink()
}

recommend_workers <- function(peak_delta_mib) {
  if (peak_delta_mib > 5000) {
    2L
  } else if (peak_delta_mib > 3000) {
    3L
  } else if (peak_delta_mib > 1500) {
    4L
  } else if (peak_delta_mib > 500) {
    5L
  } else {
    6L
  }
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
  count(Pt_number, name = "n_cells") |>
  arrange(desc(n_cells))

cat("\nPatient counts sorted descending:\n")
print(as.data.frame(pt_counts))

cells_use <- rownames(md)[as.character(md$Pt_number) == params$pilot_patient]
stopifnot(length(cells_use) > 5000)

cat("\nWorst-case pilot patient:", params$pilot_patient, "with", length(cells_use), "cells\n")

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
stopifnot(length(hvg) == params$hvg_n)

mat <- GetAssayData(sub, assay = "RNA", layer = "data")[hvg, , drop = FALSE]
stopifnot(inherits(mat, "dgCMatrix"))
stopifnot(!any(is.na(mat@x)))
stopifnot(min(mat@x) >= 0)

density_pct <- length(mat@x) / (nrow(mat) * ncol(mat)) * 100
cat(
  "\nPt15 input matrix:", nrow(mat), "genes x", ncol(mat),
  "cells, sparse density", round(density_pct, 2), "%\n"
)

set.seed(42)
prof <- peakRAM::peakRAM(
  fit <- RcppML::nmf(
    mat,
    k = params$k,
    seed = params$seed,
    L1 = c(0, 0),
    maxit = 1000,
    tol = 1e-6,
    verbose = FALSE
  )
)

stopifnot(all(fit$w >= 0))
stopifnot(all(fit$h >= 0))
stopifnot(!any(is.na(fit$w)))
stopifnot(!any(is.na(fit$h)))
stopifnot(ncol(fit$w) == params$k)
stopifnot(nrow(fit$h) == params$k)
stopifnot(ncol(fit$h) == ncol(mat))

elapsed_sec <- as.numeric(prof$Elapsed_Time_sec)
peak_total_mib <- as.numeric(prof$Total_RAM_Used_MiB)
peak_delta_mib <- as.numeric(prof$Peak_RAM_Used_MiB)
rec_workers <- recommend_workers(peak_delta_mib)

cat("\nPt15 k=7 elapsed:", round(elapsed_sec, 2), "sec\n")
cat("Pt15 k=7 peak mem total:", round(peak_total_mib, 1), "MiB\n")
cat("Pt15 k=7 peak mem delta:", round(peak_delta_mib, 1), "MiB\n")

for (workers in c(6L, 5L, 4L, 3L, 2L)) {
  estimated_worker_mem <- workers * peak_delta_mib
  budget_remaining <- params$free_memory_mib -
    params$main_process_and_os_mib -
    estimated_worker_mem
  pass_margin <- budget_remaining >= params$safety_margin_mib
  cat(
    "Memory budget with", workers, "workers:",
    "worker_peak_sum =", round(estimated_worker_mem, 1), "MiB;",
    "remaining_after_OS =", round(budget_remaining, 1), "MiB;",
    "safety_margin_pass =", pass_margin, "\n"
  )
}

cat("Recommended workers by peak delta rule:", rec_workers, "\n")

new_row <- data.frame(
  patient = params$pilot_patient,
  n_cells = ncol(mat),
  n_hvg = length(hvg),
  k = params$k,
  seed = params$seed,
  elapsed_sec = round(elapsed_sec, 2),
  mem_used_before_MB = NA_real_,
  mem_used_after_MB = round(peak_total_mib, 1),
  mem_delta_MB = round(peak_delta_mib, 1),
  sparse_density_pct = round(density_pct, 3),
  recommended_workers = rec_workers,
  est_total_min_singlecore = NA_real_,
  est_wall_min_parallel = NA_real_,
  peak_mem_total_MiB = round(peak_total_mib, 1),
  peak_mem_delta_MiB = round(peak_delta_mib, 1),
  stringsAsFactors = FALSE
)

timing_path <- file.path(params$tables_dir, "10a_pilot_timing.csv")
if (file.exists(timing_path)) {
  existing <- read.csv(timing_path, stringsAsFactors = FALSE)
  all_cols <- union(colnames(existing), colnames(new_row))
  for (col in setdiff(all_cols, colnames(existing))) existing[[col]] <- NA
  for (col in setdiff(all_cols, colnames(new_row))) new_row[[col]] <- NA
  timing <- dplyr::bind_rows(existing[, all_cols], new_row[, all_cols])
} else {
  timing <- new_row
}
write.csv(timing, timing_path, row.names = FALSE)
write_session_info(file.path(params$tables_dir, "10a_pilot_large_session_info.txt"))

cat("\nAppended pilot timing:\n")
print(new_row)

rm(list = ls())
gc()
