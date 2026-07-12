# R9 | server A-pre: install/check spacexr API and reference sanity only
# Server project root:
# /home/data/t010639/projects/GBM_R9_spatial_RCTD
# This script does not run RCTD/deconvolution.

root <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
dir.create(file.path(root, "logs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(root, "tables"), showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("spacexr", quietly = TRUE)) {
  cat("== spacexr not installed; trying classic dmcable/spacexr...\n")
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  options(timeout = 6e8)
  try(devtools::install_github("dmcable/spacexr", upgrade = "never"))
}

if (!requireNamespace("spacexr", quietly = TRUE)) {
  cat("== classic spacexr install failed.\n")
  cat("== If switching to Bioc fork, use: BiocManager::install('spacexr')\n")
  quit(save = "no", status = 1)
}

suppressPackageStartupMessages(library(spacexr))
fns <- ls(getNamespace("spacexr"))

classic_full <- all(c("create.RCTD", "run.RCTD", "Reference", "SpatialRNA") %in% fns)
classic_min <- all(c("create.RCTD", "run.RCTD") %in% fns)
bioc_fork <- any(c("createRctd", "runRctd") %in% fns)
branch <- if (classic_min) "CLASSIC" else if (bioc_fork) "BIOC_FORK" else "UNKNOWN"

cat("== spacexr version:", as.character(packageVersion("spacexr")), "\n")
cat("== classic API create.RCTD/run.RCTD/Reference/SpatialRNA:", classic_full, "\n")
cat("== Bioc fork API createRctd/runRctd:", bioc_fork, "\n")
cat("== branch:", branch, "\n")
cat("== RCTD-related functions:\n")
rctd_fns <- grep("RCTD|Rctd|Reference|SpatialRNA|SpatialExperiment", fns, value = TRUE)
print(rctd_fns)
write.csv(data.frame(function_name = rctd_fns),
          file.path(root, "tables", "Apre_spacexr_RCTD_function_list.csv"),
          row.names = FALSE)
write.csv(data.frame(
  package = "spacexr",
  version = as.character(packageVersion("spacexr")),
  branch = branch,
  classic_full = classic_full,
  bioc_fork = bioc_fork
), file.path(root, "tables", "Apre_spacexr_api_branch.csv"), row.names = FALSE)

suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(Matrix)
})

path_full <- file.path(root, "data/reference/GBM.RNA.qc_doubletfinder.infercnv_immune_reference_calls.qs2")
path_mal <- file.path(root, "data/reference/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2")
stopifnot(file.exists(path_full), file.exists(path_mal))

full <- qs2::qs_read(path_full)
mal <- qs2::qs_read(path_mal)
DefaultAssay(full) <- "RNA"
DefaultAssay(mal) <- "RNA"

cat("\n== full reference:", ncol(full), "cells x", nrow(full), "genes\n")
cat("== malignant subtype object:", ncol(mal), "cells x", nrow(mal), "genes\n")
cat("== full anno_ident counts:\n")
print(sort(table(full$anno_ident, useNA = "ifany"), decreasing = TRUE))

mal_bc <- colnames(mal)
ov <- intersect(mal_bc, colnames(full))
cat("\n== malignant barcodes mapped to full:", length(ov), "/", length(mal_bc), "\n")

ct <- as.character(full$anno_ident)
names(ct) <- colnames(full)
ct[ov] <- as.character(mal$subtype_k4[match(ov, mal_bc)])
ct_tab <- sort(table(ct), decreasing = TRUE)
cat("== RCTD provisional celltype counts:\n")
print(ct_tab)
cat("== all cell types >50:", all(ct_tab > 50), "\n")

write.csv(data.frame(celltype = names(ct_tab), n_cells = as.integer(ct_tab)),
          file.path(root, "tables", "Apre_RCTD_reference_celltype_counts.csv"),
          row.names = FALSE)

cnt <- tryCatch(
  GetAssayData(full, assay = "RNA", layer = "counts"),
  error = function(e1) tryCatch(GetAssayData(full, assay = "RNA", slot = "counts"),
                                error = function(e2) NULL)
)
stopifnot(!is.null(cnt))
nUMI <- Matrix::colSums(cnt)
write.csv(data.frame(
  metric = c("n_cells", "n_genes", "nUMI_min", "nUMI_median", "nUMI_mean", "nUMI_max"),
  value = c(ncol(cnt), nrow(cnt), min(nUMI), median(nUMI), mean(nUMI), max(nUMI))
), file.path(root, "tables", "Apre_RCTD_reference_counts_summary.csv"),
row.names = FALSE)

cat("\n[STOP A-pre server] Report spacexr version, branch, function list, reference celltype counts.\n")
