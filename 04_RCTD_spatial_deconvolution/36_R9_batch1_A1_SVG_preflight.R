#!/usr/bin/env Rscript
# =============================================================================
# R9 batch 1 | A1 SVG preflight
# Purpose:
#   Audit inputs, package availability, and planned outputs for unbiased spatial
#   variable gene analysis. No SVG statistics are computed here.
#
# STOP:
#   Do not run SPARK-X, Moran's I, or any gene-level spatial test in this script.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
})

base_dir <- getwd()
st_path <- "<DATA_ROOT>/项目/分型/分型代码/0.对象/5.ST_merge.rds"
out_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A1_SVG")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

method_docs <- data.table(
  analysis = "A1_SVG",
  source = c("SPARK official site", "Squidpy Moran example/API"),
  url = c(
    "https://xzhoulab.github.io/SPARK/",
    "https://squidpy.readthedocs.io/en/stable/notebooks/examples/graph/compute_moran.html"
  ),
  local_decision = c(
    "SPARK-X is the primary SVG method because each Visium slice has >3000 spots except one small slice and total spots are 56026.",
    "Moran's I is planned only as orthogonal SVG sensitivity, not as the main R9 proximity statistic."
  )
)
fwrite(method_docs, file.path(out_dir, "A1_SVG_method_docs.csv"))

pkg_audit <- data.table(
  package = c("SPARK", "Seurat", "Matrix", "data.table"),
  available = c(
    requireNamespace("SPARK", quietly = TRUE),
    requireNamespace("Seurat", quietly = TRUE),
    requireNamespace("Matrix", quietly = TRUE),
    requireNamespace("data.table", quietly = TRUE)
  )
)
fwrite(pkg_audit, file.path(out_dir, "A1_SVG_package_audit.csv"))

st <- readRDS(st_path)
meta <- st@meta.data
sample_col <- "orig.ident"
if (!sample_col %in% colnames(meta)) stop("Missing sample column: ", sample_col)

assay_audit <- data.table(
  assay = names(st@assays),
  class = vapply(st@assays, function(x) class(x)[1], character(1)),
  n_features = vapply(st@assays, nrow, numeric(1)),
  n_spots = vapply(st@assays, ncol, numeric(1))
)
fwrite(assay_audit, file.path(out_dir, "A1_SVG_assay_audit.csv"))

slice_audit <- as.data.table(table(meta[[sample_col]]))
setnames(slice_audit, c("slice", "n_spots"))
slice_audit[, planned_status := fifelse(n_spots >= 3000, "primary_candidate", "low_spot_count_sensitivity_only")]
fwrite(slice_audit, file.path(out_dir, "A1_SVG_slice_spot_audit.csv"))

spatial_layers <- grep("^counts", Layers(st[["Spatial"]]), value = TRUE)
input_audit <- data.table(
  item = c(
    "spatial_object", "sample_column", "n_slices", "n_spots",
    "spatial_assay_class", "sct_assay_class", "n_spatial_count_layers",
    "planned_expression_input", "planned_statistical_unit", "stop_boundary"
  ),
  value = c(
    st_path, sample_col, length(unique(meta[[sample_col]])), ncol(st),
    class(st[["Spatial"]])[1], class(st[["SCT"]])[1], length(spatial_layers),
    "Spatial raw counts for SPARK-X; SCT data only for optional ranking/visualization after STOP",
    "slice-level summaries; no pooled-spot main claim",
    "No SVG computation in preflight"
  )
)
fwrite(input_audit, file.path(out_dir, "A1_SVG_input_audit.csv"))

planned_outputs <- data.table(
  output = c(
    "A1_SVG_sparkx_per_slice.csv",
    "A1_SVG_sparkx_consensus.csv",
    "A1_SVG_moran_sensitivity.csv",
    "A1_SVG_signature_overlap.csv",
    "A1_SVG_parameters.csv"
  ),
  role = c(
    "future SPARK-X per-slice SVG table",
    "future cross-slice SVG recurrence summary",
    "future Moran's I sensitivity table",
    "future overlap with MES/vascular/FOSL1/IVY sets",
    "future method/package/threshold provenance"
  )
)
fwrite(planned_outputs, file.path(out_dir, "A1_SVG_planned_outputs.csv"))

cat("\n[STOP A1 SVG preflight] Inputs and package state audited. No SVG statistics were run.\n")
