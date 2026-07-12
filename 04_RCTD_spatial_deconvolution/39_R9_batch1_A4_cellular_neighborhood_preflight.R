#!/usr/bin/env Rscript
# =============================================================================
# R9 batch 1 | A4 cellular neighborhood preflight
# Purpose:
#   Audit the existing per-spot neighborhood composition table and define the
#   future cellular-neighborhood clustering contract. No CN clustering is run.
#
# STOP:
#   Do not run k-means, hierarchical clustering, Leiden, or any CN assignment in
#   this script.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

base_dir <- getwd()
niche_score_path <- file.path(base_dir, "tables/C3_niche_preflight/R9_C3_1_neighborhood_niche_scores.csv")
out_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A4_cellular_neighborhood")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

method_docs <- data.table(
  analysis = "A4_cellular_neighborhood",
  source = c("OSTA neighborhood analysis", "single-cell best practices neighborhood chapter"),
  url = c(
    "https://lmweber.org/OSTA/pages/img-neighborhood-analysis.html",
    "https://www.sc-best-practices.org/spatial/neighborhood.html"
  ),
  local_decision = c(
    "Cellular neighborhoods are recurring local composition patterns; in Visium/RCTD they are spot-neighborhood composition clusters, not cell-cell contact neighborhoods.",
    "R9 CN output is a landscape layer; MES-V/MES-lineage enrichment in a vascular-rich CN can only be a candidate orthogonal support."
  )
)
fwrite(method_docs, file.path(out_dir, "A4_CN_method_docs.csv"))

pkg_audit <- data.table(
  package = c("data.table", "stats", "cluster", "mclust", "igraph"),
  available = c(
    requireNamespace("data.table", quietly = TRUE),
    TRUE,
    requireNamespace("cluster", quietly = TRUE),
    requireNamespace("mclust", quietly = TRUE),
    requireNamespace("igraph", quietly = TRUE)
  )
)
fwrite(pkg_audit, file.path(out_dir, "A4_CN_package_audit.csv"))

dt <- fread(niche_score_path)
required_cols <- c(
  "spot_id", "slice", "image", "x", "y", "k", "realized_neighbors",
  "MES-lineage", "MES-V", "MES-I",
  "vascular_niche", "myeloid_niche", "microglia_niche",
  "macrophage_monocyte_niche", "neuron_control_niche"
)
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols)) stop("Missing required neighborhood columns: ", paste(missing_cols, collapse = ", "))

input_audit <- data.table(
  item = c("niche_score_path", "n_spots", "n_slices", "k_values", "required_columns_present",
           "planned_feature_matrix", "planned_primary_unit", "stop_boundary"),
  value = c(
    niche_score_path, nrow(dt), uniqueN(dt$slice), paste(sort(unique(dt$k)), collapse = ";"),
    length(missing_cols) == 0,
    "vascular_niche/myeloid_niche/microglia_niche/macrophage_monocyte_niche/neuron_control_niche plus optional malignant weights",
    "spot neighborhood for clustering; slice-level recurrence for interpretation",
    "No CN labels are generated in preflight"
  )
)
fwrite(input_audit, file.path(out_dir, "A4_CN_input_audit.csv"))

slice_audit <- dt[, .(
  n_spots = .N,
  realized_neighbors_median = median(realized_neighbors, na.rm = TRUE),
  realized_neighbors_min = min(realized_neighbors, na.rm = TRUE),
  realized_neighbors_max = max(realized_neighbors, na.rm = TRUE),
  vascular_niche_mean = mean(vascular_niche, na.rm = TRUE),
  myeloid_niche_mean = mean(myeloid_niche, na.rm = TRUE),
  neuron_control_niche_mean = mean(neuron_control_niche, na.rm = TRUE),
  mes_lineage_mean = mean(`MES-lineage`, na.rm = TRUE)
), by = slice][order(-n_spots)]
fwrite(slice_audit, file.path(out_dir, "A4_CN_slice_neighborhood_audit.csv"))

feature_plan <- data.table(
  feature = c("vascular_niche", "myeloid_niche", "microglia_niche", "macrophage_monocyte_niche", "neuron_control_niche", "MES-lineage", "MES-V", "MES-I"),
  planned_role = c(
    "primary context feature",
    "context feature; must respect locked myeloid spatial null",
    "myeloid subtype context",
    "myeloid subtype context",
    "negative-control context",
    "malignant response/context summary",
    "mechanism-subject subtype context",
    "MES-lineage companion context"
  ),
  writing_boundary = c(
    "vascular-rich CN may support local vascular proximity only if reproducible",
    "cannot revive TAM recruitment",
    "cannot revive TAM recruitment",
    "cannot revive TAM recruitment",
    "if positive with claim, downgrade to internal-audit",
    "spot-level RCTD weight, not cell fraction",
    "spot-level RCTD weight, not single-cell location",
    "spot-level RCTD weight, not single-cell location"
  )
)
fwrite(feature_plan, file.path(out_dir, "A4_CN_feature_plan.csv"))

planned_outputs <- data.table(
  output = c(
    "A4_CN_labels_per_spot.csv",
    "A4_CN_composition_summary.csv",
    "A4_CN_recurrence_by_slice.csv",
    "A4_CN_MESV_vascular_enrichment_audit.csv",
    "A4_CN_parameters.csv"
  ),
  role = c(
    "future cellular neighborhood assignment",
    "future neighborhood composition table",
    "future cross-slice recurrence audit",
    "future strict-gate candidate mainline support audit",
    "future clustering/package/seed/feature provenance"
  )
)
fwrite(planned_outputs, file.path(out_dir, "A4_CN_planned_outputs.csv"))

cat("\n[STOP A4 cellular neighborhood preflight] Inputs audited. No CN clustering was run.\n")
