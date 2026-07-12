#!/usr/bin/env Rscript
# =============================================================================
# R9 batch 1 | A2 unsupervised spatial domain / niche preflight
# Purpose:
#   Audit inputs and method choices for data-driven spatial domain/niche
#   discovery. No clustering or BayesSpace MCMC is run here.
#
# STOP:
#   Do not run BayesSpace, graph clustering, k-means, Leiden, or any domain
#   assignment in this script.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
})

base_dir <- getwd()
weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
score_path <- file.path(base_dir, "tables/R9_stage3_signature_scores/R9_stage3_IVY_hypoxia_scores_per_spot.csv")
out_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A2_unsupervised_domain")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

method_docs <- data.table(
  analysis = "A2_unsupervised_domain",
  source = c("BayesSpace Bioconductor manual", "BayesSpace GitHub"),
  url = c(
    "https://bioconductor.org/packages/release/bioc/manuals/BayesSpace/man/BayesSpace.pdf",
    "https://github.com/edward130603/BayesSpace"
  ),
  local_decision = c(
    "BayesSpace spatialPreprocess/spatialCluster/qTune are legitimate spatial-domain options after input conversion to SingleCellExperiment.",
    "For R9 batch 1, RCTD-composition graph clustering is a simpler parallel track and avoids treating IVY ROI-like labels as ground truth."
  )
)
fwrite(method_docs, file.path(out_dir, "A2_domain_method_docs.csv"))

pkg_audit <- data.table(
  package = c("BayesSpace", "SingleCellExperiment", "mclust", "igraph", "data.table", "qs2"),
  available = c(
    requireNamespace("BayesSpace", quietly = TRUE),
    requireNamespace("SingleCellExperiment", quietly = TRUE),
    requireNamespace("mclust", quietly = TRUE),
    requireNamespace("igraph", quietly = TRUE),
    requireNamespace("data.table", quietly = TRUE),
    requireNamespace("qs2", quietly = TRUE)
  )
)
fwrite(pkg_audit, file.path(out_dir, "A2_domain_package_audit.csv"))

weights <- as.data.table(qs2::qs_read(weights_path))
required_weight_cols <- c(
  "spot_id", "slice", "image", "x", "y",
  "Subtype1", "Subtype2", "Subtype3", "Subtype4",
  "Endothelial", "Mural cells", "Macrophages", "Microglial", "Monocytes",
  "Neurons"
)
missing_weight_cols <- setdiff(required_weight_cols, names(weights))
if (length(missing_weight_cols)) {
  stop("Missing required RCTD weight columns: ", paste(missing_weight_cols, collapse = ", "))
}

score_exists <- file.exists(score_path)
score_cols <- character(0)
if (score_exists) {
  score_header <- fread(score_path, nrows = 0)
  score_cols <- names(score_header)
}

input_audit <- data.table(
  item = c(
    "weights_path", "score_path", "score_path_exists", "n_spots",
    "n_slices", "n_weight_columns", "required_weight_columns_present",
    "planned_primary_features", "planned_context_features", "stop_boundary"
  ),
  value = c(
    weights_path, score_path, score_exists, nrow(weights),
    uniqueN(weights$slice), length(names(weights)), length(missing_weight_cols) == 0,
    "RCTD composition: four malignant subtypes plus vascular/myeloid/neuron context",
    if (score_exists) paste(intersect(score_cols, c("IVY_CT", "IVY_CTmvp", "IVY_CTpan", "IVY_LE", "Hypoxia_Buffa")), collapse = ";") else "not_available",
    "No unsupervised domain labels are generated in preflight"
  )
)
fwrite(input_audit, file.path(out_dir, "A2_domain_input_audit.csv"))

slice_audit <- weights[, .(
  n_spots = .N,
  x_min = min(x), x_max = max(x),
  y_min = min(y), y_max = max(y),
  mes_lineage_mean = mean(Subtype3 + Subtype4),
  vascular_mean = mean(Endothelial + `Mural cells`),
  neuron_mean = mean(Neurons)
), by = slice][order(-n_spots)]
fwrite(slice_audit, file.path(out_dir, "A2_domain_slice_audit.csv"))

planned_outputs <- data.table(
  output = c(
    "A2_domain_labels_per_spot.csv",
    "A2_domain_composition_by_slice.csv",
    "A2_domain_marker_or_weight_summary.csv",
    "A2_domain_reproducibility_audit.csv",
    "A2_domain_parameters.csv"
  ),
  role = c(
    "future unsupervised domain assignment",
    "future per-slice domain abundance and RCTD composition",
    "future domain interpretation table",
    "future cross-slice reproducibility/fragmentation audit",
    "future method/package/seed/feature-set provenance"
  )
)
fwrite(planned_outputs, file.path(out_dir, "A2_domain_planned_outputs.csv"))

cat("\n[STOP A2 domain preflight] Inputs and package state audited. No domain clustering was run.\n")
