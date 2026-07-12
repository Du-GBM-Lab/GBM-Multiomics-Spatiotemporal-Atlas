#!/usr/bin/env Rscript
# =============================================================================
# R9 batch 1 | A3 squidpy neighborhood enrichment / co-occurrence preflight
# Purpose:
#   Prepare the input audit and future AnnData conversion plan for squidpy.
#   No Python/squidpy spatial graph, enrichment, or co-occurrence is computed.
#
# STOP:
#   Do not run squidpy.gr.spatial_neighbors, nhood_enrichment, or co_occurrence
#   in this script.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
})

base_dir <- getwd()
weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
score_path <- file.path(base_dir, "tables/R9_stage3_signature_scores/R9_stage3_IVY_hypoxia_scores_per_spot.csv")
out_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A3_squidpy")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

python_available <- nzchar(Sys.which("python"))
py_import <- function(module) {
  if (!python_available) return(FALSE)
  cmd <- sprintf("import importlib.util; raise SystemExit(0 if importlib.util.find_spec('%s') else 1)", module)
  status <- suppressWarnings(system2("python", c("-c", shQuote(cmd)), stdout = FALSE, stderr = FALSE))
  identical(status, 0L)
}

method_docs <- data.table(
  analysis = "A3_squidpy",
  source = c("squidpy spatial_neighbors API", "squidpy neighborhood enrichment example", "single-cell best practices neighborhood chapter"),
  url = c(
    "https://squidpy.readthedocs.io/en/stable/api/squidpy.gr.spatial_neighbors.html",
    "https://squidpy.readthedocs.io/en/stable/notebooks/examples/graph/compute_nhood_enrichment.html",
    "https://www.sc-best-practices.org/spatial/neighborhood.html"
  ),
  local_decision = c(
    "spatial_neighbors will be built per slice/library using spatial coordinates; n_neighs=6 is the planned Visium-scale primary graph.",
    "nhood_enrichment compares observed graph proximity against permutations and returns z-scores; this is landscape support, not the locked v3-B null.",
    "co-occurrence/enrichment outputs require neuron control and curve-shape review before any R9 writing."
  )
)
fwrite(method_docs, file.path(out_dir, "A3_squidpy_method_docs.csv"))

pkg_audit <- data.table(
  runtime_or_package = c("python", "anndata", "scanpy", "squidpy", "pandas", "numpy"),
  available = c(
    python_available,
    py_import("anndata"),
    py_import("scanpy"),
    py_import("squidpy"),
    py_import("pandas"),
    py_import("numpy")
  )
)
fwrite(pkg_audit, file.path(out_dir, "A3_squidpy_package_audit.csv"))

weights <- as.data.table(qs2::qs_read(weights_path))
required_cols <- c("spot_id", "slice", "image", "x", "y", "Subtype1", "Subtype2", "Subtype3", "Subtype4",
                   "Endothelial", "Mural cells", "Neurons")
missing_cols <- setdiff(required_cols, names(weights))
if (length(missing_cols)) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

class_plan <- data.table(
  planned_obs_column = c("malignant_subtype_high", "context_high", "neuron_control_high"),
  planned_rule = c(
    "per-slice top-decile RCTD weight for NPC-P/OPC-M/MES-V/MES-I or mutually exclusive max-label sensitivity",
    "per-slice top-decile vascular/myeloid/MVP-like/hypoxia-like contextual score, primary focus vascular",
    "per-slice top-decile neuron RCTD weight as negative-control class"
  ),
  boundary = c(
    "spot-level deconvolution class, not single-cell label",
    "landscape screen only; formal proximity remains v3-B where needed",
    "if neuron follows the claimed pattern, result is internal-audit"
  )
)
fwrite(class_plan, file.path(out_dir, "A3_squidpy_class_definition_plan.csv"))

input_audit <- data.table(
  item = c("weights_path", "score_path", "score_path_exists", "n_spots", "n_slices",
           "planned_spatial_key", "planned_library_key", "planned_primary_graph",
           "planned_outputs_boundary", "stop_boundary"),
  value = c(weights_path, score_path, file.exists(score_path), nrow(weights), uniqueN(weights$slice),
            "obsm['spatial'] from x/y", "obs['slice']", "squidpy.gr.spatial_neighbors(..., library_key='slice', coord_type='generic', n_neighs=6)",
            "z-scores and co-occurrence curves are supplementary landscape evidence unless confirmed by neuron/null gates",
            "No AnnData object or squidpy result is generated in preflight")
)
fwrite(input_audit, file.path(out_dir, "A3_squidpy_input_audit.csv"))

planned_outputs <- data.table(
  output = c(
    "A3_squidpy_input_obs.csv",
    "A3_squidpy_input_X_placeholder_or_scores.h5ad",
    "A3_nhood_enrichment_zscore.csv",
    "A3_co_occurrence_curves.csv",
    "A3_neuron_control_audit.csv",
    "A3_squidpy_parameters.csv"
  ),
  role = c(
    "future AnnData obs table",
    "future AnnData container for squidpy",
    "future neighborhood enrichment matrix",
    "future distance-dependent co-occurrence curves",
    "future negative-control decision table",
    "future Python/package/version/procedure provenance"
  )
)
fwrite(planned_outputs, file.path(out_dir, "A3_squidpy_planned_outputs.csv"))

cat("\n[STOP A3 squidpy preflight] Inputs and Python package state audited. No squidpy graph or enrichment was run.\n")
