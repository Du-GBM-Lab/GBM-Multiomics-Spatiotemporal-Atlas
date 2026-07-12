#!/usr/bin/env Rscript
# =============================================================================
# R9-5 | FOSL1 regulon score preflight
# Purpose:
#   Build a spatial score for the R6-source FOSL1(+) SCENIC target set before
#   any co-enrichment test.
#
# Scope:
#   - Source set: R6 scenic_grn_edgelist.csv, TF == FOSL1
#   - Remove PLAUR / PLAU to avoid circularity with R9-5 PLAUR/PLAU analyses
#   - Score per Visium spot using the same SCT/AddModuleScore口径 as Stage 3
#   - Report target count, SCT/Spatial coverage, and score/target-detection
#     density. No Ripley/random-labeling co-enrichment is run here.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(Matrix)
  library(qs2)
})

set.seed(1)

base_dir <- getwd()
st_path <- "<DATA_ROOT>/项目/分型/分型代码/0.对象/5.ST_merge.rds"
weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
edge_path <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/08_发育时间_TF_ATAC验证/02_TF_regulon_SCENIC/outputs/scenic_grn_edgelist.csv"
out_dir <- file.path(base_dir, "tables/R9_5_spatial_gene_coenrichment")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

exclude_genes <- c("PLAUR", "PLAU")
score_prefix <- "R9FOSL1reg"

## ---- Build R6-source target set -------------------------------------------
edges <- fread(edge_path)
stopifnot(all(c("TF", "target", "importance") %in% names(edges)))
fosl1_edges <- edges[TF == "FOSL1"]
if (!nrow(fosl1_edges)) stop("No FOSL1 targets found in scenic_grn_edgelist.csv")

targets_raw <- unique(as.character(fosl1_edges$target))
targets_no_plaur_plau <- setdiff(targets_raw, exclude_genes)
removed <- intersect(targets_raw, exclude_genes)

target_dt <- data.table(gene = targets_no_plaur_plau)

## ---- Load spatial object and coverage -------------------------------------
st <- readRDS(st_path)
weights <- as.data.table(qs2::qs_read(weights_path))
stopifnot(all(c("spot_id", "slice", "image", "x", "y", "Subtype3", "Subtype4",
                "Endothelial", "Mural cells", "Microglial", "Macrophages",
                "Monocytes", "Neurons") %in% names(weights)))

sct_genes <- rownames(st[["SCT"]])
spatial_genes <- rownames(st[["Spatial"]])
targets_sct <- intersect(targets_no_plaur_plau, sct_genes)
targets_spatial <- intersect(targets_no_plaur_plau, spatial_genes)
if (length(targets_sct) < 10) stop("Too few FOSL1 target genes in SCT assay.")

target_dt[, in_SCT := gene %in% sct_genes]
target_dt[, in_Spatial := gene %in% spatial_genes]
target_dt <- merge(
  target_dt,
  fosl1_edges[, .(gene = target, importance)][, .(importance = max(importance, na.rm = TRUE)), by = gene],
  by = "gene",
  all.x = TRUE
)
setorder(target_dt, -importance, gene)

coverage <- data.table(
  metric = c("raw_targets_from_R6_grn", "removed_PLAUR_PLAU", "targets_after_exclusion",
             "targets_in_SCT", "targets_in_Spatial", "SCT_coverage_fraction",
             "Spatial_coverage_fraction"),
  value = c(length(targets_raw), paste(removed, collapse = ";"), length(targets_no_plaur_plau),
            length(targets_sct), length(targets_spatial),
            length(targets_sct) / length(targets_no_plaur_plau),
            length(targets_spatial) / length(targets_no_plaur_plau))
)

## ---- Score with SCT AddModuleScore ----------------------------------------
DefaultAssay(st) <- "SCT"
st <- AddModuleScore(
  object = st,
  features = list(targets_sct),
  assay = "SCT",
  name = score_prefix,
  seed = 1,
  nbin = 24,
  ctrl = 100,
  slot = "data"
)
score_col <- paste0(score_prefix, "1")
stopifnot(score_col %in% colnames(st@meta.data))

## ---- Raw target detection density -----------------------------------------
sp_layers <- grep("^counts", Layers(st[["Spatial"]]), value = TRUE)
if (!length(sp_layers)) stop("No Spatial counts layers found.")
sp_counts <- do.call(cbind, lapply(sp_layers, function(ly) {
  LayerData(st[["Spatial"]], layer = ly)[targets_spatial, , drop = FALSE]
}))
missing_spots <- setdiff(weights$spot_id, colnames(sp_counts))
if (length(missing_spots)) {
  stop("Spatial counts missing spots: ", paste(head(missing_spots, 10), collapse = ", "))
}
sp_counts <- sp_counts[targets_spatial, weights$spot_id, drop = FALSE]

det_n <- Matrix::colSums(sp_counts > 0)
det_frac <- det_n / length(targets_spatial)

meta <- as.data.table(st@meta.data, keep.rownames = "spot_id")
score_dt <- data.table(
  spot_id = weights$spot_id,
  slice = weights$slice,
  image = weights$image,
  x = weights$x,
  y = weights$y,
  `MES-V` = weights$Subtype3,
  `MES-I` = weights$Subtype4,
  `MES-lineage` = weights$Subtype3 + weights$Subtype4,
  vascular = weights$Endothelial + weights[["Mural cells"]],
  myeloid = weights$Macrophages + weights$Microglial + weights$Monocytes,
  neuron_control = weights$Neurons,
  FOSL1_regulon_score = meta[[score_col]][match(weights$spot_id, meta$spot_id)],
  FOSL1_regulon_detected_targets_n = as.numeric(det_n),
  FOSL1_regulon_detected_targets_frac = as.numeric(det_frac)
)

summarize_density <- function(dt) {
  score <- dt$FOSL1_regulon_score
  dens <- dt$FOSL1_regulon_detected_targets_frac
  top10_cut <- as.numeric(quantile(score, 0.90, na.rm = TRUE))
  top10 <- score >= top10_cut
  data.table(
    n_spots = nrow(dt),
    score_min = min(score, na.rm = TRUE),
    score_q25 = as.numeric(quantile(score, 0.25, na.rm = TRUE)),
    score_median = median(score, na.rm = TRUE),
    score_q75 = as.numeric(quantile(score, 0.75, na.rm = TRUE)),
    score_q90 = as.numeric(quantile(score, 0.90, na.rm = TRUE)),
    score_q95 = as.numeric(quantile(score, 0.95, na.rm = TRUE)),
    score_q99 = as.numeric(quantile(score, 0.99, na.rm = TRUE)),
    score_top10_cutoff = top10_cut,
    top10_n = sum(top10),
    target_detect_frac_median = median(dens, na.rm = TRUE),
    target_detect_frac_q90 = as.numeric(quantile(dens, 0.90, na.rm = TRUE)),
    target_detect_frac_q95 = as.numeric(quantile(dens, 0.95, na.rm = TRUE)),
    top10_target_detect_frac_median = median(dens[top10], na.rm = TRUE),
    top10_target_detect_frac_q25 = as.numeric(quantile(dens[top10], 0.25, na.rm = TRUE)),
    top10_target_detect_frac_q75 = as.numeric(quantile(dens[top10], 0.75, na.rm = TRUE)),
    score_top10_cutoff_ties = sum(score == top10_cut, na.rm = TRUE)
  )
}

overall_density <- summarize_density(score_dt)
overall_density[, slice := "ALL"]
setcolorder(overall_density, c("slice", setdiff(names(overall_density), "slice")))

per_slice_density <- score_dt[, summarize_density(.SD), by = slice]

flags <- per_slice_density[, .(
  n_slices = .N,
  median_score_top10_cutoff = median(score_top10_cutoff, na.rm = TRUE),
  min_score_top10_cutoff = min(score_top10_cutoff, na.rm = TRUE),
  max_score_top10_cutoff = max(score_top10_cutoff, na.rm = TRUE),
  n_slices_top10_cutoff_nonpositive = sum(score_top10_cutoff <= 0, na.rm = TRUE),
  median_target_detect_frac = median(target_detect_frac_median, na.rm = TRUE),
  median_top10_target_detect_frac = median(top10_target_detect_frac_median, na.rm = TRUE),
  min_top10_target_detect_frac = min(top10_target_detect_frac_median, na.rm = TRUE),
  max_top10_target_detect_frac = max(top10_target_detect_frac_median, na.rm = TRUE)
)]

params <- data.table(
  parameter = c("target_source", "target_rule", "excluded_genes", "scoring_method",
                "assay", "slot", "nbin", "ctrl", "seed", "stop_rule"),
  value = c(edge_path, "TF == FOSL1 from R6 SCENIC GRN edgelist",
            paste(exclude_genes, collapse = ";"), "Seurat::AddModuleScore",
            "SCT", "data", "24", "100", "1",
            "No Ripley/random-labeling co-enrichment run in this script")
)

fwrite(target_dt, file.path(out_dir, "R9_5_FOSL1_regulon_targets_R6_no_PLAUR_PLAU.csv"))
fwrite(coverage, file.path(out_dir, "R9_5_FOSL1_regulon_target_coverage.csv"))
fwrite(score_dt, file.path(out_dir, "R9_5_FOSL1_regulon_spot_scores.csv"))
fwrite(overall_density, file.path(out_dir, "R9_5_FOSL1_regulon_score_density_overall.csv"))
fwrite(per_slice_density, file.path(out_dir, "R9_5_FOSL1_regulon_score_density_per_slice.csv"))
fwrite(flags, file.path(out_dir, "R9_5_FOSL1_regulon_score_density_flags.csv"))
fwrite(params, file.path(out_dir, "R9_5_FOSL1_regulon_score_parameters.csv"))

cat("\n== FOSL1 regulon target coverage:\n")
print(coverage)

cat("\n== Overall FOSL1 regulon score density:\n")
print(overall_density[, .(
  n_spots, score_median = signif(score_median, 3),
  score_q90 = signif(score_q90, 3),
  score_q95 = signif(score_q95, 3),
  score_top10_cutoff = signif(score_top10_cutoff, 3),
  top10_n,
  target_detect_frac_median = round(target_detect_frac_median, 4),
  top10_target_detect_frac_median = round(top10_target_detect_frac_median, 4)
)])

cat("\n== Per-slice density flags:\n")
print(flags[, .(
  n_slices,
  median_score_top10_cutoff = signif(median_score_top10_cutoff, 3),
  n_slices_top10_cutoff_nonpositive,
  median_target_detect_frac = round(median_target_detect_frac, 4),
  median_top10_target_detect_frac = round(median_top10_target_detect_frac, 4),
  min_top10_target_detect_frac = round(min_top10_target_detect_frac, 4)
)])

cat("\n[STOP R9-5 FOSL1 regulon score] Report target counts and density before any proximity test.\n")
