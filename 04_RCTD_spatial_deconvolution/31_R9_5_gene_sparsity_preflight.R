#!/usr/bin/env Rscript
# =============================================================================
# R9-5 | Step 1: FOSL1 / PLAUR / PLAU spatial sparsity preflight
# Purpose:
#   Before any co-enrichment test, audit single-gene detectability in Visium.
#   Detection is calculated from raw Spatial counts. Continuous per-spot scores
#   for later high-spot ranking use the same SCT data slot used in Stage 3.
#
# STOP:
#   No Ripley/random-labeling co-enrichment is run in this script.
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
out_dir <- file.path(base_dir, "tables/R9_5_spatial_gene_coenrichment")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

genes <- c("FOSL1", "PLAUR", "PLAU")

st <- readRDS(st_path)
weights <- as.data.table(qs2::qs_read(weights_path))
stopifnot(all(c("spot_id", "slice", "image", "x", "y", "Subtype3", "Subtype4",
                "Endothelial", "Mural cells", "Neurons") %in% names(weights)))

present <- data.table(
  gene = genes,
  in_Spatial = genes %in% rownames(st[["Spatial"]]),
  in_SCT = genes %in% rownames(st[["SCT"]])
)
if (any(!present$in_Spatial | !present$in_SCT)) {
  fwrite(present, file.path(out_dir, "R9_5_gene_presence.csv"))
  stop("One or more genes missing from Spatial or SCT assay; inspect R9_5_gene_presence.csv")
}

sp_layers <- grep("^counts", Layers(st[["Spatial"]]), value = TRUE)
if (!length(sp_layers)) stop("No Spatial counts layers found.")
sp_counts <- do.call(cbind, lapply(sp_layers, function(ly) {
  LayerData(st[["Spatial"]], layer = ly)[genes, , drop = FALSE]
}))
missing_spots <- setdiff(weights$spot_id, colnames(sp_counts))
if (length(missing_spots)) {
  stop("Spatial counts missing spots: ", paste(head(missing_spots, 10), collapse = ", "))
}
sp_counts <- sp_counts[genes, weights$spot_id, drop = FALSE]

sct_data <- LayerData(st[["SCT"]], layer = "data")[genes, weights$spot_id, drop = FALSE]

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
  neuron_control = weights$Neurons
)

for (g in genes) {
  score_dt[, paste0(g, "_raw_count") := as.numeric(sp_counts[g, ])]
  score_dt[, paste0(g, "_detected") := get(paste0(g, "_raw_count")) > 0]
  score_dt[, paste0(g, "_SCT_data") := as.numeric(sct_data[g, ])]
}

gene_summary_one <- function(dt, g) {
  raw <- dt[[paste0(g, "_raw_count")]]
  det <- raw > 0
  sct <- dt[[paste0(g, "_SCT_data")]]
  top10_raw_cut <- as.numeric(quantile(raw, 0.90, na.rm = TRUE))
  top10_sct_cut <- as.numeric(quantile(sct, 0.90, na.rm = TRUE))
  top10_sct <- sct >= top10_sct_cut
  data.table(
    gene = g,
    n_spots = nrow(dt),
    n_detected = sum(det),
    pct_detected = mean(det),
    mean_raw_count = mean(raw),
    median_raw_count = median(raw),
    q90_raw_count = as.numeric(quantile(raw, 0.90, na.rm = TRUE)),
    q95_raw_count = as.numeric(quantile(raw, 0.95, na.rm = TRUE)),
    q99_raw_count = as.numeric(quantile(raw, 0.99, na.rm = TRUE)),
    mean_SCT_data = mean(sct),
    median_SCT_data = median(sct),
    q90_SCT_data = as.numeric(quantile(sct, 0.90, na.rm = TRUE)),
    q95_SCT_data = as.numeric(quantile(sct, 0.95, na.rm = TRUE)),
    q99_SCT_data = as.numeric(quantile(sct, 0.99, na.rm = TRUE)),
    top10_raw_cutoff = top10_raw_cut,
    top10_SCT_cutoff = top10_sct_cut,
    top10_SCT_n = sum(top10_sct),
    top10_SCT_detected_fraction = if (sum(top10_sct) > 0) mean(det[top10_sct]) else NA_real_,
    raw_top10_cutoff_zero = top10_raw_cut <= 0
  )
}

per_slice <- rbindlist(lapply(sort(unique(score_dt$slice)), function(sl) {
  dt <- score_dt[slice == sl]
  rbindlist(lapply(genes, function(g) gene_summary_one(dt, g)))[, slice := sl][]
}), use.names = TRUE)
setcolorder(per_slice, c("slice", setdiff(names(per_slice), "slice")))

overall <- rbindlist(lapply(genes, function(g) gene_summary_one(score_dt, g)), use.names = TRUE)
overall[, slice := "ALL"]
setcolorder(overall, c("slice", setdiff(names(overall), "slice")))

flags <- per_slice[, .(
  n_slices = .N,
  median_pct_detected = median(pct_detected),
  min_pct_detected = min(pct_detected),
  max_pct_detected = max(pct_detected),
  n_slices_detect_lt_0p05 = sum(pct_detected < 0.05),
  n_slices_detect_lt_0p10 = sum(pct_detected < 0.10),
  n_slices_raw_top10_cutoff_zero = sum(raw_top10_cutoff_zero),
  median_top10_SCT_detected_fraction = median(top10_SCT_detected_fraction, na.rm = TRUE)
), by = gene]

params <- data.table(
  parameter = c("detection_assay", "detection_slot", "continuous_score_assay",
                "continuous_score_slot", "genes", "stop_rule"),
  value = c("Spatial", "counts", "SCT", "data", paste(genes, collapse = ";"),
            "No co-enrichment run; report sparsity before R9-5 proximity.")
)

fwrite(present, file.path(out_dir, "R9_5_gene_presence.csv"))
fwrite(score_dt, file.path(out_dir, "R9_5_gene_spot_scores_SCT_and_raw.csv"))
fwrite(per_slice, file.path(out_dir, "R9_5_gene_sparsity_per_slice.csv"))
fwrite(overall, file.path(out_dir, "R9_5_gene_sparsity_overall.csv"))
fwrite(flags, file.path(out_dir, "R9_5_gene_sparsity_flags.csv"))
fwrite(params, file.path(out_dir, "R9_5_gene_sparsity_parameters.csv"))

cat("\n== R9-5 gene presence:\n")
print(present)

cat("\n== Overall sparsity:\n")
print(overall[, .(
  gene, n_spots, n_detected,
  pct_detected = round(pct_detected, 4),
  mean_raw_count = signif(mean_raw_count, 3),
  q90_raw_count, q95_raw_count, q99_raw_count,
  q90_SCT_data = signif(q90_SCT_data, 3),
  top10_SCT_detected_fraction = round(top10_SCT_detected_fraction, 4),
  raw_top10_cutoff_zero
)])

cat("\n== Per-slice sparsity flags:\n")
print(flags[, .(
  gene, n_slices,
  median_pct_detected = round(median_pct_detected, 4),
  min_pct_detected = round(min_pct_detected, 4),
  max_pct_detected = round(max_pct_detected, 4),
  n_slices_detect_lt_0p05,
  n_slices_detect_lt_0p10,
  n_slices_raw_top10_cutoff_zero,
  median_top10_SCT_detected_fraction = round(median_top10_SCT_detected_fraction, 4)
)])

cat("\n[STOP R9-5 Step 1] Report expression sparsity before any spatial co-enrichment test.\n")
