#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-0: ecological niche analysis preflight
# Purpose: audit inputs and lock implementation parameters before running C3.
# No formal niche enrichment, no spatial null model, no biological conclusion.
# =============================================================================

suppressPackageStartupMessages({
  library(qs2)
  library(data.table)
  library(Seurat)
})

base_dir <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
weights_path <- file.path(base_dir, "outputs/R9_A2_RCTD_weights_allslices_long.qs2")
st_path <- file.path(base_dir, "data/spatial/5.ST_merge.rds")
ivy_path <- file.path(base_dir, "data/spatial/IVY_gap_signatures.csv")
out_dir <- file.path(base_dir, "outputs/R9_C3_niche_preflight")
tab_dir <- file.path(base_dir, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

big <- as.data.table(qs2::qs_read(weights_path))
required <- c("spot_id", "slice", "x", "y",
              "Subtype3", "Subtype4",
              "Endothelial", "Mural cells",
              "Macrophages", "Microglial", "Monocytes",
              "Oligodendrocytes", "Neurons")
miss <- setdiff(required, names(big))
if (length(miss)) stop("missing required columns in A2 weights: ", paste(miss, collapse = ", "))

meta_cols <- intersect(c("spot_id", "slice", "image", "x", "y"), names(big))
weight_cols <- setdiff(names(big), meta_cols)
weight_cols <- weight_cols[vapply(big[, ..weight_cols], is.numeric, logical(1))]

## 1) weights integrity
row_sum <- rowSums(as.matrix(big[, ..weight_cols]))
weights_qc <- data.table(
  n_spots = nrow(big),
  n_slices = uniqueN(big$slice),
  n_weight_cols = length(weight_cols),
  n_na_weights = sum(is.na(as.matrix(big[, ..weight_cols]))),
  row_sum_min = min(row_sum),
  row_sum_q01 = as.numeric(quantile(row_sum, 0.01)),
  row_sum_median = median(row_sum),
  row_sum_q99 = as.numeric(quantile(row_sum, 0.99)),
  row_sum_max = max(row_sum)
)

slice_counts <- big[, .(n_spots = .N), by = slice][order(slice)]

## 2) coordinate pitch audit using existing x/y coordinates
nearest_dist <- function(dt) {
  xy <- as.matrix(dt[, .(x, y)])
  n <- nrow(xy)
  if (n < 3) return(rep(NA_real_, n))
  out <- numeric(n)
  for (i in seq_len(n)) {
    d <- sqrt((xy[, 1] - xy[i, 1])^2 + (xy[, 2] - xy[i, 2])^2)
    d[i] <- Inf
    out[i] <- min(d, na.rm = TRUE)
  }
  out
}

pitch_list <- lapply(split(big, big$slice), function(dt) {
  nd <- nearest_dist(dt)
  data.table(
    slice = dt$slice[1],
    n_spots = nrow(dt),
    nn_dist_median = median(nd, na.rm = TRUE),
    nn_dist_q05 = as.numeric(quantile(nd, 0.05, na.rm = TRUE)),
    nn_dist_q95 = as.numeric(quantile(nd, 0.95, na.rm = TRUE)),
    distance_threshold_1p5x = 1.5 * median(nd, na.rm = TRUE)
  )
})
pitch_summary <- rbindlist(pitch_list)

global_pitch <- median(pitch_summary$nn_dist_median, na.rm = TRUE)
distance_threshold <- 1.5 * global_pitch

## 3) negative-control candidates
big[, `MES-lineage` := Subtype3 + Subtype4]
neg_candidates <- c("Oligodendrocytes", "Neurons", "Astrocytes", "OPCs")
neg_candidates <- intersect(neg_candidates, weight_cols)

candidate_rows <- lapply(neg_candidates, function(ct) {
  x <- big[[ct]]
  per_slice_cor <- big[, .(
    cor_mes = suppressWarnings(cor(`MES-lineage`, get(ct), method = "spearman")),
    mean_weight = mean(get(ct)),
    q90_weight = as.numeric(quantile(get(ct), 0.90)),
    frac_gt_005 = mean(get(ct) > 0.05)
  ), by = slice]
  data.table(
    candidate = ct,
    global_mean = mean(x),
    global_q90 = as.numeric(quantile(x, 0.90)),
    global_frac_gt_005 = mean(x > 0.05),
    median_slice_mean = median(per_slice_cor$mean_weight),
    median_slice_q90 = median(per_slice_cor$q90_weight),
    median_slice_frac_gt_005 = median(per_slice_cor$frac_gt_005),
    median_slice_cor_with_mes_lineage = median(per_slice_cor$cor_mes, na.rm = TRUE),
    n_slice_cor_positive = sum(per_slice_cor$cor_mes > 0, na.rm = TRUE),
    n_slice_cor_negative = sum(per_slice_cor$cor_mes < 0, na.rm = TRUE)
  )
})
negative_control_audit <- rbindlist(candidate_rows)

## 4) IVY signature audit and gene overlap
ivy <- fread(ivy_path)
setnames(ivy, old = names(ivy), new = make.names(names(ivy)))
if (!all(c("NAME", "Assigned") %in% names(ivy))) {
  stop("IVY signature CSV must contain NAME and Assigned columns after make.names.")
}
ivy[, NAME := as.character(NAME)]
ivy[, Assigned := as.character(Assigned)]
ivy_counts <- ivy[, .(n_genes = uniqueN(NAME)), by = Assigned][order(Assigned)]

st <- readRDS(st_path)
spatial_genes <- rownames(st[["Spatial"]])
sct_genes <- rownames(st[["SCT"]])
ivy_overlap <- ivy[, .(
  n_genes = uniqueN(NAME),
  n_in_spatial = uniqueN(NAME[NAME %in% spatial_genes]),
  n_in_sct = uniqueN(NAME[NAME %in% sct_genes]),
  frac_in_spatial = uniqueN(NAME[NAME %in% spatial_genes]) / uniqueN(NAME),
  frac_in_sct = uniqueN(NAME[NAME %in% sct_genes]) / uniqueN(NAME)
), by = Assigned][order(Assigned)]

## 5) patient/slice hint from slice names
slice_patient <- data.table(slice = unique(big$slice))
slice_patient[, patient_like_id := sub("^#([^_]+).*", "\\1", slice)]
patient_counts <- slice_patient[, .(n_slices = .N), by = patient_like_id][order(-n_slices, patient_like_id)]

recommended_negative <- negative_control_audit[
  order(abs(median_slice_cor_with_mes_lineage), -median_slice_frac_gt_005)
][1, candidate]

param_recommendation <- data.table(
  parameter = c("coordinate_source", "main_k", "sensitivity_k",
                "distance_threshold", "negative_control_niche",
                "ivy_pan_name", "ivy_mvp_name",
                "statistical_unit"),
  value = c("A2 x/y coordinates; array_row/array_col not present in A2 table",
            "6", "12",
            sprintf("%.3f (1.5 x median nearest-neighbor distance %.3f)", distance_threshold, global_pitch),
            recommended_negative,
            "CTpan -> reconstructed PAN-like score",
            "CTmvp -> reconstructed MVP-like score",
            ifelse(any(patient_counts$n_slices > 1), "patient-aware required", "slice"))
)

qs2::qs_save(
  list(
    weights_qc = weights_qc,
    slice_counts = slice_counts,
    pitch_summary = pitch_summary,
    distance_threshold = distance_threshold,
    negative_control_audit = negative_control_audit,
    ivy_counts = ivy_counts,
    ivy_overlap = ivy_overlap,
    patient_counts = patient_counts,
    param_recommendation = param_recommendation
  ),
  file.path(out_dir, "R9_C3_0_niche_preflight.qs2")
)

fwrite(weights_qc, file.path(tab_dir, "R9_C3_0_weights_integrity.csv"))
fwrite(slice_counts, file.path(tab_dir, "R9_C3_0_slice_spot_counts.csv"))
fwrite(pitch_summary, file.path(tab_dir, "R9_C3_0_spot_pitch_by_slice.csv"))
fwrite(negative_control_audit, file.path(tab_dir, "R9_C3_0_negative_control_candidates.csv"))
fwrite(ivy_counts, file.path(tab_dir, "R9_C3_0_IVY_signature_counts.csv"))
fwrite(ivy_overlap, file.path(tab_dir, "R9_C3_0_IVY_signature_gene_overlap.csv"))
fwrite(patient_counts, file.path(tab_dir, "R9_C3_0_patient_slice_counts.csv"))
fwrite(param_recommendation, file.path(tab_dir, "R9_C3_0_parameter_recommendations.csv"))

cat("\n== weights QC:\n")
print(weights_qc)
cat("\n== spot pitch summary:\n")
print(pitch_summary[, .(
  n_slices = .N,
  median_nn_dist = median(nn_dist_median),
  q05_nn_dist = quantile(nn_dist_median, 0.05),
  q95_nn_dist = quantile(nn_dist_median, 0.95),
  recommended_threshold = distance_threshold
)])
cat("\n== negative-control candidates:\n")
print(negative_control_audit)
cat("\n== IVY signature overlap:\n")
print(ivy_overlap)
cat("\n== patient/slice counts:\n")
print(patient_counts)
cat("\n== parameter recommendations:\n")
print(param_recommendation)
cat("\n[STOP C3-0] Confirm distance threshold and negative-control niche before C3-1.\n")
