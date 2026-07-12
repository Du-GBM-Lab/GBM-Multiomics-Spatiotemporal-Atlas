#!/usr/bin/env Rscript
# =============================================================================
# R9 batch 1 | A4 cellular neighborhood on server
# Purpose:
#   Cluster recurrent local composition patterns from RCTD-derived neighborhood
#   scores. This is a Visium spot-neighborhood analysis, not single-cell contact.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cluster)
})

set.seed(1)

base_dir <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
niche_score_path <- file.path(base_dir, "tables/R9_C3_1_neighborhood_niche_scores.csv")
out_dir <- file.path(base_dir, "outputs/R9_batch1_unbiased_landscape/A4_cellular_neighborhood")
tab_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A4_cellular_neighborhood")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

k_values <- as.integer(strsplit(Sys.getenv("R9_A4_CN_K", "5,6,7,8,9,10"), ",")[[1]])
primary_k <- as.integer(Sys.getenv("R9_A4_PRIMARY_K", "8"))
max_silhouette_n <- as.integer(Sys.getenv("R9_A4_SILHOUETTE_N", "8000"))

dt <- fread(niche_score_path)
required <- c("spot_id", "slice", "image", "x", "y", "k", "realized_neighbors",
              "MES-lineage", "MES-V", "MES-I",
              "vascular_niche", "myeloid_niche", "microglia_niche",
              "macrophage_monocyte_niche", "neuron_control_niche")
miss <- setdiff(required, names(dt))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
dt <- dt[k == 6]

features <- c("vascular_niche", "myeloid_niche", "microglia_niche",
              "macrophage_monocyte_niche", "neuron_control_niche",
              "MES-lineage", "MES-V", "MES-I")

scale_vec <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

work <- copy(dt)
for (ff in features) {
  work[, paste0(ff, "_z") := scale_vec(get(ff)), by = slice]
}
X <- as.matrix(work[, paste0(features, "_z"), with = FALSE])
X[!is.finite(X)] <- 0

sil_mean <- function(cl) {
  n <- length(cl)
  idx <- if (n > max_silhouette_n) sample.int(n, max_silhouette_n) else seq_len(n)
  mean(cluster::silhouette(cl[idx], dist(X[idx, , drop = FALSE]))[, "sil_width"])
}

k_audit <- rbindlist(lapply(k_values, function(kk) {
  km <- kmeans(X, centers = kk, nstart = 50, iter.max = 100)
  data.table(
    k_cn = kk,
    tot_withinss = km$tot.withinss,
    betweenss = km$betweenss,
    silhouette_sample_mean = sil_mean(km$cluster)
  )
}))
fwrite(k_audit, file.path(tab_dir, "A4_CN_k_selection_audit.csv"))

if (!primary_k %in% k_values) primary_k <- k_audit[which.max(silhouette_sample_mean), k_cn]
km <- kmeans(X, centers = primary_k, nstart = 100, iter.max = 100)

labels <- work[, .(spot_id, slice, image, x, y, realized_neighbors)]
labels[, CN_id := paste0("CN", sprintf("%02d", km$cluster))]
labels[, CN_k := primary_k]
fwrite(labels, file.path(tab_dir, "A4_CN_labels_per_spot.csv"))

ann <- cbind(work[, .(spot_id, slice, image, x, y)], CN_id = labels$CN_id, dt[, ..features])
cn_summary <- ann[, lapply(.SD, mean, na.rm = TRUE), by = CN_id, .SDcols = features][order(CN_id)]
cn_counts <- ann[, .(n_spots = .N, n_slices = uniqueN(slice)), by = CN_id][order(CN_id)]
cn_summary <- merge(cn_counts, cn_summary, by = "CN_id")
cn_summary[, vascular_rank := frank(-vascular_niche, ties.method = "dense")]
cn_summary[, mesv_rank := frank(-`MES-V`, ties.method = "dense")]
cn_summary[, neuron_rank := frank(-neuron_control_niche, ties.method = "dense")]
cn_summary[, role_suggestion := fifelse(
  vascular_rank == 1 & mesv_rank <= 3 & neuron_rank > 1,
  "vascular-MES candidate",
  fifelse(neuron_rank == 1, "neuron-control-rich candidate", "landscape candidate")
)]
fwrite(cn_summary, file.path(tab_dir, "A4_CN_composition_summary.csv"))

slice_cn <- ann[, .(n_spots = .N), by = .(slice, CN_id)]
slice_total <- ann[, .(slice_spots = .N), by = slice]
slice_cn <- merge(slice_cn, slice_total, by = "slice")
slice_cn[, fraction := n_spots / slice_spots]
fwrite(slice_cn, file.path(tab_dir, "A4_CN_recurrence_by_slice.csv"))

candidate_audit <- merge(slice_cn, cn_summary[, .(CN_id, role_suggestion)], by = "CN_id")
candidate_audit <- candidate_audit[, .(
  n_slices_present = sum(fraction > 0),
  median_fraction = median(fraction),
  max_fraction = max(fraction)
), by = .(CN_id, role_suggestion)][order(role_suggestion, -median_fraction)]
fwrite(candidate_audit, file.path(tab_dir, "A4_CN_recurrence_summary.csv"))

params <- data.table(
  parameter = c("input", "features", "scaling", "clustering", "candidate_k", "primary_k", "boundary"),
  value = c(niche_score_path, paste(features, collapse = ";"), "within-slice z-score",
            "kmeans on k=6 RCTD-derived neighborhood composition", paste(k_values, collapse = ";"), primary_k,
            "cellular neighborhood is spot-neighborhood composition; no contact/recruitment/interaction claim")
)
fwrite(params, file.path(tab_dir, "A4_CN_parameters.csv"))

cat("\n[STOP A4 CN] Cellular neighborhood labels and summaries written. Review recurrence and neuron control before figure decisions.\n")
