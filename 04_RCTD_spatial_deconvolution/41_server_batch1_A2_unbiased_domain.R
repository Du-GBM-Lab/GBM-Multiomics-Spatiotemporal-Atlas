#!/usr/bin/env Rscript
# =============================================================================
# R9 batch 1 | A2 unbiased spatial domain / niche clustering on server
# Purpose:
#   Build a data-driven spatial domain layer from existing RCTD-derived
#   neighborhood composition scores. This is independent from IVY ROI-like labels.
#
# Boundaries:
#   - Domains are spot-level, RCTD/neighborhood-derived, not histology regions.
#   - This is a landscape layer; it does not replace v3-B proximity statistics.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(cluster)
})

set.seed(1)

base_dir <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
score_path <- file.path(base_dir, "tables/R9_C3_1_neighborhood_niche_scores.csv")
out_dir <- file.path(base_dir, "outputs/R9_batch1_unbiased_landscape/A2_unbiased_domain")
tab_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A2_unbiased_domain")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

k_domain <- as.integer(strsplit(Sys.getenv("R9_A2_DOMAIN_K", "4,5,6,7,8"), ",")[[1]])
primary_k <- as.integer(Sys.getenv("R9_A2_PRIMARY_K", "6"))

dt <- fread(score_path)
required <- c("spot_id", "slice", "image", "x", "y", "k", "realized_neighbors",
              "MES-lineage", "MES-V", "MES-I", "vascular_niche", "myeloid_niche",
              "microglia_niche", "macrophage_monocyte_niche", "neuron_control_niche")
miss <- setdiff(required, names(dt))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
dt <- dt[k == 6]

features <- c("MES-lineage", "MES-V", "MES-I", "vascular_niche", "myeloid_niche",
              "microglia_niche", "macrophage_monocyte_niche", "neuron_control_niche")

scale_by_slice <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

feat_dt <- copy(dt)
for (ff in features) {
  feat_dt[, paste0(ff, "_z") := scale_by_slice(get(ff)), by = slice]
}
feat_cols <- paste0(features, "_z")
X <- as.matrix(feat_dt[, ..feat_cols])
X[!is.finite(X)] <- 0

sample_silhouette <- function(X, cl, max_n = 8000L) {
  n <- nrow(X)
  idx <- if (n > max_n) sample.int(n, max_n) else seq_len(n)
  mean(cluster::silhouette(cl[idx], dist(X[idx, , drop = FALSE]))[, "sil_width"])
}

fit_rows <- lapply(k_domain, function(kk) {
  km <- kmeans(X, centers = kk, nstart = 50, iter.max = 100)
  data.table(
    k_domain = kk,
    tot_withinss = km$tot.withinss,
    betweenss = km$betweenss,
    silhouette_sample_mean = sample_silhouette(X, km$cluster)
  )
})
fit_audit <- rbindlist(fit_rows)
fwrite(fit_audit, file.path(tab_dir, "A2_domain_k_selection_audit.csv"))

if (!primary_k %in% k_domain) primary_k <- fit_audit[which.max(silhouette_sample_mean), k_domain]
km <- kmeans(X, centers = primary_k, nstart = 100, iter.max = 100)

labels <- feat_dt[, .(spot_id, slice, image, x, y)]
labels[, domain_k := primary_k]
labels[, domain_id := paste0("D", sprintf("%02d", km$cluster))]
fwrite(labels, file.path(tab_dir, "A2_domain_labels_per_spot.csv"))

annot_dt <- cbind(feat_dt[, .(spot_id, slice, image, x, y, realized_neighbors)], domain_id = labels$domain_id, dt[, ..features])
domain_summary <- annot_dt[, lapply(.SD, mean, na.rm = TRUE), by = domain_id, .SDcols = features][order(domain_id)]
domain_counts <- annot_dt[, .(n_spots = .N, n_slices = uniqueN(slice)), by = domain_id][order(domain_id)]
domain_summary <- merge(domain_counts, domain_summary, by = "domain_id")
domain_summary[, domain_role_suggestion := fifelse(
  vascular_niche == max(vascular_niche), "vascular-rich candidate",
  fifelse(neuron_control_niche == max(neuron_control_niche), "neuron-control-rich candidate", "landscape candidate")
)]
fwrite(domain_summary, file.path(tab_dir, "A2_domain_composition_summary.csv"))

slice_domain <- annot_dt[, .(n_spots = .N), by = .(slice, domain_id)]
slice_total <- annot_dt[, .(slice_spots = .N), by = slice]
slice_domain <- merge(slice_domain, slice_total, by = "slice")
slice_domain[, fraction := n_spots / slice_spots]
fwrite(slice_domain, file.path(tab_dir, "A2_domain_recurrence_by_slice.csv"))

params <- data.table(
  parameter = c("input", "feature_set", "scaling", "clustering", "candidate_k", "primary_k", "boundary"),
  value = c(score_path, paste(features, collapse = ";"), "within-slice z-score",
            "kmeans on RCTD-derived neighborhood composition", paste(k_domain, collapse = ";"), primary_k,
            "unbiased landscape layer; not IVY Location and not a proximity proof")
)
fwrite(params, file.path(tab_dir, "A2_domain_parameters.csv"))

cat("\n[STOP A2 domain] Unbiased domain labels and composition summaries written. Review before figure/claim decisions.\n")
