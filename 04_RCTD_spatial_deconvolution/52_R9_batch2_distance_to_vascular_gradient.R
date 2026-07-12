#!/usr/bin/env Rscript
# =============================================================================
# R9 Batch 2 project 2 | distance-to-vascular gradient composition
#
# Nature:
#   Descriptive spatial landscape and overview candidate. The vascular axis is
#   defined independently from MES by dominant RCTD Endothelial/Mural spots.
#   Every spot is retained; distance bins are geometric summaries, not score-
#   based spot filtering. This script does not make relation-level p-value
#   claims; any inferential proximity claim must return to locked v3-B.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
  library(RANN)
  library(ggplot2)
})

base_dir <- getwd()
out_dir <- file.path(base_dir, "tables/R9_batch2_landscape_gradient/vascular_distance_gradient")
fig_dir <- file.path(base_dir, "figures/R9_batch2_landscape_gradient/vascular_distance_gradient")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
subtype_map <- c("Subtype1" = "NPC-P", "Subtype2" = "OPC-M", "Subtype3" = "MES-V", "Subtype4" = "MES-I")
metadata_cols <- c("spot_id", "slice", "image", "x", "y")

dominant_label <- function(dt, cols) {
  mat <- as.matrix(dt[, ..cols])
  cols[max.col(mat, ties.method = "first")]
}

rank_bins <- function(x, n_bins = 5L) {
  if (all(is.na(x))) return(rep(NA_integer_, length(x)))
  r <- frank(x, ties.method = "average", na.last = "keep")
  out <- floor(((r - 0.5) / sum(!is.na(x))) * n_bins) + 1L
  pmin(n_bins, pmax(1L, as.integer(out)))
}

w <- as.data.table(qs2::qs_read(weights_path))
required <- c(metadata_cols, names(subtype_map), "Endothelial", "Mural cells", "Macrophages", "Microglial", "Monocytes", "Neurons")
missing <- setdiff(required, names(w))
if (length(missing)) stop("Missing columns in RCTD weights: ", paste(missing, collapse = ", "))

cell_cols <- setdiff(names(w), metadata_cols)
w[, dominant_rctd_label := dominant_label(.SD, cell_cols), .SDcols = cell_cols]
w[, vascular_dominant := dominant_rctd_label %in% c("Endothelial", "Mural cells")]
w[, TAM := Macrophages + Microglial + Monocytes]
w[, `NPC-P` := Subtype1]
w[, `OPC-M` := Subtype2]
w[, `MES-V` := Subtype3]
w[, `MES-I` := Subtype4]

distance_rows <- list()
slice_qc_rows <- list()
for (sl in sort(unique(w$slice))) {
  dt <- copy(w[slice == sl])
  vascular_coords <- as.matrix(dt[vascular_dominant == TRUE, .(x, y)])
  n_vasc <- nrow(vascular_coords)
  if (!n_vasc) {
    dt[, `:=`(
      distance_to_nearest_vascular = NA_real_,
      vascular_distance_bin = NA_integer_,
      vascular_distance_available = FALSE
    )]
  } else {
    coords <- as.matrix(dt[, .(x, y)])
    d <- RANN::nn2(data = vascular_coords, query = coords, k = 1)$nn.dists[, 1]
    dt[, `:=`(
      distance_to_nearest_vascular = as.numeric(d),
      vascular_distance_bin = rank_bins(d, 5L),
      vascular_distance_available = TRUE
    )]
  }
  distance_rows[[sl]] <- dt[, .(
    spot_id, slice, image, x, y,
    dominant_rctd_label, vascular_dominant,
    distance_to_nearest_vascular, vascular_distance_bin, vascular_distance_available
  )]
  slice_qc_rows[[sl]] <- data.table(
    slice = sl,
    n_spots = nrow(dt),
    vascular_dominant_spots = n_vasc,
    vascular_distance_available = n_vasc > 0,
    median_distance_to_vascular = median(dt$distance_to_nearest_vascular, na.rm = TRUE)
  )
}

spot_distance <- rbindlist(distance_rows, use.names = TRUE, fill = TRUE)
slice_qc <- rbindlist(slice_qc_rows, use.names = TRUE, fill = TRUE)

w2 <- merge(w, spot_distance[, .(spot_id, distance_to_nearest_vascular, vascular_distance_bin, vascular_distance_available)], by = "spot_id", all.x = TRUE)

group_cols <- c(
  "NPC-P", "OPC-M", "MES-V", "MES-I",
  "TAM", "Endothelial", "Mural cells", "Neurons",
  "Oligodendrocytes", "OPCs", "T cells", "B cells", "NK cells", "cDCs", "pDCs"
)
present_group_cols <- intersect(group_cols, names(w2))

long <- melt(
  w2[vascular_distance_available == TRUE],
  id.vars = c("spot_id", "slice", "vascular_distance_bin", "distance_to_nearest_vascular"),
  measure.vars = present_group_cols,
  variable.name = "group",
  value.name = "mean_weight"
)

per_slice_bin <- long[, .(
  n_spots = .N,
  mean_weight = mean(mean_weight, na.rm = TRUE),
  median_raw_distance = median(distance_to_nearest_vascular, na.rm = TRUE)
), by = .(slice, vascular_distance_bin, group)]

summary_by_bin <- per_slice_bin[, .(
  n_slices = .N,
  median_mean_weight = median(mean_weight, na.rm = TRUE),
  q25_mean_weight = quantile(mean_weight, 0.25, na.rm = TRUE),
  q75_mean_weight = quantile(mean_weight, 0.75, na.rm = TRUE)
), by = .(vascular_distance_bin, group)]

slopes <- per_slice_bin[, {
  fit <- lm(mean_weight ~ vascular_distance_bin)
  data.table(
    slope_per_bin = unname(coef(fit)[2]),
    near_minus_far = mean_weight[vascular_distance_bin == 1][1] - mean_weight[vascular_distance_bin == 5][1]
  )
}, by = .(slice, group)]

slope_summary <- slopes[, .(
  n_slices = .N,
  median_slope_per_bin = median(slope_per_bin, na.rm = TRUE),
  positive_slope_slices = sum(slope_per_bin > 0, na.rm = TRUE),
  negative_slope_slices = sum(slope_per_bin < 0, na.rm = TRUE),
  median_near_minus_far = median(near_minus_far, na.rm = TRUE),
  near_greater_than_far_slices = sum(near_minus_far > 0, na.rm = TRUE)
), by = group]

params <- data.table(
  parameter = c("input_weights", "vascular_definition", "distance_axis", "binning", "all_spots_retained", "statistics"),
  value = c(
    weights_path,
    "dominant RCTD label is Endothelial or Mural cells; independent of MES markers",
    "Euclidean distance from each spot to nearest vascular-dominant spot within the same slice",
    "within-slice rank bins, 1=nearest and 5=farthest; no spots excluded when vascular spots exist",
    "TRUE; slices with zero vascular-dominant spots have undefined distance and are reported as low-vascular caveat",
    "descriptive per-slice composition only; no relation-level p-value from gradient curve"
  )
)

fwrite(spot_distance, file.path(out_dir, "vascular_distance_per_spot.csv"))
fwrite(slice_qc, file.path(out_dir, "vascular_distance_slice_qc.csv"))
fwrite(per_slice_bin, file.path(out_dir, "vascular_distance_composition_perslice_bins.csv"))
fwrite(summary_by_bin, file.path(out_dir, "vascular_distance_composition_summary_bins.csv"))
fwrite(slopes, file.path(out_dir, "vascular_distance_group_slopes_perslice.csv"))
fwrite(slope_summary, file.path(out_dir, "vascular_distance_group_slopes_summary.csv"))
fwrite(params, file.path(out_dir, "vascular_distance_gradient_parameters_source.csv"))

qs2::qs_save(
  list(
    spot_distance = spot_distance,
    slice_qc = slice_qc,
    per_slice_bin = per_slice_bin,
    summary_by_bin = summary_by_bin,
    slopes = slopes,
    slope_summary = slope_summary,
    parameters = params
  ),
  file.path(out_dir, "vascular_distance_gradient_results.qs2")
)

plot_groups <- c("NPC-P", "OPC-M", "MES-V", "MES-I", "TAM", "Endothelial", "Mural cells", "Neurons")
plot_dt <- summary_by_bin[group %in% plot_groups]
p <- ggplot(plot_dt, aes(x = vascular_distance_bin, y = median_mean_weight, color = group, group = group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin = q25_mean_weight, ymax = q75_mean_weight, fill = group), alpha = 0.08, color = NA) +
  scale_x_continuous(breaks = 1:5, labels = c("near", "2", "3", "4", "far")) +
  labs(x = "Distance-to-vascular rank bin per slice", y = "Median RCTD weight across slices", color = NULL, fill = NULL) +
  theme_classic(base_size = 9)
ggsave(file.path(fig_dir, "距血管梯度_成分曲线.pdf"), p, width = 6.6, height = 4.2)
ggsave(file.path(fig_dir, "距血管梯度_成分曲线.png"), p, width = 6.6, height = 4.2, dpi = 300)

message("Wrote project 2 vascular distance gradient outputs to: ", out_dir)
