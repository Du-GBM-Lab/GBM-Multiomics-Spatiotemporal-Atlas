#!/usr/bin/env Rscript
# =============================================================================
# R9 | IVY ROI-like spatial domains + malignant subtype composition
# Purpose:
#   Build continuous, smoothed IVY ROI-like regions from IVY signature scores,
#   then summarize RCTD malignant subtype composition inside each ROI.
#
# Boundaries:
#   - This is descriptive integration, not a new proximity/statistical test.
#   - ROI labels are signature-based "-like" spatial domains, not recovered
#     categorical IVY-GAP Location annotations.
#   - Do not infer contact, recruitment, interaction, hypoxia causality, or
#     mechanism from this script.
#
# STOP:
#   Report ROI segmentation audit + ROI composition matrix + preview figures.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
  library(ggplot2)
})

set.seed(1)

base_dir <- getwd()
score_path <- file.path(base_dir, "tables/R9_stage3_signature_scores/R9_stage3_IVY_hypoxia_scores_per_spot.csv")
weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
out_dir <- file.path(base_dir, "tables/R9_IVY_ROI_composition")
fig_dir <- file.path(base_dir, "figures/R9_IVY_ROI_composition")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

ivy_cols <- c("IVY_CT", "IVY_CTmvp", "IVY_CTpan", "IVY_LE")
roi_levels <- c("CT-like", "MVP-like", "PAN-like", "LE-like")
names(roi_levels) <- ivy_cols

subtype_cols <- c("NPC-P", "OPC-M", "MES-V", "MES-I")
subtype_palette <- c(
  "NPC-P" = "#4E79A7",
  "OPC-M" = "#59A14F",
  "MES-V" = "#E15759",
  "MES-I" = "#F28E2B"
)
roi_palette <- c(
  "CT-like" = "#2F6DB3",
  "MVP-like" = "#E66B00",
  "PAN-like" = "#6A3D9A",
  "LE-like" = "#009B77"
)

## ---- Helpers ---------------------------------------------------------------
scale_vec <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

row_argmax <- function(mat) {
  colnames(mat)[max.col(mat, ties.method = "first")]
}

mode_label <- function(x, current) {
  tab <- sort(table(x), decreasing = TRUE)
  if (!length(tab)) return(current)
  top <- names(tab)[tab == max(tab)]
  if (current %in% top) return(current)
  top[1]
}

build_neighbors <- function(coord, k = 6L, radius_multiplier = 1.5) {
  D <- as.matrix(dist(coord))
  diag(D) <- Inf
  kth <- apply(D, 1, function(z) sort(z, partial = k)[k])
  radius <- median(kth, na.rm = TRUE) * radius_multiplier
  neigh <- lapply(seq_len(nrow(D)), function(i) which(D[i, ] <= radius))
  list(neigh = neigh, radius = radius, realized_k = lengths(neigh))
}

smooth_scores <- function(mat, neigh) {
  out <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  colnames(out) <- colnames(mat)
  for (i in seq_len(nrow(mat))) {
    idx <- unique(c(i, neigh[[i]]))
    out[i, ] <- colMeans(mat[idx, , drop = FALSE], na.rm = TRUE)
  }
  out
}

majority_smooth_labels <- function(labels, smooth_mat, neigh, n_iter = 2L) {
  labs <- labels
  for (iter in seq_len(n_iter)) {
    old <- labs
    for (i in seq_along(old)) {
      idx <- unique(c(i, neigh[[i]]))
      labs[i] <- mode_label(old[idx], old[i])
    }
    ## In the rare case majority smoothing creates an NA-like issue, fall back to
    ## the highest smoothed signature for that spot.
    bad <- is.na(labs) | !nzchar(labs)
    if (any(bad)) labs[bad] <- row_argmax(smooth_mat[bad, , drop = FALSE])
  }
  labs
}

connected_components_one <- function(labels, neigh, target_label) {
  idx <- which(labels == target_label)
  if (!length(idx)) {
    return(data.table(n_components = 0L, largest_component_n = 0L,
                      largest_component_fraction = NA_real_))
  }
  idx_set <- rep(FALSE, length(labels))
  idx_set[idx] <- TRUE
  visited <- rep(FALSE, length(labels))
  comp_sizes <- integer(0)
  for (start in idx) {
    if (visited[start]) next
    queue <- start
    visited[start] <- TRUE
    size <- 0L
    while (length(queue)) {
      v <- queue[1]
      queue <- queue[-1]
      size <- size + 1L
      nb <- neigh[[v]]
      nb <- nb[idx_set[nb] & !visited[nb]]
      if (length(nb)) {
        visited[nb] <- TRUE
        queue <- c(queue, nb)
      }
    }
    comp_sizes <- c(comp_sizes, size)
  }
  data.table(
    n_components = length(comp_sizes),
    largest_component_n = max(comp_sizes),
    largest_component_fraction = max(comp_sizes) / length(idx)
  )
}

safe_name <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

## ---- Load inputs -----------------------------------------------------------
scores <- fread(score_path)
weights <- as.data.table(qs2::qs_read(weights_path))

needed_scores <- c("spot_id", "slice", "image", "x", "y", ivy_cols,
                   "Hypoxia_Buffa", "vascular")
stopifnot(all(needed_scores %in% names(scores)))
stopifnot(all(c("spot_id", "Subtype1", "Subtype2", "Subtype3", "Subtype4",
                "Endothelial", "Mural cells") %in% names(weights)))

feat <- merge(
  scores[, ..needed_scores],
  weights[, .(
    spot_id,
    `NPC-P` = Subtype1,
    `OPC-M` = Subtype2,
    `MES-V` = Subtype3,
    `MES-I` = Subtype4,
    `MES-lineage` = Subtype3 + Subtype4,
    malignant_total = Subtype1 + Subtype2 + Subtype3 + Subtype4,
    vascular_RCTD = Endothelial + `Mural cells`
  )],
  by = "spot_id",
  all.x = TRUE
)
stopifnot(nrow(feat) == nrow(scores))

## ---- Build smoothed ROI-like domains --------------------------------------
message("== Building smoothed IVY ROI-like domains.")
roi_list <- list()
audit_list <- list()
component_list <- list()

for (sl in sort(unique(feat$slice))) {
  dt <- copy(feat[slice == sl])
  setorder(dt, spot_id)
  coord <- as.matrix(dt[, .(x, y)])
  nb <- build_neighbors(coord, k = 6L, radius_multiplier = 1.5)

  z <- as.matrix(dt[, lapply(.SD, scale_vec), .SDcols = ivy_cols])
  colnames(z) <- ivy_cols
  raw_argmax <- row_argmax(z)

  sm <- smooth_scores(z, nb$neigh)
  smooth_argmax <- row_argmax(sm)
  roi <- majority_smooth_labels(smooth_argmax, sm, nb$neigh, n_iter = 2L)
  roi <- unname(roi_levels[roi])

  ord <- t(apply(sm, 1, sort, decreasing = TRUE))
  margin <- ord[, 1] - ord[, 2]

  out <- data.table(
    spot_id = dt$spot_id,
    slice = dt$slice,
    image = dt$image,
    x = dt$x,
    y = dt$y,
    raw_argmax = unname(roi_levels[raw_argmax]),
    smooth_argmax = unname(roi_levels[smooth_argmax]),
    IVY_ROI = factor(roi, levels = roi_palette |> names()),
    smooth_margin = margin,
    k6_radius = nb$radius,
    realized_k = nb$realized_k
  )
  out <- cbind(out, as.data.table(sm))
  setnames(out, ivy_cols, paste0(ivy_cols, "_smooth_z"))
  roi_list[[sl]] <- out

  audit_list[[sl]] <- data.table(
    slice = sl,
    n_spots = nrow(dt),
    k6_radius = nb$radius,
    realized_k_median = median(nb$realized_k),
    realized_k_min = min(nb$realized_k),
    realized_k_max = max(nb$realized_k),
    zero_neighbor_fraction = mean(nb$realized_k == 0),
    raw_vs_smoothed_roi_agreement = mean(out$raw_argmax == out$IVY_ROI),
    median_smooth_margin = median(out$smooth_margin, na.rm = TRUE)
  )

  comp <- rbindlist(lapply(names(roi_palette), function(r) {
    cbind(data.table(slice = sl, IVY_ROI = r),
          connected_components_one(as.character(out$IVY_ROI), nb$neigh, r))
  }))
  component_list[[sl]] <- comp
  message("  ROI done slice=", sl, " n=", nrow(dt), " radius=", round(nb$radius, 2))
}

roi_assign <- rbindlist(roi_list, use.names = TRUE)
roi_audit <- rbindlist(audit_list, use.names = TRUE)
roi_components <- rbindlist(component_list, use.names = TRUE)

feat_roi <- merge(feat, roi_assign[, .(spot_id, raw_argmax, smooth_argmax, IVY_ROI,
                                       smooth_margin, k6_radius, realized_k,
                                       IVY_CT_smooth_z, IVY_CTmvp_smooth_z,
                                       IVY_CTpan_smooth_z, IVY_LE_smooth_z)],
                  by = "spot_id", all.x = TRUE)

## ---- Composition summaries -------------------------------------------------
roi_area_slice <- feat_roi[, .(
  n_spots = .N,
  spot_fraction_in_slice = .N / feat_roi[slice == .BY$slice, .N]
), by = .(slice, IVY_ROI)]

roi_area_global <- feat_roi[, .(
  n_spots = .N,
  n_slices_present = uniqueN(slice),
  spot_fraction_global = .N / nrow(feat_roi),
  median_smooth_margin = median(smooth_margin, na.rm = TRUE)
), by = IVY_ROI]
setorder(roi_area_global, IVY_ROI)

abs_cols <- c(subtype_cols, "MES-lineage", "malignant_total", "vascular_RCTD",
              "vascular", "Hypoxia_Buffa")
roi_abs_slice <- feat_roi[, lapply(.SD, mean, na.rm = TRUE),
                          by = .(slice, IVY_ROI), .SDcols = abs_cols]
roi_abs_global <- feat_roi[, lapply(.SD, mean, na.rm = TRUE),
                           by = IVY_ROI, .SDcols = abs_cols]

roi_malignant_comp <- feat_roi[, {
  sums <- colSums(.SD, na.rm = TRUE)
  total <- sum(sums)
  as.list(sums / total)
}, by = IVY_ROI, .SDcols = subtype_cols]
roi_malignant_comp[, malignant_total_weight_sum := feat_roi[, sum(malignant_total, na.rm = TRUE), by = IVY_ROI]$V1]

abs_long <- melt(
  roi_abs_global[, c("IVY_ROI", subtype_cols, "MES-lineage", "malignant_total"), with = FALSE],
  id.vars = "IVY_ROI",
  variable.name = "feature",
  value.name = "mean_RCTD_weight"
)
comp_long <- melt(
  roi_malignant_comp[, c("IVY_ROI", subtype_cols), with = FALSE],
  id.vars = "IVY_ROI",
  variable.name = "subtype",
  value.name = "fraction_within_malignant"
)

subtype_roi_distribution <- feat_roi[, lapply(.SD, sum, na.rm = TRUE),
                                     by = IVY_ROI, .SDcols = subtype_cols]
subtype_roi_distribution_long <- melt(
  subtype_roi_distribution,
  id.vars = "IVY_ROI",
  variable.name = "subtype",
  value.name = "subtype_weight_sum"
)
subtype_roi_distribution_long[, fraction_of_subtype_total :=
                                subtype_weight_sum / sum(subtype_weight_sum, na.rm = TRUE),
                              by = subtype]

## ---- Output tables ---------------------------------------------------------
fwrite(roi_assign, file.path(out_dir, "R9_IVY_ROI_smoothed_domain_assignment.csv"))
fwrite(roi_audit, file.path(out_dir, "R9_IVY_ROI_segmentation_audit_by_slice.csv"))
fwrite(roi_components, file.path(out_dir, "R9_IVY_ROI_connected_components_audit.csv"))
fwrite(roi_area_slice, file.path(out_dir, "R9_IVY_ROI_area_by_slice.csv"))
fwrite(roi_area_global, file.path(out_dir, "R9_IVY_ROI_area_global.csv"))
fwrite(roi_abs_slice, file.path(out_dir, "R9_IVY_ROI_composition_absolute_by_slice.csv"))
fwrite(roi_abs_global, file.path(out_dir, "R9_IVY_ROI_composition_absolute_global.csv"))
fwrite(roi_malignant_comp, file.path(out_dir, "R9_IVY_ROI_malignant_subtype_composition_global.csv"))
fwrite(abs_long, file.path(out_dir, "R9_IVY_ROI_composition_absolute_global_long.csv"))
fwrite(comp_long, file.path(out_dir, "R9_IVY_ROI_malignant_subtype_composition_global_long.csv"))
fwrite(subtype_roi_distribution_long,
       file.path(out_dir, "R9_IVY_ROI_subtype_across_ROI_distribution_long.csv"))

params <- data.table(
  parameter = c("roi_input_scores", "roi_score_transform", "roi_smoothing",
                "roi_label_rule", "majority_smoothing_iterations",
                "k_for_spatial_smoothing", "radius_rule",
                "composition_input", "interpretation_boundary"),
  value = c("IVY_CT, IVY_CTmvp, IVY_CTpan, IVY_LE AddModuleScore scores",
            "within-slice z-score per signature before smoothing",
            "mean of focal spot plus neighbors within per-slice k6 radius",
            "highest smoothed IVY z-score after two local majority passes",
            "2",
            "6",
            "per-slice median 6th-neighbor distance x 1.5",
            "RCTD full-mode spot weights",
            "descriptive ROI-like region composition; not a new proximity test or recovered Location annotation")
)
fwrite(params, file.path(out_dir, "R9_IVY_ROI_method_parameters.csv"))

caption_notes <- data.table(
  item = c("figure_role", "roi_definition", "fragmentation_caveat",
           "subtype_scope", "y_axis", "allowed_claim", "forbidden_claim"),
  note = c(
    "Candidate main-text spatial landscape overview; descriptive integration, not a new proximity test.",
    "ROI-like domains are defined by spatially smoothed IVY signature scores; they are not recovered categorical IVY-GAP Location labels.",
    "MVP-like ROI-like domains are spatially fragmented in the segmentation audit; this must be stated if the map is shown.",
    "The current malignant subtype system contains NPC-P, OPC-M, MES-V and MES-I; there is no separate AC-like/Astro-like subtype in the revised classification.",
    "Main bar plot reports mean malignant subtype RCTD weights within each ROI-like domain; these values are not proportions and do not sum to one because non-malignant weights are not shown.",
    "MVP-like ROI-like domains show the highest descriptive MES-lineage/MES-V mean weight, consistent with the vascular/MVP-like proximity result.",
    "Do not write recovered Location, histology-grade anatomical regions, physical contact, recruitment, interaction, or hypoxia/PAN causality."
  )
)
fwrite(caption_notes, file.path(out_dir, "R9_IVY_ROI_figure_caption_notes.csv"))

## ---- Figures ---------------------------------------------------------------
theme_r9 <- theme_bw(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

map_df <- copy(feat_roi)
map_df[, IVY_ROI := factor(IVY_ROI, levels = names(roi_palette))]

p_map <- ggplot(map_df, aes(x = x, y = y, color = IVY_ROI)) +
  geom_point(size = 0.18, alpha = 0.95) +
  scale_color_manual(values = roi_palette, drop = FALSE) +
  scale_y_reverse() +
  facet_wrap(~ slice, scales = "free", ncol = 6) +
  labs(
    title = "Smoothed IVY signature-based ROI-like domains",
    subtitle = "Within-slice z-scored IVY scores, k=6 spatial smoothing, local majority refinement",
    x = NULL, y = NULL, color = "ROI-like domain"
  ) +
  theme_r9 +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "R9_IVY_ROI_smoothed_domain_maps_all_slices.pdf"),
       p_map, width = 11, height = 8.2, useDingbats = FALSE)
ggsave(file.path(fig_dir, "R9_IVY_ROI_smoothed_domain_maps_all_slices.png"),
       p_map, width = 11, height = 8.2, dpi = 300)

p_map_supp <- ggplot(map_df, aes(x = x, y = y, color = IVY_ROI)) +
  geom_point(size = 0.34, alpha = 1) +
  scale_color_manual(values = roi_palette, drop = FALSE) +
  scale_y_reverse() +
  facet_wrap(~ slice, scales = "free", ncol = 3) +
  labs(
    x = NULL, y = NULL, color = "ROI-like domain"
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey96", color = "grey65", linewidth = 0.25),
    strip.text = element_text(size = 7.5, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "grey45", linewidth = 0.25),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.margin = margin(4, 4, 4, 4)
  )
ggsave(file.path(fig_dir, "R9_IVY_ROI_smoothed_domain_maps_all_slices_supplement.pdf"),
       p_map_supp, width = 8.8, height = 14.2, useDingbats = FALSE)
ggsave(file.path(fig_dir, "R9_IVY_ROI_smoothed_domain_maps_all_slices_supplement.png"),
       p_map_supp, width = 8.8, height = 14.2, dpi = 300)

p_abs <- ggplot(abs_long[feature %in% subtype_cols],
                aes(x = IVY_ROI, y = mean_RCTD_weight, fill = feature)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = subtype_palette, drop = FALSE) +
  labs(
    title = "Mean malignant subtype RCTD weights by IVY ROI-like domain",
    subtitle = "Descriptive mean weights from RCTD; not proportions and not recovered IVY-GAP Location",
    x = NULL, y = "Mean RCTD weight", fill = "Subtype"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(fig_dir, "R9_IVY_ROI_absolute_subtype_composition_bar.pdf"),
       p_abs, width = 6.4, height = 4.2, useDingbats = FALSE)
ggsave(file.path(fig_dir, "R9_IVY_ROI_absolute_subtype_composition_bar.png"),
       p_abs, width = 6.4, height = 4.2, dpi = 300)

p_abs_main <- ggplot(abs_long[feature %in% subtype_cols],
                     aes(x = IVY_ROI, y = mean_RCTD_weight, fill = feature)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = subtype_palette, drop = FALSE) +
  labs(
    x = NULL,
    y = "Mean malignant subtype RCTD weight",
    fill = "Subtype"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title.y = element_text(size = 8)
  )
ggsave(file.path(fig_dir, "R9_IVY_ROI_mean_subtype_weights_main_candidate.pdf"),
       p_abs_main, width = 4.8, height = 3.4, useDingbats = FALSE)
ggsave(file.path(fig_dir, "R9_IVY_ROI_mean_subtype_weights_main_candidate.png"),
       p_abs_main, width = 4.8, height = 3.4, dpi = 300)

p_subtype_dist <- ggplot(subtype_roi_distribution_long,
                         aes(x = subtype, y = fraction_of_subtype_total, fill = IVY_ROI)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = roi_palette, drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  labs(
    title = "Distribution of each malignant subtype across IVY ROI-like domains",
    subtitle = "Each bar sums to 100% of that subtype's total RCTD weight",
    x = NULL, y = "Fraction of subtype weight", fill = "ROI-like domain"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(fig_dir, "R9_IVY_ROI_subtype_across_ROI_distribution.pdf"),
       p_subtype_dist, width = 6.2, height = 4.2, useDingbats = FALSE)
ggsave(file.path(fig_dir, "R9_IVY_ROI_subtype_across_ROI_distribution.png"),
       p_subtype_dist, width = 6.2, height = 4.2, dpi = 300)

p_comp <- ggplot(comp_long, aes(x = IVY_ROI, y = fraction_within_malignant, fill = subtype)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = subtype_palette, drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  labs(
    title = "Malignant-subtype composition within IVY ROI-like domains",
    subtitle = "Subtype fractions normalized within the malignant RCTD weight in each ROI",
    x = NULL, y = "Fraction within malignant weight", fill = "Subtype"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(fig_dir, "R9_IVY_ROI_malignant_subtype_composition_bar.pdf"),
       p_comp, width = 6.4, height = 4.2, useDingbats = FALSE)
ggsave(file.path(fig_dir, "R9_IVY_ROI_malignant_subtype_composition_bar.png"),
       p_comp, width = 6.4, height = 4.2, dpi = 300)

heat_dt <- abs_long[feature %in% c(subtype_cols, "MES-lineage", "malignant_total")]
p_heat <- ggplot(heat_dt, aes(x = IVY_ROI, y = feature, fill = mean_RCTD_weight)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "white", high = "#B2182B") +
  labs(
    title = "ROI composition heatmap",
    subtitle = "Mean RCTD weight per ROI-like domain",
    x = NULL, y = NULL, fill = "Mean weight"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank())
ggsave(file.path(fig_dir, "R9_IVY_ROI_composition_heatmap.pdf"),
       p_heat, width = 6.2, height = 3.8, useDingbats = FALSE)
ggsave(file.path(fig_dir, "R9_IVY_ROI_composition_heatmap.png"),
       p_heat, width = 6.2, height = 3.8, dpi = 300)

p_area <- ggplot(roi_area_global, aes(x = IVY_ROI, y = spot_fraction_global, fill = IVY_ROI)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = roi_palette, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  labs(
    title = "Global area fraction of smoothed IVY ROI-like domains",
    x = NULL, y = "Spot fraction"
  ) +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "R9_IVY_ROI_area_fraction.pdf"),
       p_area, width = 5.6, height = 3.6, useDingbats = FALSE)
ggsave(file.path(fig_dir, "R9_IVY_ROI_area_fraction.png"),
       p_area, width = 5.6, height = 3.6, dpi = 300)

## ---- STOP report -----------------------------------------------------------
cat("\n[STOP IVY ROI]\n")
cat("ROI assignment spots:", nrow(roi_assign), "\n")
cat("Segmentation audit: median realized k =",
    median(roi_audit$realized_k_median), "; max zero-neighbor fraction =",
    round(max(roi_audit$zero_neighbor_fraction), 4), "\n")
cat("Raw argmax vs smoothed ROI agreement: median =",
    round(median(roi_audit$raw_vs_smoothed_roi_agreement), 3), "\n")
cat("\nROI area global:\n")
print(roi_area_global)
cat("\nAbsolute subtype composition by ROI:\n")
print(roi_abs_global[, c("IVY_ROI", subtype_cols, "MES-lineage", "malignant_total", "vascular_RCTD"), with = FALSE])
cat("\nMalignant-normalized subtype composition by ROI:\n")
print(roi_malignant_comp)
cat("\nOutputs:\n")
cat("  tables:", out_dir, "\n")
cat("  figures:", fig_dir, "\n")
cat("\nInterpretation boundary: descriptive ROI-like composition only; do not treat these labels as recovered IVY Location or a new statistical test.\n")
