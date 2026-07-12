#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-4 v3-B supplementary/secondary panels
# Purpose:
#   1) Nearest-distance orthogonal support for MES-lineage.
#   2) MES-V / MES-I secondary crossK panel.
#   3) C3-3 continuous-association baseline panel.
#
# No new statistics are computed. All panels read finalized source tables.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- getwd()
tab_c34 <- file.path(base_dir, "tables/C3_4_local_niche")
tab_c33 <- file.path(base_dir, "tables/C3_niche_preflight")
fig_dir <- file.path(base_dir, "figures/C3_4_local_proximity")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

perslice <- fread(file.path(tab_c34, "R9_C3_4_v3B_Ripley_sensitivity_perslice.csv"))
summ <- fread(file.path(tab_c34, "R9_C3_4_v3B_Ripley_sensitivity_summary.csv"))

top_lab <- c(top05 = "Top 5%", top10 = "Top 10%", top20 = "Top 20%")
pred_lab <- c(vascular = "Vascular", neuron_control = "Neuron control")
resp_lab <- c("MES-lineage" = "MES-lineage", "MES-V" = "MES-V", "MES-I" = "MES-I")

base_theme <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5),
    plot.caption = element_text(size = 7, hjust = 0),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text = element_text(size = 8),
    legend.position = "bottom",
    plot.margin = margin(8, 8, 8, 8)
  )

write_plot <- function(plot, stub, width = 9.5, height = 6.2) {
  pdf_file <- file.path(fig_dir, paste0(stub, ".pdf"))
  png_file <- file.path(fig_dir, paste0(stub, ".png"))
  ggsave(pdf_file, plot, width = width, height = height, device = cairo_pdf)
  ggsave(png_file, plot, width = width, height = height, dpi = 300)
  cat("== figure written:\n", pdf_file, "\n", png_file, "\n")
}

## ---- 1) nearest-distance orthogonal support -------------------------------
nn_ps <- perslice[
  response == "MES-lineage" &
    predictor %in% c("vascular", "neuron_control") &
    top_rule %in% c("top05", "top10", "top20") &
    k %in% c(6, 12)
]
nn_sm <- summ[
  response == "MES-lineage" &
    predictor %in% c("vascular", "neuron_control") &
    top_rule %in% c("top05", "top10", "top20") &
    k %in% c(6, 12) &
    sensitivity == "all_slices"
]
stopifnot(nrow(nn_ps) == 18 * 2 * 3 * 2)
stopifnot(nrow(nn_sm) == 12)

nn_ps[, `:=`(
  k_label = paste0("k=", k),
  predictor_label = factor(pred_lab[predictor], levels = pred_lab),
  top_label = factor(top_lab[top_rule], levels = top_lab),
  emp_sig = p_emp_nn_less < 0.05
)]
nn_sm[, `:=`(
  k_label = paste0("k=", k),
  predictor_label = factor(pred_lab[predictor], levels = pred_lab),
  top_label = factor(top_lab[top_rule], levels = top_lab)
)]

nn_source <- rbindlist(list(
  nn_ps[, .(
    row_type = "slice", slice, k, k_label, response, predictor, predictor_label,
    top_rule, top_label, top_fraction, delta_nn_closer, p_emp_nn_less, emp_sig,
    n_response_hot, n_predictor_hot, count_preserved_all_perm
  )],
  nn_sm[, .(
    row_type = "summary", slice = "Summary", k, k_label, response, predictor,
    predictor_label, top_rule, top_label, top_fraction,
    delta_nn_closer = median_delta_nn_closer, p_emp_nn_less = median_p_nn,
    emp_sig = n_slice_nn_p_lt_0p05 > 0,
    n_response_hot = median_response_hot, n_predictor_hot = median_predictor_hot,
    count_preserved_all_perm = count_preserved_all
  )]
), use.names = TRUE, fill = TRUE)
fwrite(nn_source, file.path(tab_c34, "R9_C3_4_v3B_nearest_distance_orthogonal_source.csv"))

p_nn <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey45") +
  geom_point(
    data = nn_ps,
    aes(x = delta_nn_closer, y = top_label, shape = emp_sig),
    position = position_jitter(height = 0.08, width = 0),
    size = 1.8, alpha = 0.72, color = "grey35"
  ) +
  geom_segment(
    data = nn_sm,
    aes(x = q25_delta_nn_closer, xend = q75_delta_nn_closer,
        y = top_label, yend = top_label),
    linewidth = 1.0, color = "#b2182b"
  ) +
  geom_point(
    data = nn_sm,
    aes(x = median_delta_nn_closer, y = top_label),
    shape = 18, size = 3.7, color = "#b2182b"
  ) +
  facet_grid(k_label ~ predictor_label, scales = "free_x") +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("FALSE" = "empirical p >= 0.05", "TRUE" = "empirical p < 0.05"),
    name = "Per-slice random-labeling null"
  ) +
  labs(
    x = "nearest-distance delta vs random-labeling null",
    y = NULL,
    title = "Nearest-distance supports local vascular proximity as a conservative orthogonal metric",
    subtitle = "MES-lineage only; nearest-distance is secondary to crossK and is most stable for top10/top20.",
    caption = "Positive delta: MES-lineage-high spots are closer to predictor-high spots than random expectation.\nThis conservative metric does not independently carry MES-V/MES-I or top5 claims; primary inference remains crossK."
  ) +
  base_theme
write_plot(p_nn, "R9_C3_4_v3B_nearest_distance_orthogonal_MESlineage", width = 9.5, height = 6.0)

## ---- 2) MES-V / MES-I secondary crossK panel ------------------------------
sub_ps <- perslice[
  response %in% c("MES-V", "MES-I") &
    predictor %in% c("vascular", "neuron_control") &
    top_rule == "top10" &
    k %in% c(6, 12)
]
sub_sm <- summ[
  response %in% c("MES-V", "MES-I") &
    predictor %in% c("vascular", "neuron_control") &
    top_rule == "top10" &
    k %in% c(6, 12) &
    sensitivity == "all_slices"
]
stopifnot(nrow(sub_ps) == 18 * 2 * 2 * 2)
stopifnot(nrow(sub_sm) == 8)

sub_ps[, `:=`(
  k_label = paste0("k=", k),
  predictor_label = factor(pred_lab[predictor], levels = pred_lab),
  response_label = factor(resp_lab[response], levels = c("MES-V", "MES-I")),
  emp_sig = p_emp_crossK_greater < 0.05
)]
sub_sm[, `:=`(
  k_label = paste0("k=", k),
  predictor_label = factor(pred_lab[predictor], levels = pred_lab),
  response_label = factor(resp_lab[response], levels = c("MES-V", "MES-I"))
)]

sub_source <- rbindlist(list(
  sub_ps[, .(
    row_type = "slice", slice, k, k_label, response, response_label, predictor,
    predictor_label, top_rule, top_fraction, delta_crossK,
    p_emp_crossK_greater, emp_sig, n_response_hot, n_predictor_hot,
    count_preserved_all_perm
  )],
  sub_sm[, .(
    row_type = "summary", slice = "Summary", k, k_label, response, response_label,
    predictor, predictor_label, top_rule, top_fraction,
    delta_crossK = median_delta_crossK, p_emp_crossK_greater = median_p_crossK,
    emp_sig = n_slice_crossK_p_lt_0p05 > 0,
    n_response_hot = median_response_hot, n_predictor_hot = median_predictor_hot,
    count_preserved_all_perm = count_preserved_all
  )]
), use.names = TRUE, fill = TRUE)
fwrite(sub_source, file.path(tab_c34, "R9_C3_4_v3B_secondary_MESV_MESI_crossK_source.csv"))

p_sub <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey45") +
  geom_point(
    data = sub_ps,
    aes(x = delta_crossK, y = response_label, shape = emp_sig),
    position = position_jitter(height = 0.08, width = 0),
    size = 1.9, alpha = 0.72, color = "grey35"
  ) +
  geom_segment(
    data = sub_sm,
    aes(x = q25_delta_crossK, xend = q75_delta_crossK,
        y = response_label, yend = response_label),
    linewidth = 1.0, color = "#0072B5"
  ) +
  geom_point(
    data = sub_sm,
    aes(x = median_delta_crossK, y = response_label),
    shape = 18, size = 3.7, color = "#0072B5"
  ) +
  facet_grid(k_label ~ predictor_label, scales = "free_x") +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("FALSE" = "empirical p >= 0.05", "TRUE" = "empirical p < 0.05"),
    name = "Per-slice random-labeling null"
  ) +
  labs(
    x = "crossK count delta vs random-labeling null",
    y = NULL,
    title = "MES-V and MES-I secondary proximity signals",
    subtitle = "Top 10%; subtype-level panels provide context only, with MES-lineage retained as the primary spatial result.",
    caption = "Subtype panels are secondary because MES-V/MES-I weights are low-abundance at Visium spot resolution.\nUse these panels as subtype context only; MES-lineage crossK remains the primary spatial result."
  ) +
  base_theme
write_plot(p_sub, "R9_C3_4_v3B_secondary_MESV_MESI_crossK", width = 9.5, height = 5.3)

## ---- 3) C3-3 continuous-association baseline ------------------------------
c33 <- fread(file.path(tab_c33, "R9_C3_3_forest_k6_vascular_vs_neuron_control_source.csv"))
c33 <- c33[
  response %in% c("MES-lineage", "MES-V", "MES-I") &
    predictor %in% c("vascular_niche", "neuron_control_niche")
]
c33[, predictor_label := fifelse(predictor == "vascular_niche", "Vascular-neighborhood", "Neuron-control neighborhood")]
c33[, predictor_label := factor(predictor_label, levels = c("Vascular-neighborhood", "Neuron-control neighborhood"))]
c33[, response_label := factor(response, levels = c("MES-lineage", "MES-V", "MES-I"))]
c33_slice <- c33[slice != "Summary median"]
c33_sum <- c33[slice == "Summary median"]

baseline_source <- copy(c33)
fwrite(baseline_source, file.path(tab_c34, "R9_C3_4_baseline_C3_3_continuous_association_source.csv"))

p_base <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey45") +
  geom_point(
    data = c33_slice,
    aes(x = observed_rho, y = response_label,
        shape = vector_shift_emp_lt_0p05 == "empirical p < 0.05"),
    position = position_jitter(height = 0.08, width = 0),
    size = 1.8, alpha = 0.72, color = "grey35"
  ) +
  geom_segment(
    data = c33_sum,
    aes(x = xmin, xend = xmax, y = response_label, yend = response_label),
    linewidth = 1.0, color = "#b2182b"
  ) +
  geom_point(
    data = c33_sum,
    aes(x = observed_rho, y = response_label),
    shape = 18, size = 3.7, color = "#b2182b"
  ) +
  facet_wrap(~ predictor_label, nrow = 1, scales = "free_x") +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("FALSE" = "empirical p >= 0.05", "TRUE" = "empirical p < 0.05"),
    name = "Per-slice vector-shift null"
  ) +
  labs(
    x = "per-slice Spearman rho",
    y = NULL,
    title = "Continuous neighborhood association provides baseline spatial context",
    subtitle = "C3-3 k=6 continuous association; retained as background because v3-B proximity is the primary statistic.",
    caption = "This panel shows full-slice continuous association only. It supports directional consistency but does not carry the local proximity claim."
  ) +
  base_theme
write_plot(p_base, "R9_C3_4_baseline_C3_3_continuous_association", width = 9.5, height = 4.8)

cat("== source data written:\n",
    file.path(tab_c34, "R9_C3_4_v3B_nearest_distance_orthogonal_source.csv"), "\n",
    file.path(tab_c34, "R9_C3_4_v3B_secondary_MESV_MESI_crossK_source.csv"), "\n",
    file.path(tab_c34, "R9_C3_4_baseline_C3_3_continuous_association_source.csv"), "\n")
