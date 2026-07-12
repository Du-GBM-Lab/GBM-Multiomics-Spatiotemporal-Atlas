#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-4 v3-B main figure candidate: crossK delta forest
# Scope:
#   Main response  : MES-lineage
#   Main statistic : crossK_count delta vs random-labeling null
#   Main cut       : top10
#   Predictors     : vascular vs neuron_control in the same view
#   Scales         : k=6 main, k=12 sensitivity
#
# Output:
#   PDF/PNG figure candidate + source data.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- getwd()
tab_dir <- file.path(base_dir, "tables/C3_4_local_niche")
fig_dir <- file.path(base_dir, "figures/C3_4_local_proximity")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

perslice_path <- file.path(tab_dir, "R9_C3_4_v3B_Ripley_sensitivity_perslice.csv")
summary_path <- file.path(tab_dir, "R9_C3_4_v3B_Ripley_sensitivity_summary.csv")

ps <- fread(perslice_path)
sm <- fread(summary_path)

main_ps <- ps[
  response == "MES-lineage" &
    predictor %in% c("vascular", "neuron_control") &
    top_rule == "top10" &
    k %in% c(6, 12)
]
main_sm <- sm[
  response == "MES-lineage" &
    predictor %in% c("vascular", "neuron_control") &
    top_rule == "top10" &
    k %in% c(6, 12) &
    sensitivity == "all_slices"
]

stopifnot(nrow(main_ps) == 18 * 2 * 2)
stopifnot(nrow(main_sm) == 4)
stopifnot(all(main_ps$count_preserved_all_perm))
stopifnot(all(main_sm$count_preserved_all))

slice_order <- main_ps[k == 6 & predictor == "vascular"][order(-delta_crossK)]$slice
y_levels <- c(rev(slice_order), "Summary")

plot_dt <- copy(main_ps)
plot_dt[, `:=`(
  k_label = paste0("k=", k),
  predictor_label = fifelse(predictor == "vascular", "Vascular", "Neuron control"),
  y_label = factor(slice, levels = y_levels),
  emp_sig = p_emp_crossK_greater < 0.05
)]

sum_dt <- copy(main_sm)
sum_dt[, `:=`(
  k_label = paste0("k=", k),
  predictor_label = fifelse(predictor == "vascular", "Vascular", "Neuron control"),
  y_label = factor("Summary", levels = y_levels)
)]

source_dt <- rbindlist(list(
  plot_dt[, .(
    row_type = "slice",
    slice, k, k_label, response, predictor, predictor_label, top_rule,
    top_fraction, delta_crossK, p_emp_crossK_greater, emp_sig,
    n_response_hot, n_predictor_hot, count_preserved_all_perm
  )],
  sum_dt[, .(
    row_type = "summary",
    slice = "Summary",
    k, k_label, response, predictor, predictor_label, top_rule,
    top_fraction, delta_crossK = median_delta_crossK,
    p_emp_crossK_greater = median_p_crossK,
    emp_sig = n_slice_crossK_p_lt_0p05 > 0,
    n_response_hot = median_response_hot,
    n_predictor_hot = median_predictor_hot,
    count_preserved_all_perm = count_preserved_all
  )]
), use.names = TRUE, fill = TRUE)

fwrite(source_dt, file.path(tab_dir, "R9_C3_4_v3B_main_crossK_forest_source.csv"))

pal <- c("Vascular" = "#0072B5", "Neuron control" = "#6B7280")

p <- ggplot() +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey45") +
  geom_segment(
    data = sum_dt,
    aes(x = q25_delta_crossK, xend = q75_delta_crossK, y = y_label, yend = y_label,
        color = predictor_label),
    linewidth = 1.1
  ) +
  geom_point(
    data = plot_dt,
    aes(x = delta_crossK, y = y_label, color = predictor_label, shape = emp_sig),
    size = 2.3, alpha = 0.88
  ) +
  geom_point(
    data = sum_dt,
    aes(x = median_delta_crossK, y = y_label, color = predictor_label),
    shape = 18, size = 4.2
  ) +
  facet_grid(k_label ~ predictor_label, scales = "free_x") +
  scale_color_manual(values = pal, guide = "none") +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("FALSE" = "empirical p >= 0.05", "TRUE" = "empirical p < 0.05"),
    name = "Per-slice random-labeling null"
  ) +
  labs(
    x = "crossK count delta vs random-labeling null",
    y = NULL,
    title = "MES-lineage high spots show local vascular proximity across Visium sections",
    subtitle = "Top 10%; random-labeling null. k=6: vascular median +0.458 (16/18 emp<0.05), neuron -0.219 (3/18). k=12: vascular +0.803 (15/18), neuron -0.497 (4/18).",
    caption = "Positive delta: more vascular-high spots near MES-lineage-high spots than random expectation.\nRandom-labeling fixed spot positions and point counts; class labels were permuted only. Spot-level proximity only."
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 8.5),
    plot.caption = element_text(size = 7, hjust = 0),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom",
    plot.margin = margin(8, 8, 8, 8)
  )

pdf_file <- file.path(fig_dir, "R9_C3_4_v3B_crossK_forest_MESlineage_vascular_vs_neuron.pdf")
png_file <- file.path(fig_dir, "R9_C3_4_v3B_crossK_forest_MESlineage_vascular_vs_neuron.png")
ggsave(pdf_file, p, width = 9.5, height = 7.2, device = cairo_pdf)
ggsave(png_file, p, width = 9.5, height = 7.2, dpi = 300)

cat("== figure written:\n", pdf_file, "\n", png_file, "\n")
cat("== source data written:\n", file.path(tab_dir, "R9_C3_4_v3B_main_crossK_forest_source.csv"), "\n")
