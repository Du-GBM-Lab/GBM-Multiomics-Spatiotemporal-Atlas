#!/usr/bin/env Rscript
# =============================================================================
# R9 | C3-3: forest plots + representative spatial panels for C3 association
# Input : C3-2 per-slice rho, C3-3-stat summary, C3-1 spot-level source table.
# Output: figure candidates + source data. No new statistics.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

set.seed(1)

base_dir <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
perslice_path <- file.path(base_dir, "tables/R9_C3_2_perslice_vector_shift_association.csv")
summary_path <- file.path(base_dir, "tables/R9_C3_3stat_slice_level_summary_focus_pairs.csv")
score_path <- file.path(base_dir, "tables/R9_C3_1_neighborhood_niche_scores.csv")
out_dir <- file.path(base_dir, "outputs/R9_C3_figures")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(base_dir, "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

perslice <- fread(perslice_path)
summ <- fread(summary_path)
scores <- fread(score_path)

main_responses <- c("MES-lineage", "MES-V", "MES-I")
main_predictors <- c("vascular_niche", "neuron_control_niche")
supp_predictors <- c("vascular_niche", "myeloid_niche", "neuron_control_niche")

## ---- Representative slice selection: close to median + visible signals ------
sel_pool <- perslice[k == 6 & response == "MES-lineage" & predictor == "vascular_niche"]
med_rho <- median(sel_pool$observed_rho, na.rm = TRUE)

vis_q <- scores[k == 6, .(
  mes_lineage_mean = mean(`MES-lineage`, na.rm = TRUE),
  mes_lineage_q90 = quantile(`MES-lineage`, 0.90, na.rm = TRUE),
  mes_lineage_q95 = quantile(`MES-lineage`, 0.95, na.rm = TRUE),
  mes_lineage_q99 = quantile(`MES-lineage`, 0.99, na.rm = TRUE),
  vascular_mean = mean(vascular_niche, na.rm = TRUE),
  vascular_q90 = quantile(vascular_niche, 0.90, na.rm = TRUE),
  vascular_q95 = quantile(vascular_niche, 0.95, na.rm = TRUE),
  vascular_q99 = quantile(vascular_niche, 0.99, na.rm = TRUE)
), by = slice]

sel_pool <- merge(sel_pool, vis_q, by = "slice", all.x = TRUE)
mes_q95_cut <- median(sel_pool$mes_lineage_q95, na.rm = TRUE)
vascular_q95_cut <- median(sel_pool$vascular_q95, na.rm = TRUE)
mes_q99_cut <- quantile(sel_pool$mes_lineage_q99, 0.75, na.rm = TRUE)
vascular_q99_cut <- quantile(sel_pool$vascular_q99, 0.75, na.rm = TRUE)
sel_pool[, abs_from_median := abs(observed_rho - med_rho)]
sel_pool[, enough_spots := n_spots >= 1000]
sel_pool[, mes_visible := mes_lineage_q99 >= mes_q99_cut]
sel_pool[, vascular_visible := vascular_q99 >= vascular_q99_cut]
sel_pool[, positive_rho := observed_rho > 0]
selected <- sel_pool[enough_spots == TRUE & mes_visible == TRUE &
                       vascular_visible == TRUE & positive_rho == TRUE][
                         order(abs_from_median, -n_spots)][1]
if (nrow(selected) == 0) {
  selected <- sel_pool[enough_spots == TRUE & positive_rho == TRUE][
    order(abs_from_median, -mes_lineage_q95, -vascular_q95)][1]
}
selected_slice <- selected$slice

selection_table <- sel_pool[order(abs_from_median, -n_spots),
  .(slice, n_spots, observed_rho, abs_from_median, p_emp_two_sided,
    mes_lineage_mean, mes_lineage_q90, mes_lineage_q95, mes_lineage_q99,
    vascular_mean, vascular_q90, vascular_q95, vascular_q99,
    enough_spots, mes_visible, vascular_visible, positive_rho)]
fwrite(selection_table, file.path(tab_dir, "R9_C3_3_representative_slice_selection.csv"))

## ---- Forest plot helper -----------------------------------------------------
make_forest <- function(preds, file_stub, width = 10, height = 8) {
  d <- perslice[k == 6 & response %in% main_responses & predictor %in% preds,
                .(slice, response, predictor, observed_rho, p_emp_two_sided)]
  s <- summ[k == 6 & response %in% main_responses & predictor %in% preds,
            .(response, predictor,
              observed_rho = median_rho,
              xmin = median_rho_ci_lo,
              xmax = median_rho_ci_hi,
              n_positive, n_negative, sign_p_greater, wilcox_p_greater,
              fisher_z_p_greater, median_emp_p_two_sided,
              n_slice_emp_p_lt_0p05)]
  s[, slice := "Summary median"]
  s[, emp_label := sprintf("emp<0.05: %d/18\nmed emp p=%.2f",
                           n_slice_emp_p_lt_0p05, median_emp_p_two_sided)]
  d[, `:=`(xmin = NA_real_, xmax = NA_real_, n_positive = NA_integer_,
           n_negative = NA_integer_, sign_p_greater = NA_real_,
           wilcox_p_greater = NA_real_, fisher_z_p_greater = NA_real_,
           median_emp_p_two_sided = NA_real_, emp_label = NA_character_)]
  d[, vector_shift_emp_lt_0p05 := p_emp_two_sided < 0.05]
  s[, `:=`(p_emp_two_sided = NA_real_, vector_shift_emp_lt_0p05 = NA)]
  plot_df <- rbindlist(list(d, s), fill = TRUE)
  slice_levels <- c("Summary median", sort(unique(d$slice)))
  plot_df[, slice := factor(slice, levels = rev(slice_levels))]
  plot_df[, predictor := factor(predictor, levels = preds)]
  plot_df[, response := factor(response, levels = main_responses)]
  plot_df[, vector_shift_emp_lt_0p05 := factor(vector_shift_emp_lt_0p05,
                                               levels = c(FALSE, TRUE),
                                               labels = c("empirical p >= 0.05",
                                                          "empirical p < 0.05"))]

  p <- ggplot(plot_df, aes(x = observed_rho, y = slice)) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey45") +
    geom_point(data = plot_df[as.character(slice) != "Summary median"],
               aes(shape = vector_shift_emp_lt_0p05),
               size = 1.75, color = "grey25", alpha = 0.85) +
    geom_errorbarh(data = plot_df[as.character(slice) == "Summary median"],
                   aes(xmin = xmin, xmax = xmax), height = 0.18,
                   linewidth = 0.65, color = "#b2182b") +
    geom_point(data = plot_df[as.character(slice) == "Summary median"],
               size = 2.5, color = "#b2182b") +
    geom_text(data = plot_df[as.character(slice) == "Summary median"],
              aes(x = 0.44, label = emp_label), inherit.aes = FALSE,
              y = Inf, vjust = 1.35, hjust = 1, size = 2.1, color = "grey25") +
    facet_grid(response ~ predictor, scales = "free_x") +
    scale_shape_manual(values = c("empirical p >= 0.05" = 16,
                                  "empirical p < 0.05" = 17),
                       na.translate = FALSE, name = "Vector-shift null") +
    coord_cartesian(xlim = c(-0.45, 0.45), clip = "off") +
    labs(x = "Per-slice Spearman rho", y = NULL,
         title = "Slice-level spatial association with vector-shift null",
         subtitle = "Points are slices; triangles mark slice-wise vector-shift empirical p < 0.05; red interval is bootstrap median rho 95% CI") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", color = "grey75"),
          axis.text.y = element_text(size = 6),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(fig_dir, paste0(file_stub, ".pdf")), p, width = width, height = height, units = "in")
  ggsave(file.path(fig_dir, paste0(file_stub, ".png")), p, width = width, height = height, units = "in", dpi = 300)
  fwrite(plot_df, file.path(tab_dir, paste0(file_stub, "_source.csv")))
}

make_forest(main_predictors, "R9_C3_3_forest_k6_vascular_vs_neuron_control")
make_forest(supp_predictors, "R9_C3_3_forest_k6_vascular_myeloid_neuron_supplement", width = 12, height = 8)

## ---- Representative maps ----------------------------------------------------
map_dt <- scores[slice == selected_slice & k == 6]
features <- c("MES-lineage", "MES-V", "MES-I", "vascular_niche", "neuron_control_niche")
map_source <- map_dt[, c("spot_id", "slice", "image", "x", "y", features), with = FALSE]

nearest_pitch <- function(xy) {
  dm <- as.matrix(dist(xy))
  diag(dm) <- Inf
  median(apply(dm, 1, min), na.rm = TRUE)
}
pitch <- nearest_pitch(as.matrix(map_source[, .(x, y)]))
scale_um <- 500
bar_len <- pitch * (scale_um / 100)  # Visium center-to-center pitch is 100 um.
x0 <- min(map_source$x) + 0.06 * diff(range(map_source$x))
y0 <- max(map_source$y) - 0.06 * diff(range(map_source$y))

plot_map <- function(feature, show_title = TRUE) {
  d <- copy(map_source)
  d[, value := get(feature)]
  upper <- quantile(d$value, 0.98, na.rm = TRUE)
  if (!is.finite(upper) || upper <= 0) upper <- max(d$value, na.rm = TRUE)
  p <- ggplot(d, aes(x = x, y = y, color = value)) +
    geom_point(size = 0.72, alpha = 0.95) +
    annotate("segment", x = x0, xend = x0 + bar_len, y = y0, yend = y0,
             linewidth = 0.75, color = "black") +
    annotate("text", x = x0 + bar_len / 2, y = y0 - 0.035 * diff(range(d$y)),
             label = paste0(scale_um, " um"), size = 2.8) +
    scale_y_reverse() +
    coord_fixed() +
    scale_color_gradientn(colors = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
                          limits = c(0, upper), oob = scales::squish,
                          name = "weight") +
    labs(title = if (show_title) paste0(selected_slice, " | ", feature) else feature,
         subtitle = if (show_title) "Representative slice; Visium spot-deconvolution level" else NULL,
         x = NULL, y = NULL) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 7),
          legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "white", color = NA))
  safe_feature <- gsub("[^A-Za-z0-9]+", "_", feature)
  ggsave(file.path(fig_dir, paste0("R9_C3_3_map_", safe_feature, "_", gsub("[#]", "", selected_slice), ".pdf")),
         p, width = 5.4, height = 4.4, units = "in", bg = "white")
  ggsave(file.path(fig_dir, paste0("R9_C3_3_map_", safe_feature, "_", gsub("[#]", "", selected_slice), ".png")),
         p, width = 5.4, height = 4.4, units = "in", dpi = 300, bg = "white")
  list(plot = p, feature = feature, color_upper_q98 = upper)
}
map_info <- rbindlist(lapply(features, function(ff) {
  x <- plot_map(ff)
  data.table(feature = x$feature, color_upper_q98 = x$color_upper_q98)
}))

pair_left <- plot_map("MES-lineage", show_title = FALSE)$plot
pair_right <- plot_map("vascular_niche", show_title = FALSE)$plot
pair <- gridExtra::arrangeGrob(pair_left, pair_right, ncol = 2,
                               top = paste0(selected_slice, " | paired spatial maps: MES-lineage and vascular niche"))
ggsave(file.path(fig_dir, paste0("R9_C3_3_pairmap_MESlineage_vascular_", gsub("[#]", "", selected_slice), ".pdf")),
       pair, width = 10.8, height = 4.8, units = "in", bg = "white")
ggsave(file.path(fig_dir, paste0("R9_C3_3_pairmap_MESlineage_vascular_", gsub("[#]", "", selected_slice), ".png")),
       pair, width = 10.8, height = 4.8, units = "in", dpi = 300, bg = "white")

fwrite(map_source, file.path(tab_dir, "R9_C3_3_representative_map_source.csv"))
fwrite(data.table(selected_slice = selected_slice,
                  criterion = "k=6 MES-lineage~vascular_niche observed rho closest to 18-slice median among slices with n_spots >= 1000, positive rho, MES-lineage q99 >= cross-slice upper quartile q99, and vascular q99 >= cross-slice upper quartile q99; ties resolved by larger n_spots",
                  n_spots = selected$n_spots,
                  observed_rho = selected$observed_rho,
                  cross_slice_median_rho = med_rho,
                  abs_from_median = selected$abs_from_median,
                  p_emp_two_sided = selected$p_emp_two_sided,
                  mes_lineage_q95 = selected$mes_lineage_q95,
                  mes_lineage_q95_cut = mes_q95_cut,
                  mes_lineage_q99 = selected$mes_lineage_q99,
                  mes_lineage_q99_cut = as.numeric(mes_q99_cut),
                  vascular_q95 = selected$vascular_q95,
                  vascular_q95_cut = vascular_q95_cut,
                  vascular_q99 = selected$vascular_q99,
                  vascular_q99_cut = as.numeric(vascular_q99_cut),
                  pitch_px = pitch,
                  scale_bar_um = scale_um,
                  scale_bar_px = bar_len),
       file.path(tab_dir, "R9_C3_3_figure_decision_source.csv"))
fwrite(map_info, file.path(tab_dir, "R9_C3_3_map_color_scale_source.csv"))

cat("\n== C3-3 representative slice:\n")
print(fread(file.path(tab_dir, "R9_C3_3_figure_decision_source.csv")))
cat("\n== C3-3 figures written to:", fig_dir, "\n")
cat("[STOP C3-3] Review forest plots and representative maps. Figures are candidates; final manuscript role must be decided before copying to 图片表格.\n")
