#!/usr/bin/env Rscript
# =============================================================================
# R9 Batch 2 project 3 | vascular-distance x Buffa hypoxia dual-axis landscape
#
# Nature:
#   Descriptive two-axis macro-landscape. The axes are independent sources:
#   distance to dominant endothelial/mural RCTD spots and Buffa hypoxia score.
#   TAM/myeloid is shown only as neutral background composition and is not used
#   to make a TAM hypoxia-preference, recruitment, interaction, or causality
#   claim.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
  library(ggplot2)
})

base_dir <- getwd()
out_dir <- file.path(base_dir, "tables/R9_batch2_landscape_gradient/dual_axis_vascular_hypoxia")
fig_dir <- file.path(base_dir, "figures/R9_batch2_landscape_gradient/dual_axis_vascular_hypoxia")
doc_dir <- file.path(base_dir, "docs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(doc_dir, showWarnings = FALSE, recursive = TRUE)

weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
score_path <- file.path(base_dir, "tables/R9_stage3_signature_scores/R9_stage3_IVY_hypoxia_scores_per_spot.csv")
distance_path <- file.path(base_dir, "tables/R9_batch2_landscape_gradient/vascular_distance_gradient/vascular_distance_per_spot.csv")
territory_summary_path <- file.path(base_dir, "tables/R9_batch2_landscape_gradient/territory_segregation/territory_segregation_global_summary.csv")
vascular_gradient_summary_path <- file.path(base_dir, "tables/R9_batch2_landscape_gradient/vascular_distance_gradient/vascular_distance_group_slopes_summary.csv")

rank_bins <- function(x, n_bins = 5L) {
  if (all(is.na(x))) return(rep(NA_integer_, length(x)))
  r <- frank(x, ties.method = "average", na.last = "keep")
  out <- floor(((r - 0.5) / sum(!is.na(x))) * n_bins) + 1L
  pmin(n_bins, pmax(1L, as.integer(out)))
}

w <- as.data.table(qs2::qs_read(weights_path))
scores <- fread(score_path)
distance <- fread(distance_path)

required_scores <- c("spot_id", "Hypoxia_Buffa")
missing_scores <- setdiff(required_scores, names(scores))
if (length(missing_scores)) stop("Missing score columns: ", paste(missing_scores, collapse = ", "))

required_weights <- c("spot_id", "Subtype1", "Subtype2", "Subtype3", "Subtype4", "Macrophages", "Microglial", "Monocytes", "Endothelial", "Mural cells", "Neurons")
missing_weights <- setdiff(required_weights, names(w))
if (length(missing_weights)) stop("Missing RCTD columns: ", paste(missing_weights, collapse = ", "))

w[, `:=`(
  `NPC-P` = Subtype1,
  `OPC-M` = Subtype2,
  `MES-V` = Subtype3,
  `MES-I` = Subtype4,
  TAM = Macrophages + Microglial + Monocytes
)]

dt <- merge(
  w[, .(spot_id, slice, image, x, y, `NPC-P`, `OPC-M`, `MES-V`, `MES-I`, TAM, Endothelial, `Mural cells`, Neurons)],
  distance[, .(spot_id, distance_to_nearest_vascular, vascular_distance_bin, vascular_distance_available)],
  by = "spot_id",
  all.x = TRUE
)
dt <- merge(dt, scores[, .(spot_id, Hypoxia_Buffa)], by = "spot_id", all.x = TRUE)
dt[, hypoxia_buffa_bin := rank_bins(Hypoxia_Buffa, 5L), by = slice]

group_cols <- c("NPC-P", "OPC-M", "MES-V", "MES-I", "TAM", "Endothelial", "Mural cells", "Neurons")
long <- melt(
  dt[vascular_distance_available == TRUE & !is.na(hypoxia_buffa_bin)],
  id.vars = c("spot_id", "slice", "vascular_distance_bin", "hypoxia_buffa_bin", "distance_to_nearest_vascular", "Hypoxia_Buffa"),
  measure.vars = group_cols,
  variable.name = "group",
  value.name = "weight"
)

grid_perslice <- long[, .(
  n_spots = .N,
  mean_weight = mean(weight, na.rm = TRUE),
  median_distance = median(distance_to_nearest_vascular, na.rm = TRUE),
  median_hypoxia = median(Hypoxia_Buffa, na.rm = TRUE)
), by = .(slice, vascular_distance_bin, hypoxia_buffa_bin, group)]

grid_summary <- grid_perslice[, .(
  n_slices = .N,
  median_mean_weight = median(mean_weight, na.rm = TRUE),
  q25_mean_weight = quantile(mean_weight, 0.25, na.rm = TRUE),
  q75_mean_weight = quantile(mean_weight, 0.75, na.rm = TRUE)
), by = .(vascular_distance_bin, hypoxia_buffa_bin, group)]

centers <- long[, .(
  weight_sum = sum(weight, na.rm = TRUE),
  vascular_bin_center = ifelse(sum(weight, na.rm = TRUE) > 0, weighted.mean(vascular_distance_bin, weight, na.rm = TRUE), NA_real_),
  hypoxia_bin_center = ifelse(sum(weight, na.rm = TRUE) > 0, weighted.mean(hypoxia_buffa_bin, weight, na.rm = TRUE), NA_real_)
), by = .(slice, group)]

center_summary <- centers[, .(
  n_slices = .N,
  median_vascular_bin_center = median(vascular_bin_center, na.rm = TRUE),
  q25_vascular_bin_center = quantile(vascular_bin_center, 0.25, na.rm = TRUE),
  q75_vascular_bin_center = quantile(vascular_bin_center, 0.75, na.rm = TRUE),
  median_hypoxia_bin_center = median(hypoxia_bin_center, na.rm = TRUE),
  q25_hypoxia_bin_center = quantile(hypoxia_bin_center, 0.25, na.rm = TRUE),
  q75_hypoxia_bin_center = quantile(hypoxia_bin_center, 0.75, na.rm = TRUE)
), by = group]

params <- data.table(
  parameter = c("input_weights", "input_hypoxia_scores", "input_vascular_distance", "vascular_axis", "hypoxia_axis", "binning", "TAM_handling", "statistics"),
  value = c(
    weights_path,
    score_path,
    distance_path,
    "distance to nearest dominant RCTD Endothelial/Mural spot from script 52",
    "Buffa hypoxia score from locked R9_stage3_IVY_hypoxia_scores_per_spot.csv",
    "within-slice rank bins, 1=near/low and 5=far/high",
    "neutral background composition only; no TAM hypoxia relation test or claim",
    "descriptive per-slice two-axis composition only"
  )
)

fwrite(dt[, .(spot_id, slice, image, x, y, distance_to_nearest_vascular, vascular_distance_bin, vascular_distance_available, Hypoxia_Buffa, hypoxia_buffa_bin)], file.path(out_dir, "dual_axis_per_spot_axes.csv"))
fwrite(grid_perslice, file.path(out_dir, "dual_axis_composition_perslice_grid.csv"))
fwrite(grid_summary, file.path(out_dir, "dual_axis_composition_summary_grid.csv"))
fwrite(centers, file.path(out_dir, "dual_axis_group_centers_perslice.csv"))
fwrite(center_summary, file.path(out_dir, "dual_axis_group_centers_summary.csv"))
fwrite(params, file.path(out_dir, "dual_axis_parameters_source.csv"))

qs2::qs_save(
  list(
    axes = dt,
    grid_perslice = grid_perslice,
    grid_summary = grid_summary,
    centers = centers,
    center_summary = center_summary,
    parameters = params
  ),
  file.path(out_dir, "dual_axis_vascular_hypoxia_results.qs2")
)

for (g in c("NPC-P", "OPC-M", "MES-V", "MES-I", "TAM")) {
  pd <- grid_summary[group == g]
  p <- ggplot(pd, aes(x = vascular_distance_bin, y = hypoxia_buffa_bin, fill = median_mean_weight)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_x_continuous(breaks = 1:5, labels = c("near", "2", "3", "4", "far")) +
    scale_y_continuous(breaks = 1:5, labels = c("low", "2", "3", "4", "high")) +
    scale_fill_viridis_c(option = "C") +
    labs(x = "Distance-to-vascular bin", y = "Buffa hypoxia bin", fill = "Median weight", title = g) +
    theme_classic(base_size = 9)
  safe_g <- gsub("[^A-Za-z0-9]+", "_", g)
  ggsave(file.path(fig_dir, paste0("双轴分布_", safe_g, ".pdf")), p, width = 3.8, height = 3.5)
  ggsave(file.path(fig_dir, paste0("双轴分布_", safe_g, ".png")), p, width = 3.8, height = 3.5, dpi = 300)
}

read_metric <- function(path, metric_name) {
  if (!file.exists(path)) return(NA_character_)
  x <- fread(path)
  val <- x[metric == metric_name, value]
  if (length(val)) as.character(val[1]) else NA_character_
}

territory_fdr <- read_metric(territory_summary_path, "subtype_same_neighbor_FDR_pass")
territory_pos <- read_metric(territory_summary_path, "subtype_same_neighbor_positive_vs_null_median")
neuron_fdr <- read_metric(territory_summary_path, "neuron_control_FDR_pass")

vg <- if (file.exists(vascular_gradient_summary_path)) fread(vascular_gradient_summary_path) else data.table()
mesv_line <- if (nrow(vg[group == "MES-V"])) vg[group == "MES-V"][1] else NULL
tam_line <- if (nrow(vg[group == "TAM"])) vg[group == "TAM"][1] else NULL
center_mesv <- center_summary[group == "MES-V"][1]
center_tam <- center_summary[group == "TAM"][1]

stop_lines <- c(
  "# R9 Batch 2 landscape / gradient STOP",
  "",
  "性质声明：Batch 2 三项均为 Visium spot-level 空间景观/梯度描述；不证明因果、演化、物理接触、single-cell colocalization、interaction 或 recruitment。relation-level 主张仍回到已锁 v3-B。",
  "",
  "## Project 1: malignant subtype territory / spatial segregation (51_)",
  paste0("- 输出目录：`", file.path(base_dir, "tables/R9_batch2_landscape_gradient/territory_segregation"), "`"),
  paste0("- 方法：k=6 same-neighbor / assortativity / local entropy；toroidal coordinate-shift label-field null，未使用 CSR 或独立 label permutation。"),
  paste0("- subtype same-neighbor FDR pass：", territory_fdr, "/18；positive vs null median：", territory_pos, "/18。"),
  paste0("- neuron control segregated FDR pass：", neuron_fdr, "/18；此处 neuron 作为 null 行为裁判，不支持任何 TAM/myeloid 叙事。"),
  "- 自评落位：supplement landscape。若指导判断 MES 斑块化跨片足够稳定，正文最多一句作为空间格局背景；禁止写演化方向或亚型转化。",
  "",
  "## Project 2: distance-to-vascular gradient composition (52_)",
  paste0("- 输出目录：`", file.path(base_dir, "tables/R9_batch2_landscape_gradient/vascular_distance_gradient"), "`"),
  "- 血管定义：dominant RCTD label 为 Endothelial 或 Mural cells；独立于 MES marker。每个 spot 到同片最近 vascular spot 的距离进入 5 个 within-slice rank bins，未按 score 筛 spot。",
  paste0("- MES-V median near-minus-far：", if (!is.null(mesv_line)) signif(mesv_line$median_near_minus_far, 4) else "NA",
         "；near>far slices：", if (!is.null(mesv_line)) mesv_line$near_greater_than_far_slices else "NA", "/",
         if (!is.null(mesv_line)) mesv_line$n_slices else "NA", "。"),
  paste0("- TAM median near-minus-far：", if (!is.null(tam_line)) signif(tam_line$median_near_minus_far, 4) else "NA",
         "；near>far slices：", if (!is.null(tam_line)) tam_line$near_greater_than_far_slices else "NA", "/",
         if (!is.null(tam_line)) tam_line$n_slices else "NA", "。TAM 只作为与已锁 spatial null 对齐的背景检查，不写 recruitment。"),
  "- 自评落位：main overview candidate as descriptive gradient if MES-V near-vascular preference is visually stable; inferential relation wording still relies on v3-B, not on this curve alone.",
  "",
  "## Project 3: vascular-distance x Buffa hypoxia dual-axis landscape (53_)",
  paste0("- 输出目录：`", out_dir, "`"),
  "- 方法：距血管 rank bin x Buffa hypoxia rank bin 的 per-slice 2D composition；两轴均为独立来源。",
  paste0("- MES-V median center：vascular bin ", if (nrow(center_mesv)) signif(center_mesv$median_vascular_bin_center, 4) else "NA",
         "，hypoxia bin ", if (nrow(center_mesv)) signif(center_mesv$median_hypoxia_bin_center, 4) else "NA", "。"),
  paste0("- TAM median center：vascular bin ", if (nrow(center_tam)) signif(center_tam$median_vascular_bin_center, 4) else "NA",
         "，hypoxia bin ", if (nrow(center_tam)) signif(center_tam$median_hypoxia_bin_center, 4) else "NA",
         "；仅背景呈现，不做 TAM hypoxia relation test。"),
  "- 自评落位：supplement landscape / visualization for locked vascular and hypoxia context; no TAM-specific or hypoxia-driving claim.",
  "",
  "## 禁止写法",
  "- 禁写 causal proof、physical contact、single-cell colocalization、interaction、recruitment、trajectory/evolution、vascular drives MES、hypoxia/PAN drives MES、TAM recruitment。",
  "",
  "## Source data",
  paste0("- RCTD weights：`", weights_path, "`"),
  paste0("- Buffa hypoxia score：`", score_path, "`"),
  paste0("- Project 1 source：`", file.path(base_dir, "tables/R9_batch2_landscape_gradient/territory_segregation"), "`"),
  paste0("- Project 2 source：`", file.path(base_dir, "tables/R9_batch2_landscape_gradient/vascular_distance_gradient"), "`"),
  paste0("- Project 3 source：`", out_dir, "`")
)
writeLines(stop_lines, file.path(doc_dir, "R9_batch2_landscape_gradient_STOP.md"), useBytes = TRUE)

message("Wrote project 3 dual-axis outputs to: ", out_dir)
message("Wrote combined STOP to: ", file.path(doc_dir, "R9_batch2_landscape_gradient_STOP.md"))
