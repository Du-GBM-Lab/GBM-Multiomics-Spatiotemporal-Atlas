#!/usr/bin/env Rscript
# =============================================================================
# R9 Batch 2 supplement | malignant-gradient border-core composition
#
# Nature:
#   Static composition along the RCTD malignant-total axis. This is not an
#   invasion-front, infiltration, displacement, recruitment, or causal analysis.
#   It only describes component changes along a malignant-density gradient in
#   slices that visibly contain both low- and high-malignant regions.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
  library(ggplot2)
})

base_dir <- getwd()
out_dir <- file.path(base_dir, "tables/R9_batch2_landscape_gradient/border_core_malignant_gradient")
fig_dir <- file.path(base_dir, "figures/R9_batch2_landscape_gradient/border_core_malignant_gradient")
doc_dir <- file.path(base_dir, "docs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(doc_dir, showWarnings = FALSE, recursive = TRUE)

weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")

bin_breaks <- c(-Inf, 0.05, 0.15, 0.30, 0.50, Inf)
bin_labels <- c("B1_low_malignant", "B2_low_intermediate", "B3_intermediate", "B4_high", "B5_core_like")
low_threshold <- 0.05
high_threshold <- 0.50
min_low_spots <- 10L
min_high_spots <- 10L
min_per_bin_spots <- 1L

w <- as.data.table(qs2::qs_read(weights_path))
metadata_cols <- c("spot_id", "slice", "image", "x", "y")
required <- c(
  metadata_cols,
  "Subtype1", "Subtype2", "Subtype3", "Subtype4",
  "Astrocytes", "Oligodendrocytes", "OPCs",
  "Macrophages", "Microglial", "Monocytes",
  "Endothelial", "Mural cells", "Neurons",
  "T cells", "B cells", "NK cells", "cDCs", "pDCs", "Ependymal cells", "Radial glial"
)
missing <- setdiff(required, names(w))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

w[, `:=`(
  `NPC-P` = Subtype1,
  `OPC-M` = Subtype2,
  `MES-V` = Subtype3,
  `MES-I` = Subtype4,
  malignant_total = Subtype1 + Subtype2 + Subtype3 + Subtype4,
  Astrocyte = Astrocytes,
  Oligo_lineage = Oligodendrocytes + OPCs,
  TAM = Macrophages + Microglial + Monocytes,
  Vascular = Endothelial + `Mural cells`,
  Lymphoid = `T cells` + `B cells` + `NK cells`,
  DC = cDCs + pDCs,
  Other_nonmalignant = Ambiguous + `Ependymal cells` + `Radial glial`
)]
w[, nonmalignant_total := pmax(0, 1 - malignant_total)]
w[, malignant_gradient_bin := cut(malignant_total, breaks = bin_breaks, labels = bin_labels, include.lowest = TRUE, right = TRUE)]

bin_counts <- w[, .N, by = .(slice, malignant_gradient_bin)]
bin_counts_wide <- dcast(bin_counts, slice ~ malignant_gradient_bin, value.var = "N", fill = 0)
for (bn in bin_labels) if (!bn %in% names(bin_counts_wide)) bin_counts_wide[, (bn) := 0L]

slice_qc <- w[, .(
  n_spots = .N,
  malignant_min = min(malignant_total, na.rm = TRUE),
  malignant_q05 = quantile(malignant_total, 0.05, na.rm = TRUE),
  malignant_median = median(malignant_total, na.rm = TRUE),
  malignant_q95 = quantile(malignant_total, 0.95, na.rm = TRUE),
  malignant_max = max(malignant_total, na.rm = TRUE),
  low_malignant_spots = sum(malignant_total <= low_threshold, na.rm = TRUE),
  high_malignant_spots = sum(malignant_total >= high_threshold, na.rm = TRUE)
), by = slice]
slice_qc <- merge(slice_qc, bin_counts_wide, by = "slice", all.x = TRUE)
slice_qc[, all_bins_present := Reduce(`&`, lapply(.SD, function(x) x >= min_per_bin_spots)), .SDcols = bin_labels]
slice_qc[, border_gradient_captured := low_malignant_spots >= min_low_spots & high_malignant_spots >= min_high_spots & all_bins_present]
slice_qc[, exclusion_reason := fifelse(
  border_gradient_captured, "included",
  paste(
    fifelse(low_malignant_spots < min_low_spots, "low_malignant_region_insufficient", ""),
    fifelse(high_malignant_spots < min_high_spots, "high_malignant_region_insufficient", ""),
    fifelse(!all_bins_present, "missing_one_or_more_bins", ""),
    sep = ";"
  )
)]

included_slices <- slice_qc[border_gradient_captured == TRUE, slice]

group_cols <- c(
  "malignant_total", "nonmalignant_total",
  "NPC-P", "OPC-M", "MES-V", "MES-I",
  "Astrocyte", "Oligo_lineage", "TAM", "Vascular", "Neurons", "Lymphoid", "DC", "Other_nonmalignant"
)

long <- melt(
  w[slice %in% included_slices],
  id.vars = c("spot_id", "slice", "malignant_gradient_bin", "malignant_total"),
  measure.vars = group_cols,
  variable.name = "component",
  value.name = "weight"
)

per_slice_bin <- long[, .(
  n_spots = .N,
  mean_weight = mean(weight, na.rm = TRUE),
  median_axis_malignant_total = median(malignant_total, na.rm = TRUE)
), by = .(slice, malignant_gradient_bin, component)]

bin_summary <- per_slice_bin[, .(
  n_slices = .N,
  mean_across_slices = mean(mean_weight, na.rm = TRUE),
  sd_across_slices = sd(mean_weight, na.rm = TRUE),
  min_across_slices = min(mean_weight, na.rm = TRUE),
  max_across_slices = max(mean_weight, na.rm = TRUE),
  median_across_slices = median(mean_weight, na.rm = TRUE),
  q25_across_slices = quantile(mean_weight, 0.25, na.rm = TRUE),
  q75_across_slices = quantile(mean_weight, 0.75, na.rm = TRUE),
  cv_sd_over_abs_mean = ifelse(abs(mean(mean_weight, na.rm = TRUE)) > 0, sd(mean_weight, na.rm = TRUE) / abs(mean(mean_weight, na.rm = TRUE)), NA_real_)
), by = .(component, malignant_gradient_bin)]
setorder(bin_summary, component, malignant_gradient_bin)

wide <- dcast(per_slice_bin, slice + component ~ malignant_gradient_bin, value.var = "mean_weight")
for (bn in bin_labels) if (!bn %in% names(wide)) wide[, (bn) := NA_real_]
trend_by_slice <- wide[, {
  vals <- as.numeric(.SD)
  diffs <- diff(vals)
  data.table(
    low_bin = vals[1],
    core_bin = vals[5],
    core_minus_low = vals[5] - vals[1],
    low_minus_core = vals[1] - vals[5],
    adjacent_increase_count = sum(diffs > 0, na.rm = TRUE),
    adjacent_decrease_count = sum(diffs < 0, na.rm = TRUE),
    monotonic_increasing = all(diffs >= 0, na.rm = TRUE),
    monotonic_decreasing = all(diffs <= 0, na.rm = TRUE),
    core_gt_low = vals[5] > vals[1],
    low_gt_core = vals[1] > vals[5]
  )
}, by = .(slice, component), .SDcols = bin_labels]

trend_summary <- trend_by_slice[, .(
  n_slices = .N,
  median_core_minus_low = median(core_minus_low, na.rm = TRUE),
  median_low_minus_core = median(low_minus_core, na.rm = TRUE),
  core_gt_low_slices = sum(core_gt_low, na.rm = TRUE),
  low_gt_core_slices = sum(low_gt_core, na.rm = TRUE),
  monotonic_increasing_slices = sum(monotonic_increasing, na.rm = TRUE),
  monotonic_decreasing_slices = sum(monotonic_decreasing, na.rm = TRUE),
  median_adjacent_increase_count = median(adjacent_increase_count, na.rm = TRUE),
  median_adjacent_decrease_count = median(adjacent_decrease_count, na.rm = TRUE)
), by = component]

adjacent <- per_slice_bin[order(slice, component, malignant_gradient_bin)]
adjacent <- adjacent[, .(
  adjacent_pair = paste0(bin_labels[-length(bin_labels)], "_to_", bin_labels[-1]),
  increase = diff(mean_weight) > 0,
  decrease = diff(mean_weight) < 0,
  delta_next_minus_current = diff(mean_weight)
), by = .(slice, component)]
adjacent_summary <- adjacent[, .(
  n_slices = .N,
  increase_slices = sum(increase, na.rm = TRUE),
  decrease_slices = sum(decrease, na.rm = TRUE),
  median_delta_next_minus_current = median(delta_next_minus_current, na.rm = TRUE)
), by = .(component, adjacent_pair)]

key_components <- c("malignant_total", "nonmalignant_total", "NPC-P", "OPC-M", "MES-V", "MES-I", "Astrocyte", "Oligo_lineage", "TAM")
key_flags <- merge(
  trend_summary[component %in% key_components],
  bin_summary[component %in% key_components, .(median_cv = median(cv_sd_over_abs_mean, na.rm = TRUE)), by = component],
  by = "component",
  all.x = TRUE
)
key_flags[, clean_increasing := monotonic_increasing_slices >= ceiling(0.70 * n_slices) & core_gt_low_slices >= ceiling(0.80 * n_slices) & median_cv < 1.5]
key_flags[, clean_decreasing := monotonic_decreasing_slices >= ceiling(0.70 * n_slices) & low_gt_core_slices >= ceiling(0.80 * n_slices) & median_cv < 1.5]

malignant_subtype_clean_n <- key_flags[component %in% c("NPC-P", "OPC-M", "MES-V", "MES-I"), sum(clean_increasing, na.rm = TRUE)]
overview_candidate <- length(included_slices) >= 8 &&
  key_flags[component == "malignant_total", clean_increasing] == TRUE &&
  key_flags[component == "nonmalignant_total", clean_decreasing] == TRUE &&
  malignant_subtype_clean_n >= 3
placement <- if (overview_candidate) {
  "supplement_landscape_candidate_with_possible_main_overview_if_guidance_accepts_static_composition"
} else if (length(included_slices) >= 8 &&
  key_flags[component == "malignant_total", clean_increasing] == TRUE &&
  key_flags[component == "nonmalignant_total", clean_decreasing] == TRUE) {
  "supplement_landscape_only_not_main_overview_due_to_subtype_noise"
} else {
  "discard_internal_audit_due_to_unclean_or_underpowered_malignant_gradient"
}

params <- data.table(
  parameter = c(
    "input_weights", "malignant_total_definition", "bin_breaks", "bin_labels",
    "included_slice_rule", "all_spots_retained", "statistics", "placement_rule"
  ),
  value = c(
    weights_path,
    "Subtype1 + Subtype2 + Subtype3 + Subtype4 RCTD weights",
    paste(bin_breaks, collapse = ", "),
    paste(bin_labels, collapse = ", "),
    paste0("low_malignant_spots >= ", min_low_spots, " at malignant_total <= ", low_threshold,
           "; high_malignant_spots >= ", min_high_spots, " at malignant_total >= ", high_threshold,
           "; all five bins present"),
    "TRUE; binning is descriptive axis assignment, not spot exclusion within included slices",
    "per-slice composition then cross-slice summary over included slices only; no relation-level p-value",
    "overview candidate requires N>=8, malignant_total clean increasing, nonmalignant_total clean decreasing, and >=3/4 malignant subtypes clean increasing"
  )
)

fwrite(slice_qc, file.path(out_dir, "border_core_step0_slice_gradient_qc.csv"))
fwrite(per_slice_bin, file.path(out_dir, "border_core_composition_perslice_bins.csv"))
fwrite(bin_summary, file.path(out_dir, "border_core_composition_bin_summary.csv"))
fwrite(trend_by_slice, file.path(out_dir, "border_core_trend_by_slice.csv"))
fwrite(trend_summary, file.path(out_dir, "border_core_trend_summary.csv"))
fwrite(adjacent_summary, file.path(out_dir, "border_core_adjacent_bin_consistency.csv"))
fwrite(key_flags, file.path(out_dir, "border_core_key_component_cleanliness_flags.csv"))
fwrite(params, file.path(out_dir, "border_core_parameters_source.csv"))
qs2::qs_save(
  list(
    slice_qc = slice_qc,
    per_slice_bin = per_slice_bin,
    bin_summary = bin_summary,
    trend_by_slice = trend_by_slice,
    trend_summary = trend_summary,
    adjacent_summary = adjacent_summary,
    key_flags = key_flags,
    parameters = params
  ),
  file.path(out_dir, "border_core_malignant_gradient_results.qs2")
)

plot_dt <- bin_summary[component %in% key_components]
p <- ggplot(plot_dt, aes(x = malignant_gradient_bin, y = median_across_slices, color = component, group = component)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  geom_ribbon(aes(ymin = q25_across_slices, ymax = q75_across_slices, fill = component), alpha = 0.07, color = NA) +
  labs(x = "Malignant-density bin", y = "Median RCTD weight across included slices", color = NULL, fill = NULL) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
ggsave(file.path(fig_dir, "malignant_gradient_composition_curve.pdf"), p, width = 7.2, height = 4.4)
ggsave(file.path(fig_dir, "malignant_gradient_composition_curve.png"), p, width = 7.2, height = 4.4, dpi = 300)

fmt_component <- function(comp) {
  vals <- bin_summary[component == comp][order(malignant_gradient_bin), median_across_slices]
  paste(signif(vals, 4), collapse = " -> ")
}
flag_line <- function(comp) {
  k <- key_flags[component == comp]
  if (!nrow(k)) return("NA")
  paste0(
    "core>low ", k$core_gt_low_slices, "/", k$n_slices,
    "; low>core ", k$low_gt_core_slices, "/", k$n_slices,
    "; monotonic inc ", k$monotonic_increasing_slices, "/", k$n_slices,
    "; monotonic dec ", k$monotonic_decreasing_slices, "/", k$n_slices,
    "; median CV ", signif(k$median_cv, 4)
  )
}

stop_lines <- c(
  "# R9 Batch 2 supplement border-core / malignant-gradient STOP",
  "",
  "性质声明：本分析只描述沿 RCTD malignant-total 轴的静态成分变化。它不是 invasion-front、浸润推进、驱赶、招募或因果分析；静态 Visium 切片只能支持 malignant-gradient composition。",
  "",
  "## Step 0: captured-gradient slice gate",
  paste0("- bin rule：", paste(bin_labels, c("<=0.05", "0.05-0.15", "0.15-0.30", "0.30-0.50", ">0.50"), sep = "=", collapse = "; "), "。"),
  paste0("- included slices：", length(included_slices), "/18；", paste(included_slices, collapse = ", "), "。"),
  paste0("- excluded slices：", 18 - length(included_slices), "/18；见 `border_core_step0_slice_gradient_qc.csv`。"),
  "",
  "## Step 1-2: composition along malignant-density axis",
  paste0("- malignant_total median trajectory：", fmt_component("malignant_total"), "；", flag_line("malignant_total"), "。"),
  paste0("- nonmalignant_total median trajectory：", fmt_component("nonmalignant_total"), "；", flag_line("nonmalignant_total"), "。"),
  paste0("- NPC-P median trajectory：", fmt_component("NPC-P"), "；", flag_line("NPC-P"), "。"),
  paste0("- OPC-M median trajectory：", fmt_component("OPC-M"), "；", flag_line("OPC-M"), "。"),
  paste0("- MES-V median trajectory：", fmt_component("MES-V"), "；", flag_line("MES-V"), "。"),
  paste0("- MES-I median trajectory：", fmt_component("MES-I"), "；", flag_line("MES-I"), "。"),
  paste0("- TAM median trajectory：", fmt_component("TAM"), "；", flag_line("TAM"), "。TAM 只作 non-malignant background composition，不写 recruitment。"),
  paste0("- clean increasing malignant subtypes：", malignant_subtype_clean_n, "/4。"),
  "",
  "## Placement",
  paste0("- self-assessed placement：", placement, "。"),
  "- 曲线本身只支持 `along the malignant-density axis, component weights changed` 这类描述；若要下任何 region relation 主张，必须另回 v3-B。",
  "",
  "## Forbidden wording",
  "- 禁写 infiltration front 作动态、invasion、frontier advance、tumor pushing/driving normal cells、recruitment、drives、malignant cells displace/expel、causality、interaction。",
  "",
  "## Source Data",
  paste0("- input weights：`", weights_path, "`"),
  paste0("- output dir：`", out_dir, "`")
)
writeLines(stop_lines, file.path(doc_dir, "R9_batch2_border_core_malignant_gradient_STOP.md"), useBytes = TRUE)

message("Wrote border-core malignant-gradient outputs to: ", out_dir)
message("Wrote STOP to: ", file.path(doc_dir, "R9_batch2_border_core_malignant_gradient_STOP.md"))
