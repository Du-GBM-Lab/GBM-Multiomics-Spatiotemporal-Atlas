#!/usr/bin/env Rscript
# =============================================================================
# 07_发育时间_TF_ATAC验证 / 01_mouse_developmental
# Phase 1 - Script 03: 小鼠四亚型 跨阶段占比 可视化
# -----------------------------------------------------------------------------
# 输入: Script 02 产出的注释对象 + source_data 比例表 + sanity marker 表。
# 输出: 每张图独立 panel PDF + source data。
#
# 不使用极坐标环形图：弧长会绑定到每阶段总恶性细胞数，而该数值受捕获量影响，
# 不是生物学量。主图用可审计的堆叠比例柱 + 趋势线。
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

module_dir <- normalizePath("07_发育时间_TF_ATAC验证", winslash = "/", mustWork = TRUE)
project_dir <- file.path(module_dir, "01_mouse_developmental")
tables_dir <- file.path(module_dir, "tables")
outputs_dir <- file.path(project_dir, "outputs")
figures_dir <- file.path(project_dir, "figures")
src_data_dir <- file.path(figures_dir, "source_data")
logs_dir <- file.path(project_dir, "logs")

obj_path <- file.path(outputs_dir, "GSE278511小鼠恶性模型_四亚型投射.qs2")
prop_all_path <- file.path(src_data_dir, "mouse_phase1_stage_subtype_proportions.csv")
prop_hc_path <- file.path(src_data_dir, "mouse_phase1_stage_subtype_proportions_highconf.csv")
per_cell_path <- file.path(src_data_dir, "mouse_phase1_per_cell_assigned_subtype.csv")
sanity_path <- file.path(tables_dir, "mouse_phase1_sanity_markers.csv")

for (p in c(obj_path, prop_all_path, per_cell_path)) {
  stopifnot(file.exists(p))
}

subtype_levels <- c("NPC-P", "OPC-M", "MES-V", "MES-I")
stage_levels <- c("Preneoplastic", "Early lesion", "Mid lesion", "Endpoint")
subtype_colors <- c(
  "NPC-P" = "#0072B5",
  "OPC-M" = "#E18727",
  "MES-V" = "#20854E",
  "MES-I" = "#BC3C29"
)

base_theme <- theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(color = "black"),
    legend.key.size = unit(10, "pt"),
    plot.title = element_text(size = 10, face = "bold")
  )

save_pdf <- function(plot, name, w, h) {
  ggsave(
    file.path(figures_dir, name),
    plot = plot,
    width = w,
    height = h,
    useDingbats = FALSE
  )
}

prop_all <- read_csv(prop_all_path, show_col_types = FALSE) %>%
  mutate(
    stage_inferred = factor(stage_inferred, levels = stage_levels),
    assigned_subtype = factor(assigned_subtype, levels = subtype_levels)
  )
stage_n <- prop_all %>%
  distinct(stage_inferred, stage_total) %>%
  mutate(label = paste0("n=", scales::comma(stage_total)))

per_cell <- read_csv(per_cell_path, show_col_types = FALSE)

mouse_obj <- qs2::qs_read(obj_path)
red_name <- intersect(c("umap", "UMAP", "Umap"), names(mouse_obj@reductions))[1]
if (is.na(red_name)) {
  stop("No UMAP reduction found in object.")
}

pA <- DimPlot(
  mouse_obj,
  reduction = red_name,
  group.by = "assigned_subtype",
  pt.size = 0.25,
  cols = subtype_colors
) +
  base_theme +
  labs(title = "Cross-species assigned subtype") +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  )
save_pdf(pA, "A_小鼠四亚型投射UMAP.pdf", 5.4, 4.4)

pB <- ggplot(prop_all, aes(stage_inferred, proportion, fill = assigned_subtype)) +
  geom_col(width = 0.74, color = "white", linewidth = 0.25) +
  geom_text(
    data = stage_n,
    aes(x = stage_inferred, y = 1.03, label = label),
    inherit.aes = FALSE,
    size = 2.9,
    fontface = "italic"
  ) +
  scale_fill_manual(values = subtype_colors, drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.08), expand = c(0, 0)) +
  labs(x = NULL, y = "Fraction of malignant cells", fill = NULL, title = "Subtype composition across tumorigenesis stages") +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_pdf(pB, "B_小鼠四阶段亚型堆叠比例.pdf", 5.0, 4.2)

pC <- ggplot(prop_all, aes(stage_inferred, n_cells, fill = assigned_subtype)) +
  geom_col(width = 0.74, position = position_dodge(width = 0.8), color = "white", linewidth = 0.2) +
  scale_fill_manual(values = subtype_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::comma, expand = expansion(c(0, 0.05))) +
  labs(x = NULL, y = "Malignant cells (n)", fill = NULL, title = "Absolute subtype abundance per stage") +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_pdf(pC, "C_小鼠四阶段亚型绝对数量.pdf", 5.6, 4.2)

pD <- ggplot(
  prop_all,
  aes(stage_inferred, proportion, color = assigned_subtype, group = assigned_subtype)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  scale_color_manual(values = subtype_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA), expand = expansion(c(0, 0.08))) +
  labs(x = NULL, y = "Fraction of malignant cells", color = NULL, title = "Subtype dynamics across tumorigenesis stages") +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_pdf(pD, "D_小鼠四亚型占比趋势线.pdf", 5.2, 4.0)

if (file.exists(sanity_path)) {
  sanity <- read_csv(sanity_path, show_col_types = FALSE)
  feats <- split(sanity$mouse_symbol, factor(sanity$short_label, levels = subtype_levels))
  feats <- lapply(feats, function(g) intersect(g, rownames(mouse_obj)))
  feats <- feats[vapply(feats, length, integer(1)) > 0]
  pE <- DotPlot(mouse_obj, features = feats, group.by = "assigned_subtype") +
    base_theme +
    labs(x = NULL, y = NULL, title = "Sanity check: human-derived markers per assigned subtype") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  save_pdf(pE, "E_marker_sanity_dotplot.pdf", 7.5, 3.6)
}

pF1 <- ggplot(per_cell, aes(assigned_subtype_score, fill = assigned_subtype)) +
  geom_histogram(bins = 50, color = NA, alpha = 0.85) +
  scale_fill_manual(values = subtype_colors, guide = "none") +
  facet_wrap(~assigned_subtype, scales = "free_y", ncol = 4) +
  labs(x = "Top-1 subtype score (z)", y = "Cells", title = "Assignment score distribution") +
  base_theme
save_pdf(pF1, "F1_投射score分布.pdf", 8.0, 2.6)

pF2 <- ggplot(per_cell, aes(assigned_subtype_margin, fill = assignment_confidence)) +
  geom_histogram(bins = 60, color = NA, alpha = 0.85, position = "identity") +
  scale_fill_manual(values = c("high" = "#4D4D4D", "low" = "#D62728")) +
  labs(x = "Top1 - Top2 margin (z)", y = "Cells", fill = "Confidence", title = "Assignment margin (low-confidence = bottom decile)") +
  base_theme
save_pdf(pF2, "F2_投射margin分布.pdf", 5.0, 3.4)

if (file.exists(prop_hc_path)) {
  prop_hc <- read_csv(prop_hc_path, show_col_types = FALSE) %>%
    mutate(
      stage_inferred = factor(stage_inferred, levels = stage_levels),
      assigned_subtype = factor(assigned_subtype, levels = subtype_levels)
    )
  pG <- ggplot(prop_hc, aes(stage_inferred, proportion, fill = assigned_subtype)) +
    geom_col(width = 0.74, color = "white", linewidth = 0.25) +
    scale_fill_manual(values = subtype_colors, drop = FALSE) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0, 0)) +
    labs(x = NULL, y = "Fraction (high-confidence cells)", fill = NULL, title = "Sensitivity: high-confidence cells only") +
    base_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_pdf(pG, "G_高置信子集堆叠比例.pdf", 5.0, 4.2)
}

writeLines(
  capture.output(sessionInfo()),
  file.path(logs_dir, "03_小鼠阶段占比可视化_sessionInfo.txt")
)
message("== plotting done: figures written to ", figures_dir, " ==")
