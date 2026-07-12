#!/usr/bin/env Rscript
# =============================================================================
# 07_发育时间_TF_ATAC验证 / 01_mouse_developmental
# Phase 1 - Script 04: 论文用最终图 (仅三张)
# -----------------------------------------------------------------------------
#   A. 投射亚型 UMAP                          -> 正文
#   E. 四亚型 canonical marker dotplot          -> 正文候选
#   R. 同心圆亚型组成 + 右上角占比变化曲线       -> 正文/graphical 概览
#
# 关于同心圆:
#   每个环固定整圈 360 度，扇区角度只编码亚型比例。弧长不编码每阶段细胞总数，
#   避免把 Preneoplastic 455 / Early 1818 / Mid 5647 / Endpoint 11097 这类
#   测序捕获量差异画成主要视觉变量。
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
prop_path <- file.path(src_data_dir, "mouse_phase1_stage_subtype_proportions.csv")
sanity_path <- file.path(tables_dir, "mouse_phase1_sanity_markers.csv")
stopifnot(file.exists(obj_path), file.exists(prop_path))

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

mouse_obj <- qs2::qs_read(obj_path)

prop <- read_csv(prop_path, show_col_types = FALSE) %>%
  mutate(
    stage_inferred = factor(stage_inferred, levels = stage_levels),
    assigned_subtype = factor(assigned_subtype, levels = subtype_levels)
  )

# =============================================================================
# A. 投射亚型 UMAP
# =============================================================================
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
  labs(title = "Cross-species assigned subtype (mouse malignant cells)") +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  )
save_pdf(pA, "A_小鼠四亚型投射UMAP.pdf", 5.4, 4.4)

# =============================================================================
# E. canonical marker dotplot
# =============================================================================
if (file.exists(sanity_path)) {
  sanity <- read_csv(sanity_path, show_col_types = FALSE)
  feats <- split(
    sanity$mouse_symbol,
    factor(sanity$short_label, levels = subtype_levels)
  )
  feats <- lapply(feats, function(g) intersect(g, rownames(mouse_obj)))
  feats <- feats[vapply(feats, length, integer(1)) > 0]

  pE <- DotPlot(mouse_obj, features = feats, group.by = "assigned_subtype") +
    base_theme +
    labs(x = NULL, y = NULL, title = "Human-derived subtype markers in mouse malignant cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  save_pdf(pE, "E_marker_sanity_dotplot.pdf", 7.6, 3.6)
} else {
  message("[warn] sanity marker table not found, skip dotplot.")
}

# =============================================================================
# R. 同心圆亚型组成 + 右上角占比变化曲线
# =============================================================================
inner_radius <- 1.0
ring_width <- 0.85
ring_gap <- 0.06

ring_df <- prop %>%
  arrange(stage_inferred, assigned_subtype) %>%
  group_by(stage_inferred) %>%
  mutate(
    theta_max = cumsum(proportion),
    theta_min = theta_max - proportion
  ) %>%
  ungroup() %>%
  mutate(
    ring_idx = as.integer(stage_inferred),
    r_min = inner_radius + (ring_idx - 1) * (ring_width + ring_gap),
    r_max = r_min + ring_width,
    theta_mid = (theta_min + theta_max) / 2,
    r_mid = (r_min + r_max) / 2,
    pct_lab = ifelse(proportion >= 0.08, paste0(round(proportion * 100), "%"), "")
  )

stage_lab_df <- ring_df %>%
  distinct(stage_inferred, r_mid_ring = (r_min + r_max) / 2) %>%
  group_by(stage_inferred) %>%
  summarise(r = mean(r_mid_ring), .groups = "drop") %>%
  mutate(
    stage_n = prop %>%
      distinct(stage_inferred, stage_total) %>%
      arrange(stage_inferred) %>%
      pull(stage_total),
    label = paste0(stage_inferred, "  (n=", comma(stage_n), ")")
  )

r_outer <- max(ring_df$r_max)

pRing <- ggplot(ring_df) +
  geom_rect(
    aes(
      xmin = theta_min,
      xmax = theta_max,
      ymin = r_min,
      ymax = r_max,
      fill = assigned_subtype
    ),
    color = "white",
    linewidth = 0.35
  ) +
  geom_text(
    aes(x = theta_mid, y = r_mid, label = pct_lab),
    size = 2.5,
    color = "white",
    fontface = "bold"
  ) +
  geom_text(
    data = stage_lab_df,
    aes(x = 0, y = r, label = label),
    inherit.aes = FALSE,
    hjust = 1.05,
    size = 2.8,
    fontface = "bold"
  ) +
  scale_fill_manual(values = subtype_colors, name = NULL) +
  coord_polar(theta = "x", start = 0, direction = 1) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, r_outer + 0.2), expand = c(0, 0)) +
  labs(title = "Subtype composition across tumorigenesis stages") +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    legend.position = "left",
    legend.key.size = unit(11, "pt")
  )

inset <- ggplot(
  prop,
  aes(stage_inferred, proportion, color = assigned_subtype, group = assigned_subtype)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_color_manual(values = subtype_colors, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA), expand = expansion(c(0, 0.08))) +
  scale_x_discrete(labels = c("Pre", "Early", "Mid", "End")) +
  labs(x = NULL, y = "Fraction") +
  theme_classic(base_size = 8) +
  theme(
    axis.text = element_text(color = "black", size = 7),
    axis.title.y = element_text(size = 7.5),
    plot.background = element_rect(fill = "white", color = NA)
  )

if (requireNamespace("cowplot", quietly = TRUE)) {
  pR <- cowplot::ggdraw(pRing) +
    cowplot::draw_plot(inset, x = 0.62, y = 0.62, width = 0.36, height = 0.34)
} else {
  message("[warn] 'cowplot' not installed, saving inset as separate PDF.")
  pR <- pRing
  save_pdf(inset + theme(legend.position = "none"), "R_inset_占比曲线.pdf", 3.6, 2.8)
}

save_pdf(pR, "R_同心圆亚型组成.pdf", 7.2, 6.2)

writeLines(
  capture.output(sessionInfo()),
  file.path(logs_dir, "04_论文用最终图_sessionInfo.txt")
)
message("== final figures written to ", figures_dir, " ==")
