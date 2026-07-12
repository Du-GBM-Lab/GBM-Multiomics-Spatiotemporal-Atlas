#!/usr/bin/env Rscript
# =============================================================================
# 07_发育时间_TF_ATAC验证 / 01_mouse_developmental
# Phase 1 - Script 04f: 论文用同心圆图 (右侧双图例 + 环框 + 网格折线)
# -----------------------------------------------------------------------------
# 1. subtype legend 和 stage legend 都放右侧。
# 2. 每个同心圆环加灰阶外框，对应 stage legend。
# 3. 折线图使用浅灰网格背景。
# 4. 每个阶段固定 270 度，弧长不编码细胞总数。
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

if (!requireNamespace("cowplot", quietly = TRUE)) {
  stop("需要 cowplot: install.packages('cowplot')")
}

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
  "NPC-P" = "#1F6FB4",
  "OPC-M" = "#E08214",
  "MES-V" = "#1B7837",
  "MES-I" = "#B2182B"
)
stage_colors <- c(
  "Preneoplastic" = "#FFFFFF",
  "Early lesion" = "#FFFFFF",
  "Mid lesion" = "#FFFFFF",
  "Endpoint" = "#FFFFFF"
)
stage_border_colors <- c(
  "Preneoplastic" = "#7B61FF",
  "Early lesion" = "#00A6D6",
  "Mid lesion" = "#6BBF45",
  "Endpoint" = "#F05A28"
)

save_pdf <- function(plot, name, w, h) {
  ggsave(file.path(figures_dir, name), plot = plot, width = w, height = h, useDingbats = FALSE)
}

mouse_obj <- qs2::qs_read(obj_path)
prop <- read_csv(prop_path, show_col_types = FALSE) %>%
  mutate(
    stage_inferred = factor(stage_inferred, levels = stage_levels),
    assigned_subtype = factor(assigned_subtype, levels = subtype_levels)
  )

# A. UMAP
red_name <- intersect(c("umap", "UMAP", "Umap"), names(mouse_obj@reductions))[1]
if (is.na(red_name)) {
  stop("No UMAP reduction.")
}

pA <- DimPlot(
  mouse_obj,
  reduction = red_name,
  group.by = "assigned_subtype",
  pt.size = 0.25,
  cols = subtype_colors
) +
  theme_classic(base_size = 10) +
  labs(title = "Cross-species assigned subtype (mouse malignant cells)") +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    plot.title = element_text(size = 10, face = "bold")
  )
save_pdf(pA, "A_小鼠四亚型投射UMAP.pdf", 5.4, 4.4)

# E. marker dotplot
if (file.exists(sanity_path)) {
  sanity <- read_csv(sanity_path, show_col_types = FALSE)
  feats <- split(sanity$mouse_symbol, factor(sanity$short_label, levels = subtype_levels))
  feats <- lapply(feats, function(g) intersect(g, rownames(mouse_obj)))
  feats <- feats[vapply(feats, length, integer(1)) > 0]
  pE <- DotPlot(mouse_obj, features = feats, group.by = "assigned_subtype") +
    theme_classic(base_size = 10) +
    labs(x = NULL, y = NULL, title = "Human-derived subtype markers in mouse malignant cells") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      plot.title = element_text(size = 10, face = "bold")
    )
  save_pdf(pE, "E_marker_sanity_dotplot.pdf", 7.6, 3.6)
}

# R. 同心圆
arc_frac <- 0.75
inner_radius <- 1.35
ring_width <- 0.80
ring_gap <- 0.07

ring_df <- prop %>%
  arrange(stage_inferred, assigned_subtype) %>%
  group_by(stage_inferred) %>%
  mutate(
    cum_max = cumsum(proportion),
    cum_min = cum_max - proportion,
    theta_min = cum_min * arc_frac,
    theta_max = cum_max * arc_frac,
    theta_mid = (theta_min + theta_max) / 2
  ) %>%
  ungroup() %>%
  mutate(
    ring_idx = as.integer(stage_inferred),
    r_min = inner_radius + (ring_idx - 1) * (ring_width + ring_gap),
    r_max = r_min + ring_width,
    r_mid = (r_min + r_max) / 2,
    pct_lab = ifelse(proportion >= 0.09, paste0(round(proportion * 100), "%"), "")
  )

outline_df <- ring_df %>%
  distinct(stage_inferred, r_min, r_max) %>%
  mutate(theta_min = 0, theta_max = arc_frac)

r_outer <- max(ring_df$r_max)

pRing_base <- ggplot(ring_df) +
  geom_rect(
    aes(xmin = theta_min, xmax = theta_max, ymin = r_min, ymax = r_max, fill = assigned_subtype),
    color = "white",
    linewidth = 0.65
  ) +
  geom_rect(
    data = outline_df,
    aes(xmin = theta_min, xmax = theta_max, ymin = r_min, ymax = r_max, color = stage_inferred),
    fill = NA,
    linewidth = 0.95,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = theta_mid, y = r_mid, label = pct_lab),
    size = 2.7,
    color = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(values = subtype_colors, name = NULL) +
  scale_color_manual(values = stage_border_colors) +
  coord_polar(theta = "x", start = 0, direction = -1, clip = "off") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, r_outer + 0.34), expand = c(0, 0)) +
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  theme_void(base_size = 10) +
  theme(
    plot.margin = margin(4, 4, 4, 4),
    legend.position = "none"
  )

pRing <- pRing_base

inset <- ggplot(
  prop,
  aes(stage_inferred, proportion, color = assigned_subtype, group = assigned_subtype)
) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.7) +
  scale_color_manual(values = subtype_colors, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA), expand = expansion(c(0, 0.10))) +
  scale_x_discrete(labels = c("Pre", "Early", "Mid", "End")) +
  labs(x = NULL, y = "Fraction") +
  theme_bw(base_size = 7.5) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey94", linewidth = 0.25),
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.45),
    axis.text = element_text(color = "black", size = 7),
    axis.title.y = element_text(size = 7.5),
    axis.line = element_line(linewidth = 0.45, color = "black"),
    axis.ticks = element_line(linewidth = 0.45, color = "black"),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(5, 5, 5, 5)
  )

stage_legend_df <- prop %>%
  distinct(stage_inferred, stage_total) %>%
  arrange(stage_inferred) %>%
  mutate(
    stage_inferred = factor(stage_inferred, levels = stage_levels),
    stage_label = paste0(stage_inferred, " (n=", comma(stage_total), ")"),
    y = 4.7 - seq_along(stage_levels),
    border_color = unname(stage_border_colors[as.character(stage_inferred)])
  )

subtype_legend_df <- tibble::tibble(
  label = subtype_levels,
  fill_color = unname(subtype_colors[subtype_levels]),
  x = 1,
  y = 10.2 - seq_along(subtype_levels)
)

pRightLegend <- ggplot() +
  geom_text(
    data = tibble::tibble(x = 1, y = 10.1, label = "Subtype"),
    aes(x = x, y = y, label = label),
    hjust = 0,
    size = 3.6,
    fontface = "bold",
    color = "grey15"
  ) +
  geom_tile(
    data = subtype_legend_df,
    aes(x = x, y = y, fill = label),
    width = 0.18,
    height = 0.36,
    color = NA
  ) +
  geom_text(
    data = subtype_legend_df,
    aes(x = 1.23, y = y, label = label),
    hjust = 0,
    size = 4.2,
    fontface = "bold",
    color = "black"
  ) +
  geom_text(
    data = tibble::tibble(x = 1, y = 5.25, label = "Stage (inner to outer)"),
    aes(x = x, y = y, label = label),
    hjust = 0,
    size = 3.4,
    fontface = "bold",
    color = "grey15"
  ) +
  geom_tile(
    data = stage_legend_df,
    aes(x = 1, y = y, fill = stage_inferred, color = stage_inferred),
    width = 0.18,
    height = 0.40,
    linewidth = 0.9
  ) +
  geom_text(
    data = stage_legend_df,
    aes(x = 1.23, y = y, label = stage_label),
    hjust = 0,
    size = 3.1,
    fontface = "bold",
    color = "grey20"
  ) +
  scale_fill_manual(values = c(subtype_colors, stage_colors), guide = "none") +
  scale_color_manual(values = stage_border_colors, guide = "none") +
  scale_x_continuous(limits = c(0.75, 3.85), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.2, 10.6), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = 10) +
  theme(plot.margin = margin(0, 0, 0, 0))

title_gg <- cowplot::ggdraw() +
  cowplot::draw_label(
    "Subtype composition across tumorigenesis stages",
    fontface = "bold",
    size = 12,
    x = 0.5,
    hjust = 0.5
  )

ring_with_inset <- cowplot::ggdraw(pRing) +
  cowplot::draw_plot(inset, x = 0.60, y = 0.58, width = 0.35, height = 0.32)

body <- cowplot::plot_grid(
  ring_with_inset,
  pRightLegend,
  nrow = 1,
  rel_widths = c(1, 0.36)
)

pR <- cowplot::plot_grid(title_gg, body, ncol = 1, rel_heights = c(0.08, 1))
save_pdf(pR, "R_同心圆亚型组成.pdf", 9.2, 6.8)

writeLines(capture.output(sessionInfo()), file.path(logs_dir, "04f_论文用同心圆图_sessionInfo.txt"))
message("== figures written to ", figures_dir, " ==")
