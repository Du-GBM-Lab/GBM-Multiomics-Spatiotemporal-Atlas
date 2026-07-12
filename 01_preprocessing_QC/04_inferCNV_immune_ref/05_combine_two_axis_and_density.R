suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

project_root <- "<DATA_ROOT>/项目/分型/修稿杠生信"
work_dir <- file.path(project_root, "重新分析/04_inferCNV_免疫参考验证")
source_dir <- file.path(
  project_root,
  "图片表格/02_恶性筛选（没有过一遍）/源数据/补图_inferCNV双指标与CNV负荷分布"
)
out_pdf <- file.path(
  project_root,
  "图片表格/02_恶性筛选（没有过一遍）/补充图/补图_inferCNV双指标与CNV负荷分布.pdf"
)
out_png <- file.path(source_dir, "补图_inferCNV双指标与CNV负荷分布_preview.png")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

score_path <- file.path(
  work_dir,
  "tables/infercnv_immune_reference_cell_cnv_burden_and_calls.csv"
)
threshold_path <- file.path(
  work_dir,
  "tables/infercnv_immune_reference_samplewise_thresholds.csv"
)

scores <- fread(score_path)
thresholds <- fread(threshold_path)

sample_levels <- paste0("Pt", c(1:24))
scores[, sample := factor(sample, levels = sample_levels)]
thresholds[, sample := factor(sample, levels = sample_levels)]

call_levels <- c(
  "non_malignant_reference",
  "non_malignant_like_CNV_low",
  "malignant_like_CNV_burden_only",
  "malignant_like_CNV_high_confidence"
)
call_labels <- c(
  "Reference",
  "CNV-low",
  "CNV-burden only",
  "High-confidence malignant"
)
call_palette <- c(
  "non_malignant_reference" = "#6B8798",
  "non_malignant_like_CNV_low" = "#91B8C9",
  "malignant_like_CNV_burden_only" = "#E58A57",
  "malignant_like_CNV_high_confidence" = "#B9182B"
)

scores[, infercnv_call := factor(infercnv_call, levels = call_levels)]

theme_small <- theme_bw(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.18, color = "#E4E4E4"),
    strip.background = element_rect(fill = "#F0F0F0", color = "#BEBEBE", linewidth = 0.25),
    strip.text = element_text(face = "bold", size = 6.7, margin = margin(1.5, 0, 1.5, 0)),
    axis.text = element_text(size = 5.8, color = "#303030"),
    axis.title = element_text(size = 7.4),
    panel.spacing = unit(0.16, "cm"),
    plot.margin = margin(2, 2, 2, 2)
  )

p_two_axis <- ggplot(
  scores[!is.na(cnv_burden_z) & !is.na(cnv_correlation_ref_z)],
  aes(cnv_burden_z, cnv_correlation_ref_z, color = infercnv_call)
) +
  geom_point(size = 0.18, alpha = 0.48, stroke = 0) +
  geom_vline(xintercept = 3, linetype = "dashed", linewidth = 0.25, color = "#202020") +
  geom_hline(yintercept = 3, linetype = "dashed", linewidth = 0.25, color = "#202020") +
  facet_wrap(~sample, scales = "free", ncol = 6, drop = FALSE) +
  scale_color_manual(
    values = call_palette,
    breaks = call_levels,
    labels = call_labels,
    drop = FALSE
  ) +
  labs(
    x = "CNV burden z-score",
    y = "CNV correlation z-score",
    color = NULL
  ) +
  theme_small +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6.5),
    legend.key.width = unit(0.35, "cm"),
    legend.spacing.x = unit(0.08, "cm")
  ) +
  guides(color = guide_legend(
    nrow = 1,
    override.aes = list(size = 1.8, alpha = 1)
  ))

density_df <- scores[
  !is.na(cnv_burden_z),
  .(
    cnv_burden_z,
    density_group = ifelse(is_reference_group, "Reference", "Non-reference"),
    sample
  )
]

density_palette <- c(
  "Reference" = "#8EA6B3",
  "Non-reference" = "#C66A67"
)

p_density <- ggplot(
  density_df,
  aes(cnv_burden_z, fill = density_group, color = density_group)
) +
  geom_density(alpha = 0.38, linewidth = 0.28, adjust = 0.85) +
  geom_vline(xintercept = 3, linetype = "dashed", linewidth = 0.25, color = "#202020") +
  facet_wrap(~sample, scales = "free", ncol = 6, drop = FALSE) +
  scale_fill_manual(values = density_palette) +
  scale_color_manual(values = density_palette) +
  labs(
    x = "CNV burden z-score",
    y = "Density",
    fill = NULL,
    color = NULL
  ) +
  theme_small +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6.5),
    legend.key.width = unit(0.42, "cm")
  ) +
  guides(
    fill = guide_legend(nrow = 1),
    color = "none"
  )

combined <- p_two_axis / p_density +
  plot_layout(heights = c(1.55, 1), guides = "keep") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 11),
      plot.tag.position = c(0.005, 0.995)
    )
  )

ggsave(out_pdf, combined, width = 11.5, height = 10.2, units = "in", device = cairo_pdf, bg = "white")
ggsave(out_png, combined, width = 11.5, height = 10.2, units = "in", dpi = 220, bg = "white")

fwrite(scores, file.path(source_dir, "infercnv_immune_reference_cell_cnv_burden_and_calls.csv"))
fwrite(thresholds, file.path(source_dir, "infercnv_immune_reference_samplewise_thresholds.csv"))
capture.output(sessionInfo(), file = file.path(source_dir, "sessionInfo.txt"))
