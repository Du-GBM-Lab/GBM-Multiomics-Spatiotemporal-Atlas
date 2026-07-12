#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(grid)
})

archive_dir <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/06_细胞通讯"
source_dir <- file.path(archive_dir, "source_data")
main_fig_dir <- file.path(archive_dir, "正文图")

dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(main_fig_dir, recursive = TRUE, showWarnings = FALSE)

expr_path <- file.path(source_dir, "PLAU_PLAUR表达per_cell_source.csv")
stopifnot(file.exists(expr_path))

subtype_levels <- c("NPC-P", "OPC-M", "MES-V", "MES-I")
subtype_colors <- c(
  "NPC-P" = "#0072B5",
  "OPC-M" = "#E18727",
  "MES-V" = "#20854E",
  "MES-I" = "#BC3C29"
)

expr <- read_csv(expr_path, show_col_types = FALSE) %>%
  filter(gene %in% c("PLAU", "PLAUR")) %>%
  mutate(
    subtype = case_when(
      comm_group == "Malignant_NPC-P" ~ "NPC-P",
      comm_group == "Malignant_OPC-M" ~ "OPC-M",
      comm_group == "Malignant_MES-V" ~ "MES-V",
      comm_group == "Malignant_MES-I" ~ "MES-I",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(subtype)) %>%
  mutate(subtype = factor(subtype, levels = subtype_levels))

summary_tbl <- expr %>%
  group_by(gene, subtype) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(Pt_number),
    pct_detected = mean(detected) * 100,
    mean_expr = mean(expr),
    median_expr = median(expr),
    q75_expr = quantile(expr, 0.75),
    .groups = "drop"
  )

plaur_summary <- summary_tbl %>%
  filter(gene == "PLAUR") %>%
  mutate(subtype = factor(subtype, levels = subtype_levels))

write_csv(
  plaur_summary,
  file.path(source_dir, "四亚型PLAUR表达比较_source.csv")
)

p_plaur <- ggplot(
  plaur_summary,
  aes(x = subtype, y = mean_expr, color = subtype, fill = subtype)
) +
  geom_col(
    width = 0.58,
    alpha = 0.18,
    color = NA
  ) +
  geom_point(
    aes(size = pct_detected),
    shape = 21,
    stroke = 0.35
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct_detected)),
    y = plaur_summary$mean_expr + 0.08,
    size = 2.6,
    color = "grey20"
  ) +
  scale_color_manual(values = subtype_colors, drop = FALSE) +
  scale_fill_manual(values = subtype_colors, drop = FALSE) +
  scale_size_continuous(
    range = c(2.4, 6.4),
    breaks = c(15, 30, 45, 60),
    limits = c(0, 60),
    name = "Detected cells (%)"
  ) +
  scale_y_continuous(
    breaks = pretty_breaks(n = 4),
    expand = expansion(mult = c(0, 0.10))
  ) +
  labs(
    x = NULL,
    y = "Mean log-normalized expression",
    title = "PLAUR expression across malignant subtypes"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.key.height = unit(0.32, "cm"),
    legend.key.width = unit(0.32, "cm"),
    plot.title = element_text(face = "bold", size = 10, hjust = 0),
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.25)
  )

ggsave(
  filename = file.path(main_fig_dir, "待审核_PLAUR四亚型表达比较.pdf"),
  plot = p_plaur,
  width = 3.25,
  height = 2.45,
  useDingbats = FALSE
)
