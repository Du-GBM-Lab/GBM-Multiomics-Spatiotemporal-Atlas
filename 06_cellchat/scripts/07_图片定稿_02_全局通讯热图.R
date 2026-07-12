#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(cowplot)
  library(scales)
})

module_dir <- "07_细胞通讯"
figures_dir <- file.path(module_dir, "figures")
source_data_dir <- file.path(figures_dir, "source_data")
archive_dir <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/06_细胞通讯"
archive_supp_dir <- file.path(archive_dir, "补充候选")
archive_source_dir <- file.path(archive_dir, "source_data")

dir.create(archive_supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(archive_source_dir, recursive = TRUE, showWarnings = FALSE)

count_path <- file.path(source_data_dir, "CellChat纯TME_global_count_matrix.csv")
weight_path <- file.path(source_data_dir, "CellChat纯TME_global_weight_matrix.csv")

count_mat <- read_csv(count_path, show_col_types = FALSE)
weight_mat <- read_csv(weight_path, show_col_types = FALSE)

group_order <- c(
  "Malignant_NPC-P", "Malignant_OPC-M", "Malignant_MES-V", "Malignant_MES-I",
  "Macrophages", "Microglial", "Monocytes", "cDCs", "pDCs",
  "T cells", "NK cells", "B cells",
  "Endothelial", "Mural cells"
)

to_long <- function(tbl, value_name) {
  tbl %>%
    pivot_longer(-source, names_to = "target", values_to = value_name) %>%
    mutate(
      source = factor(source, levels = rev(group_order)),
      target = factor(target, levels = group_order)
    ) %>%
    filter(!is.na(source), !is.na(target))
}

count_long <- to_long(count_mat, "count")
weight_long <- to_long(weight_mat, "weight")

write_csv(count_long, file.path(source_data_dir, "全局通讯数量热图_source.csv"))
write_csv(weight_long, file.path(source_data_dir, "全局通讯强度热图_source.csv"))
write_csv(count_long, file.path(archive_source_dir, "全局通讯数量热图_source.csv"))
write_csv(weight_long, file.path(archive_source_dir, "全局通讯强度热图_source.csv"))

heat_theme <- theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6.5, color = "black"),
    axis.text.y = element_text(size = 6.5, color = "black"),
    axis.title = element_text(size = 8),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(size = 10, face = "bold", hjust = 0),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.5),
    plot.margin = margin(4, 4, 4, 4)
  )

p_count <- ggplot(count_long, aes(target, source, fill = count)) +
  geom_tile(color = "white", linewidth = 0.22) +
  scale_fill_gradient(low = "white", high = "#9ECAE1", name = "Count") +
  labs(title = "Number of interactions", x = "Receiver", y = "Sender") +
  heat_theme

p_weight <- ggplot(weight_long, aes(target, source, fill = weight)) +
  geom_tile(color = "white", linewidth = 0.22) +
  scale_fill_gradient(
    low = "white",
    high = "#B2182B",
    name = "Weight",
    labels = label_number(accuracy = 0.001)
  ) +
  labs(title = "Interaction strength", x = "Receiver", y = NULL) +
  heat_theme +
  theme(axis.text.y = element_blank())

combined <- cowplot::plot_grid(
  p_count,
  p_weight,
  nrow = 1,
  rel_widths = c(1.05, 1)
)

out_pdf <- file.path(figures_dir, "[补充定稿]全局通讯数量与强度热图.pdf")
archive_pdf <- file.path(archive_supp_dir, "全局通讯数量与强度热图.pdf")

ggsave(out_pdf, combined, width = 8.6, height = 4.8, useDingbats = FALSE)
ggsave(archive_pdf, combined, width = 8.6, height = 4.8, useDingbats = FALSE)

file.copy(
  file.path(figures_dir, "CellChat纯TME_总体通讯数量热图.pdf"),
  file.path(archive_supp_dir, "待审核_CellChat纯TME_总体通讯数量热图.pdf"),
  overwrite = TRUE
)
file.copy(
  file.path(figures_dir, "CellChat纯TME_总体通讯强度热图.pdf"),
  file.path(archive_supp_dir, "待审核_CellChat纯TME_总体通讯强度热图.pdf"),
  overwrite = TRUE
)

message("Finalized global communication heatmap: ", archive_pdf)
