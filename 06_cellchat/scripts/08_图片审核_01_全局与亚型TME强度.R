#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(circlize)
  library(scales)
})

archive_dir <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/06_细胞通讯"
main_dir <- file.path(archive_dir, "正文图")
supp_dir <- file.path(archive_dir, "补充图")
source_dir <- file.path(archive_dir, "source_data")

dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

weight_path <- file.path(source_dir, "CellChat纯TME_global_weight_matrix.csv")
count_path <- file.path(source_dir, "CellChat纯TME_global_count_matrix.csv")
stopifnot(file.exists(weight_path), file.exists(count_path))

malignant_groups <- c(
  "Malignant_NPC-P", "Malignant_OPC-M", "Malignant_MES-V", "Malignant_MES-I"
)
tme_groups <- c(
  "Macrophages", "Microglial", "Monocytes", "cDCs", "pDCs",
  "T cells", "NK cells", "B cells", "Endothelial", "Mural cells"
)
group_order <- c(malignant_groups, tme_groups)

display_labels <- c(
  "Malignant_NPC-P" = "NPC-P",
  "Malignant_OPC-M" = "OPC-M",
  "Malignant_MES-V" = "MES-V",
  "Malignant_MES-I" = "MES-I",
  "Macrophages" = "Macrophages",
  "Microglial" = "Microglia",
  "Monocytes" = "Monocytes",
  "cDCs" = "cDCs",
  "pDCs" = "pDCs",
  "T cells" = "T cells",
  "NK cells" = "NK cells",
  "B cells" = "B cells",
  "Endothelial" = "Endothelial",
  "Mural cells" = "Mural"
)

group_colors <- c(
  "Malignant_NPC-P" = "#0072B5",
  "Malignant_OPC-M" = "#E18727",
  "Malignant_MES-V" = "#20854E",
  "Malignant_MES-I" = "#BC3C29",
  "Macrophages" = "#7E3F98",
  "Microglial" = "#9A6BB2",
  "Monocytes" = "#B58AC8",
  "cDCs" = "#C9A9D8",
  "pDCs" = "#DCC6E6",
  "T cells" = "#4DBBD5",
  "NK cells" = "#00A087",
  "B cells" = "#91D1C2",
  "Endothelial" = "#8491B4",
  "Mural cells" = "#3C5488"
)

read_matrix_long <- function(path, value_name) {
  read_csv(path, show_col_types = FALSE) %>%
    rename(source = 1) %>%
    pivot_longer(-source, names_to = "target", values_to = value_name)
}

weight_long <- read_matrix_long(weight_path, "weight")
count_long <- read_matrix_long(count_path, "count")

global_summary <- weight_long %>%
  filter(source %in% group_order, target %in% group_order) %>%
  group_by(group = source) %>%
  summarise(outgoing_weight = sum(weight, na.rm = TRUE), .groups = "drop") %>%
  full_join(
    weight_long %>%
      filter(source %in% group_order, target %in% group_order) %>%
      group_by(group = target) %>%
      summarise(incoming_weight = sum(weight, na.rm = TRUE), .groups = "drop"),
    by = "group"
  ) %>%
  mutate(
    total_weight = outgoing_weight + incoming_weight,
    group_label = unname(display_labels[group])
  ) %>%
  arrange(desc(total_weight))

write_csv(global_summary, file.path(source_dir, "全局交互强度汇总.csv"))

mal_tme_edges <- weight_long %>%
  filter(
    (source %in% malignant_groups & target %in% tme_groups) |
      (source %in% tme_groups & target %in% malignant_groups)
  ) %>%
  mutate(
    source_label = unname(display_labels[source]),
    target_label = unname(display_labels[target]),
    direction_class = case_when(
      source %in% malignant_groups ~ "malignant_to_TME",
      TRUE ~ "TME_to_malignant"
    )
  ) %>%
  arrange(desc(weight)) %>%
  mutate(rank = row_number())

write_csv(mal_tme_edges, file.path(source_dir, "四亚型TME交互强度弦图_source.csv"))

chord_df <- mal_tme_edges %>%
  transmute(from = source, to = target, value = weight)

link_cols <- alpha(group_colors[chord_df$from], rescale(chord_df$value, to = c(0.22, 0.82)))

pdf(file.path(main_dir, "待审核_四亚型TME交互强度弦图.pdf"), width = 8.2, height = 8.0, useDingbats = FALSE)
circos.clear()
circos.par(
  start.degree = 90,
  gap.after = c(rep(2, length(malignant_groups) - 1), 8, rep(2, length(tme_groups) - 1), 8),
  track.margin = c(0.01, 0.01),
  points.overflow.warning = FALSE
)
chordDiagram(
  x = chord_df,
  order = group_order,
  grid.col = group_colors[group_order],
  col = link_cols,
  transparency = 0.35,
  directional = 1,
  direction.type = c("arrows", "diffHeight"),
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.largest.ontop = TRUE,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.13)
)
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    circos.text(
      x = mean(xlim),
      y = ylim[1] + 0.1,
      labels = display_labels[[sector]],
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.72
    )
  },
  bg.border = NA
)
title("Interaction strength between malignant subtypes and TME cells", cex.main = 1)
circos.clear()
dev.off()

top_malignant_tme <- mal_tme_edges %>%
  slice_head(n = 25)
write_csv(top_malignant_tme, file.path(source_dir, "四亚型TME交互强度_top25.csv"))
