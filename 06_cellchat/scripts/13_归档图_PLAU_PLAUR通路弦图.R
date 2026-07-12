#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(qs2)
  library(dplyr)
  library(readr)
})

archive_dir <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/06_细胞通讯"
supp_fig_dir <- file.path(archive_dir, "补充图")
source_dir <- file.path(archive_dir, "source_data")
dir.create(supp_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

cellchat <- qs2::qs_read("07_细胞通讯/outputs/CellChat结果_四亚型加纯TME.qs2")

subtype_colors <- c(
  "Malignant_NPC-P" = "#0072B5",
  "Malignant_OPC-M" = "#E18727",
  "Malignant_MES-V" = "#20854E",
  "Malignant_MES-I" = "#BC3C29"
)
tme_colors <- c(
  "Macrophages" = "#6A3D9A",
  "Microglial" = "#8E63B0",
  "Monocytes" = "#CAB2D6",
  "cDCs" = "#FB9A99",
  "pDCs" = "#FDBF6F",
  "T cells" = "#1F78B4",
  "NK cells" = "#A6CEE3",
  "B cells" = "#33A02C",
  "Endothelial" = "#B15928",
  "Mural cells" = "#FF7F00"
)
cell_colors <- c(subtype_colors, tme_colors)

plau_lr <- extractEnrichedLR(
  cellchat,
  signaling = "PLAU",
  geneLR.return = FALSE,
  enriched.only = TRUE
)
write_csv(
  as_tibble(plau_lr),
  file.path(source_dir, "PLAU_PLAUR通路LR构成_source.csv")
)

plau_comm <- read_csv(
  "07_细胞通讯/tables/CellChat纯TME_PLAUR相关通讯_按prob排序.csv",
  show_col_types = FALSE
) %>%
  filter(ligand == "PLAU", receptor == "PLAUR")
write_csv(
  plau_comm,
  file.path(source_dir, "PLAU_PLAUR通路通讯边_source.csv")
)

pdf(
  file.path(supp_fig_dir, "PLAU_PLAUR通路弦图.pdf"),
  width = 7,
  height = 7,
  useDingbats = FALSE
)
netVisual_aggregate(
  cellchat,
  signaling = "PLAU",
  signaling.name = "PLAU-PLAUR pathway",
  layout = "chord",
  color.use = cell_colors[levels(cellchat@idents)]
)
dev.off()

