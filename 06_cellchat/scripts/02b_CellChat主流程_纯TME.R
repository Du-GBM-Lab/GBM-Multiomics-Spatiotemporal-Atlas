#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(CellChat)
  library(dplyr)
  library(readr)
})

module_dir <- "07_细胞通讯"
outputs_dir <- file.path(module_dir, "outputs")
tables_dir <- file.path(module_dir, "tables")
figures_dir <- file.path(module_dir, "figures")
source_data_dir <- file.path(figures_dir, "source_data")
logs_dir <- file.path(module_dir, "logs")

dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

obj_path <- file.path(outputs_dir, "CellChat输入_四亚型加纯TME.qs2")
stopifnot(file.exists(obj_path))

checkpoint_overexpressed <- file.path(outputs_dir, "CellChat中间对象_overexpressed_纯TME.qs2")

if (file.exists(checkpoint_overexpressed)) {
  cellchat <- qs2::qs_read(checkpoint_overexpressed)
} else {
  obj <- qs2::qs_read(obj_path)
  DefaultAssay(obj) <- "RNA"
  obj$samples <- factor(obj$Pt_number)

  cellchat <- createCellChat(
    object = obj,
    group.by = "comm_group",
    assay = "RNA"
  )

  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  qs2::qs_save(cellchat, checkpoint_overexpressed)
}

cellchat <- computeCommunProb(
  cellchat,
  type = "triMean",
  raw.use = TRUE,
  population.size = TRUE
)
min_cells <- 50
cellchat <- filterCommunication(cellchat, min.cells = min_cells)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

qs2::qs_save(cellchat, file.path(outputs_dir, "CellChat结果_四亚型加纯TME.qs2"))

comm_all <- subsetCommunication(cellchat)
write_csv(comm_all, file.path(tables_dir, "CellChat纯TME_全部通讯结果.csv"))

pathway_rank <- comm_all %>%
  group_by(pathway_name) %>%
  summarise(
    n_pairs = n(),
    total_prob = sum(prob, na.rm = TRUE),
    mean_prob = mean(prob, na.rm = TRUE),
    min_pval = min(pval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_prob), min_pval)
write_csv(pathway_rank, file.path(tables_dir, "CellChat纯TME_pathway无偏排序.csv"))

pair_pathway_rank <- comm_all %>%
  group_by(source, target, pathway_name) %>%
  summarise(
    n_pairs = n(),
    total_prob = sum(prob, na.rm = TRUE),
    mean_prob = mean(prob, na.rm = TRUE),
    min_pval = min(pval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(source, target, desc(total_prob), min_pval)
write_csv(pair_pathway_rank, file.path(tables_dir, "CellChat纯TME_source_target_pathway无偏排序.csv"))

group_size <- as.data.frame(table(cellchat@idents), stringsAsFactors = FALSE)
colnames(group_size) <- c("comm_group", "n_cells")
write_csv(group_size, file.path(tables_dir, "CellChat纯TME_分组细胞数.csv"))

parameter_log <- data.frame(
  parameter = c("database", "computeCommunProb_type", "raw.use", "population.size", "min.cells"),
  value = c("CellChatDB.human_all_categories", "triMean", "TRUE", "TRUE", as.character(min_cells))
)
write_csv(parameter_log, file.path(tables_dir, "CellChat纯TME_主流程参数.csv"))

plaur_comm <- comm_all %>%
  filter(ligand %in% c("PLAU", "PLAUR") | receptor %in% c("PLAU", "PLAUR"))
write_csv(plaur_comm, file.path(tables_dir, "CellChat纯TME_PLAUR相关通讯.csv"))

pdf(file.path(figures_dir, "CellChat纯TME_总体通讯强度热图.pdf"), width = 7, height = 6)
netVisual_heatmap(cellchat, measure = "weight")
dev.off()

pdf(file.path(figures_dir, "CellChat纯TME_总体通讯数量热图.pdf"), width = 7, height = 6)
netVisual_heatmap(cellchat, measure = "count")
dev.off()

writeLines(
  capture.output(sessionInfo()),
  file.path(logs_dir, "02b_CellChat主流程_纯TME_sessionInfo.txt")
)
