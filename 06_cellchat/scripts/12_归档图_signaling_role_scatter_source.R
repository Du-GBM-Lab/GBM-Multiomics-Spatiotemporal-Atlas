#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(qs2)
  library(CellChat)
  library(dplyr)
  library(readr)
})

archive_dir <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/06_细胞通讯"
source_dir <- file.path(archive_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

cellchat <- qs2::qs_read("07_细胞通讯/outputs/CellChat结果_四亚型加纯TME.qs2")

centr <- cellchat@netP$centr
outgoing <- matrix(0, nrow = nlevels(cellchat@idents), ncol = length(centr))
incoming <- matrix(0, nrow = nlevels(cellchat@idents), ncol = length(centr))
dimnames(outgoing) <- list(levels(cellchat@idents), names(centr))
dimnames(incoming) <- dimnames(outgoing)

for (i in seq_along(centr)) {
  outgoing[, i] <- centr[[i]]$outdeg
  incoming[, i] <- centr[[i]]$indeg
}

num_link <- CellChat::aggregateNet(
  cellchat,
  signaling = NULL,
  return.object = FALSE,
  remove.isolate = FALSE
)$count
num_link <- rowSums(num_link) + colSums(num_link) - diag(num_link)

role_tbl <- tibble(
  comm_group = names(num_link),
  outgoing_interaction_strength = rowSums(outgoing)[names(num_link)],
  incoming_interaction_strength = rowSums(incoming)[names(num_link)],
  significant_interaction_count = as.numeric(num_link)
) %>%
  arrange(desc(outgoing_interaction_strength + incoming_interaction_strength))

write_csv(
  role_tbl,
  file.path(source_dir, "细胞通讯角色散点图_source.csv")
)

