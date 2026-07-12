#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(CellChat)
  library(qs2)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

module_dir <- "07_细胞通讯"
outputs_dir <- file.path(module_dir, "outputs")
tables_dir <- file.path(module_dir, "tables")
figures_dir <- file.path(module_dir, "figures")
source_data_dir <- file.path(figures_dir, "source_data")
logs_dir <- file.path(module_dir, "logs")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

pure_obj_path <- file.path(outputs_dir, "CellChat结果_四亚型加纯TME.qs2")
full_obj_path <- file.path(outputs_dir, "CellChat结果_四亚型加TME.qs2")
stopifnot(file.exists(pure_obj_path), file.exists(full_obj_path))

cellchat <- qs2::qs_read(pure_obj_path)
cellchat_full <- qs2::qs_read(full_obj_path)

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

group_size <- as.numeric(table(cellchat@idents))
names(group_size) <- names(table(cellchat@idents))
write_csv(
  tibble(comm_group = names(group_size), n_cells = as.integer(group_size)),
  file.path(source_data_dir, "CellChat纯TME_group_size.csv")
)

save_pdf <- function(expr, filename, width = 7, height = 6) {
  pdf(file.path(figures_dir, filename), width = width, height = height, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  force(expr)
}

# Global interaction network.
write_csv(
  as.data.frame(cellchat@net$count) %>%
    tibble::rownames_to_column("source"),
  file.path(source_data_dir, "CellChat纯TME_global_count_matrix.csv")
)
write_csv(
  as.data.frame(cellchat@net$weight) %>%
    tibble::rownames_to_column("source"),
  file.path(source_data_dir, "CellChat纯TME_global_weight_matrix.csv")
)

save_pdf(
  netVisual_circle(
    cellchat@net$count,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    color.use = cell_colors[names(group_size)],
    title.name = "Number of interactions"
  ),
  "[补充]CellChat纯TME_聚合互作数量圈图.pdf",
  7,
  7
)

save_pdf(
  netVisual_circle(
    cellchat@net$weight,
    vertex.weight = group_size,
    weight.scale = TRUE,
    label.edge = FALSE,
    color.use = cell_colors[names(group_size)],
    title.name = "Interaction strength"
  ),
  "[补充]CellChat纯TME_聚合互作强度圈图.pdf",
  7,
  7
)

# Signaling role.
save_pdf(
  print(netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.heatmap = "Reds", title = "Outgoing signaling")),
  "[正文]CellChat纯TME_outgoing_signaling_role_heatmap.pdf",
  8,
  7
)

save_pdf(
  print(netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.heatmap = "Blues", title = "Incoming signaling")),
  "[正文]CellChat纯TME_incoming_signaling_role_heatmap.pdf",
  8,
  7
)

save_pdf(
  print(netAnalysis_signalingRole_scatter(cellchat)),
  "[补充]CellChat纯TME_signaling_role_scatter.pdf",
  6.5,
  5.8
)

# Pathway rankings.
pathway_rank <- read_csv(file.path(tables_dir, "CellChat纯TME_pathway无偏排序.csv"), show_col_types = FALSE) %>%
  mutate(rank = row_number(), is_PLAU = pathway_name == "PLAU")
write_csv(pathway_rank, file.path(source_data_dir, "CellChat纯TME_global_pathway_rank.csv"))

p_global_rank <- pathway_rank %>%
  slice_head(n = 30) %>%
  mutate(pathway_name = reorder(pathway_name, total_prob)) %>%
  ggplot(aes(x = pathway_name, y = total_prob, fill = is_PLAU)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#BC3C29"), guide = "none") +
  labs(x = NULL, y = "Total communication probability", title = "Global pathway ranking") +
  theme_classic(base_size = 10)
ggsave(file.path(figures_dir, "[正文]CellChat纯TME_rankNet全局_PLAU标注.pdf"), p_global_rank, width = 5.2, height = 6.2, useDingbats = FALSE)

source_target_rank <- read_csv(file.path(tables_dir, "CellChat纯TME_source_target_pathway无偏排序.csv"), show_col_types = FALSE)
mesv_outgoing_rank <- source_target_rank %>%
  filter(source == "Malignant_MES-V", !grepl("^Malignant_", target)) %>%
  group_by(pathway_name) %>%
  summarise(
    n_pairs = sum(n_pairs),
    total_prob = sum(total_prob, na.rm = TRUE),
    min_pval = min(min_pval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_prob), min_pval) %>%
  mutate(rank = row_number(), is_PLAU = pathway_name == "PLAU")
write_csv(mesv_outgoing_rank, file.path(source_data_dir, "CellChat纯TME_MESV_outgoing_pathway_rank.csv"))

p_mesv_rank <- mesv_outgoing_rank %>%
  slice_head(n = 30) %>%
  mutate(pathway_name = reorder(pathway_name, total_prob)) %>%
  ggplot(aes(x = pathway_name, y = total_prob, fill = is_PLAU)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#BC3C29"), guide = "none") +
  labs(x = NULL, y = "Total probability from MES-V to TME", title = "MES-V outgoing pathway ranking") +
  theme_classic(base_size = 10)
ggsave(file.path(figures_dir, "[正文]CellChat纯TME_rankNet_MESV_outgoing_PLAU标注.pdf"), p_mesv_rank, width = 5.2, height = 6.2, useDingbats = FALSE)

# PLAU pathway plots.
plau_lr <- extractEnrichedLR(cellchat, signaling = "PLAU", geneLR.return = FALSE, enriched.only = TRUE)
write_csv(as_tibble(plau_lr), file.path(source_data_dir, "CellChat纯TME_PLAU通路_enriched_LR.csv"))

save_pdf(
  netVisual_aggregate(
    cellchat,
    signaling = "PLAU",
    layout = "circle",
    vertex.weight = group_size,
    color.use = cell_colors[names(group_size)],
    weight.scale = TRUE,
    label.edge = FALSE
  ),
  "[正文]CellChat纯TME_PLAU通路圈图.pdf",
  7,
  7
)

save_pdf(
  netVisual_aggregate(
    cellchat,
    signaling = "PLAU",
    layout = "chord",
    color.use = cell_colors[names(group_size)]
  ),
  "[正文]CellChat纯TME_PLAU通路弦图.pdf",
  7,
  7
)

save_pdf(
  netVisual_heatmap(cellchat, signaling = "PLAU", color.heatmap = "Reds"),
  "[补充]CellChat纯TME_PLAU通路热图.pdf",
  6.2,
  5.5
)

save_pdf(
  netVisual_individual(
    cellchat,
    signaling = "PLAU",
    pairLR.use = plau_lr,
    layout = "circle",
    vertex.weight = group_size,
    color.use = cell_colors[names(group_size)],
    weight.scale = TRUE,
    label.edge = FALSE
  ),
  "[正文]CellChat纯TME_PLAU_PLAUR单LR圈图.pdf",
  7,
  7
)

save_pdf(
  print(netAnalysis_signalingRole_heatmap(cellchat, signaling = "PLAU", pattern = "outgoing", color.heatmap = "Reds", title = "PLAU outgoing")),
  "[正文]CellChat纯TME_PLAU_outgoing_role_heatmap.pdf",
  6.5,
  5.5
)

save_pdf(
  print(netAnalysis_signalingRole_heatmap(cellchat, signaling = "PLAU", pattern = "incoming", color.heatmap = "Blues", title = "PLAU incoming")),
  "[正文]CellChat纯TME_PLAU_incoming_role_heatmap.pdf",
  6.5,
  5.5
)

malignant_groups <- c("Malignant_NPC-P", "Malignant_OPC-M", "Malignant_MES-V", "Malignant_MES-I")
tme_groups <- setdiff(names(group_size), malignant_groups)

save_pdf(
  print(netVisual_bubble(
    cellchat,
    sources.use = malignant_groups,
    targets.use = tme_groups,
    signaling = "PLAU",
    remove.isolate = TRUE,
    angle.x = 45,
    title.name = "PLAU signaling from malignant subtypes to TME"
  )),
  "[正文]CellChat纯TME_PLAU_malignant_to_TME_bubble.pdf",
  7.5,
  4.8
)

save_pdf(
  print(netVisual_bubble(
    cellchat,
    sources.use = tme_groups,
    targets.use = malignant_groups,
    signaling = "PLAU",
    remove.isolate = TRUE,
    angle.x = 45,
    title.name = "PLAU signaling from TME to malignant subtypes"
  )),
  "[正文]CellChat纯TME_PLAU_TME_to_malignant_bubble.pdf",
  7.5,
  4.8
)

save_pdf(
  print(netVisual_bubble(
    cellchat,
    sources.use = malignant_groups,
    targets.use = tme_groups,
    remove.isolate = TRUE,
    angle.x = 45,
    title.name = "All significant malignant-to-TME interactions"
  )),
  "[补充]CellChat纯TME_全显著互作_malignant_to_TME_bubble.pdf",
  10,
  7.5
)

# Sensitivity: pureTME vs full. Absolute probabilities are not comparable.
full_rank <- read_csv(file.path(tables_dir, "CellChat_pathway无偏排序.csv"), show_col_types = FALSE) %>%
  mutate(rank = row_number(), analysis = "fullTME")
pure_rank <- read_csv(file.path(tables_dir, "CellChat纯TME_pathway无偏排序.csv"), show_col_types = FALSE) %>%
  mutate(rank = row_number(), analysis = "pureTME")

rank_compare <- bind_rows(full_rank, pure_rank) %>%
  filter(pathway_name %in% c("PLAU", "SPP1", "MIF", "COLLAGEN", "FN1")) %>%
  select(analysis, pathway_name, rank, n_pairs, total_prob, mean_prob, min_pval) %>%
  arrange(pathway_name, analysis)
write_csv(rank_compare, file.path(tables_dir, "CellChat_sensitivity_full_vs_pureTME_pathway_rank.csv"))

get_top_plau <- function(path, analysis_label) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(ligand == "PLAU", receptor == "PLAUR") %>%
    arrange(desc(prob), pval) %>%
    mutate(analysis = analysis_label, rank = row_number()) %>%
    select(analysis, rank, source, target, ligand, receptor, prob, pval, pathway_name) %>%
    slice_head(n = 20)
}

top_plau_compare <- bind_rows(
  get_top_plau(file.path(tables_dir, "CellChat_PLAUR相关通讯_按prob排序.csv"), "fullTME"),
  get_top_plau(file.path(tables_dir, "CellChat纯TME_PLAUR相关通讯_按prob排序.csv"), "pureTME")
)
write_csv(top_plau_compare, file.path(tables_dir, "CellChat_sensitivity_full_vs_pureTME_PLAU_top_edges.csv"))

writeLines(
  c(
    "CellChat pureTME figure and sensitivity outputs generated.",
    "Note: fullTME and pureTME CellChat probabilities are not directly comparable; compare ranking and dominant direction only."
  ),
  file.path(logs_dir, "06_CellChat纯TME出图与sensitivity_note.txt")
)

writeLines(capture.output(sessionInfo()), file.path(logs_dir, "06_CellChat纯TME出图与sensitivity_sessionInfo.txt"))
