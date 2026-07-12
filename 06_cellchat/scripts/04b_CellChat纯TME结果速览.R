#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

module_dir <- "07_细胞通讯"
tables_dir <- file.path(module_dir, "tables")
logs_dir <- file.path(module_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

comm_path <- file.path(tables_dir, "CellChat纯TME_全部通讯结果.csv")
stopifnot(file.exists(comm_path))

comm_all <- read_csv(comm_path, show_col_types = FALSE)

plaur_ranked <- comm_all %>%
  filter(ligand %in% c("PLAU", "PLAUR") | receptor %in% c("PLAU", "PLAUR")) %>%
  arrange(desc(prob), pval, source, target)
write_csv(plaur_ranked, file.path(tables_dir, "CellChat纯TME_PLAUR相关通讯_按prob排序.csv"))

malignant_groups <- c("Malignant_NPC-P", "Malignant_OPC-M", "Malignant_MES-V", "Malignant_MES-I")

plaur_malignant_tme <- plaur_ranked %>%
  filter(
    (source %in% malignant_groups & !target %in% malignant_groups) |
      (!source %in% malignant_groups & target %in% malignant_groups)
  ) %>%
  mutate(direction_class = case_when(
    source %in% malignant_groups ~ "malignant_to_TME",
    target %in% malignant_groups ~ "TME_to_malignant",
    TRUE ~ "other"
  )) %>%
  arrange(desc(prob), direction_class, source, target)
write_csv(plaur_malignant_tme, file.path(tables_dir, "CellChat纯TME_PLAUR恶性亚型_TME通讯.csv"))

plaur_pair_summary <- plaur_malignant_tme %>%
  group_by(direction_class, source, target, interaction_name, pathway_name) %>%
  summarise(
    n_pairs = n(),
    total_prob = sum(prob, na.rm = TRUE),
    max_prob = max(prob, na.rm = TRUE),
    min_pval = min(pval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(direction_class, desc(total_prob), min_pval)
write_csv(plaur_pair_summary, file.path(tables_dir, "CellChat纯TME_PLAUR恶性亚型_TME通讯汇总.csv"))

pathway_rank <- read_csv(file.path(tables_dir, "CellChat纯TME_pathway无偏排序.csv"), show_col_types = FALSE) %>%
  mutate(rank = row_number())
pathway_focus <- pathway_rank %>%
  filter(pathway_name %in% c("PLAU", "VTN", "PARs", "THBS", "SPP1", "MIF", "FN1", "COLLAGEN"))
write_csv(pathway_focus, file.path(tables_dir, "CellChat纯TME_重点pathway无偏排名.csv"))

writeLines(
  capture.output(sessionInfo()),
  file.path(logs_dir, "04b_CellChat纯TME结果速览_sessionInfo.txt")
)
