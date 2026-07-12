#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(readr)
  library(tibble)
})

module_dir <- "07_细胞通讯"
outputs_dir <- file.path(module_dir, "outputs")
tables_dir <- file.path(module_dir, "tables")
logs_dir <- file.path(module_dir, "logs")
dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

full_obj_path <- "04_inferCNV_免疫参考验证/outputs/GBM.RNA.qc_doubletfinder.infercnv_immune_reference_calls.qs2"
malignant_obj_path <- "05_恶性细胞分亚群与Neftel对照/outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"

short_label_map <- c(
  "Proliferative-NPC" = "NPC-P",
  "OPC-Myelination" = "OPC-M",
  "Vascular-niche MES" = "MES-V",
  "MES-Antigen-presenting" = "MES-I"
)

full_obj <- qs2::qs_read(full_obj_path)
mal_obj <- qs2::qs_read(malignant_obj_path)

mal_subtype <- mal_obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    cell,
    malignant_subtype = unname(short_label_map[subtype_label_final])
  )

if (any(is.na(mal_subtype$malignant_subtype))) {
  stop("Unmapped malignant subtype labels found.")
}

full_md <- full_obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  left_join(mal_subtype, by = "cell") %>%
  mutate(
    comm_group = case_when(
      malignant_call_status_immune == "malignant_high_confidence" & !is.na(malignant_subtype) ~
        paste0("Malignant_", malignant_subtype),
      malignant_call_status_immune == "non_malignant" ~ as.character(anno_ident),
      TRUE ~ NA_character_
    ),
    comm_keep = !is.na(comm_group) & comm_group != "Ambiguous"
  )

full_obj$malignant_subtype <- full_md$malignant_subtype[match(colnames(full_obj), full_md$cell)]
full_obj$comm_group <- full_md$comm_group[match(colnames(full_obj), full_md$cell)]
full_obj$comm_keep <- full_md$comm_keep[match(colnames(full_obj), full_md$cell)]

group_counts <- full_obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  count(comm_group, comm_keep, name = "n_cells") %>%
  arrange(desc(n_cells))
write_csv(group_counts, file.path(tables_dir, "CellChat输入_comm_group计数.csv"))

patient_group_counts <- full_obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  mutate(Pt_number = as.character(Pt_number)) %>%
  count(Pt_number, comm_group, comm_keep, name = "n_cells") %>%
  arrange(Pt_number, desc(n_cells))
write_csv(patient_group_counts, file.path(tables_dir, "CellChat输入_patient_x_comm_group计数.csv"))

excluded_counts <- full_obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  filter(!comm_keep) %>%
  count(malignant_call_status_immune, anno_ident, comm_group, name = "n_cells") %>%
  arrange(desc(n_cells))
write_csv(excluded_counts, file.path(tables_dir, "CellChat输入_排除细胞审计.csv"))

cellchat_obj <- subset(full_obj, subset = comm_keep)
cellchat_obj$comm_group <- factor(cellchat_obj$comm_group)

qs2::qs_save(
  cellchat_obj,
  file.path(outputs_dir, "CellChat输入_四亚型加TME.qs2")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(logs_dir, "01_构建CellChat输入_sessionInfo.txt")
)

message("CellChat input object written to ", outputs_dir)
