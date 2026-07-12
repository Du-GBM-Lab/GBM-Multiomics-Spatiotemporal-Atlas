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

pure_tme_groups <- c(
  "Macrophages", "Microglial", "Monocytes", "cDCs", "pDCs",
  "T cells", "NK cells", "B cells", "Endothelial", "Mural cells"
)

excluded_glial_like_groups <- c(
  "Radial glial", "OPCs", "Oligodendrocytes", "Astrocytes",
  "Neurons", "Ependymal cells"
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
      malignant_call_status_immune == "non_malignant" & anno_ident %in% pure_tme_groups ~
        as.character(anno_ident),
      TRUE ~ NA_character_
    ),
    exclusion_reason = case_when(
      malignant_call_status_immune == "ambiguous_burden_only" ~ "ambiguous_burden_only",
      anno_ident == "Ambiguous" ~ "annotation_Ambiguous",
      malignant_call_status_immune == "non_malignant" & anno_ident %in% excluded_glial_like_groups ~
        "excluded_glial_neural_like_TME",
      malignant_call_status_immune == "non_malignant" & !anno_ident %in% pure_tme_groups ~
        "excluded_other_non_malignant",
      TRUE ~ NA_character_
    ),
    comm_keep_pure_tme = !is.na(comm_group)
  )

full_obj$malignant_subtype <- full_md$malignant_subtype[match(colnames(full_obj), full_md$cell)]
full_obj$comm_group <- full_md$comm_group[match(colnames(full_obj), full_md$cell)]
full_obj$comm_keep_pure_tme <- full_md$comm_keep_pure_tme[match(colnames(full_obj), full_md$cell)]
full_obj$comm_exclusion_reason_pure_tme <- full_md$exclusion_reason[match(colnames(full_obj), full_md$cell)]

group_counts <- full_obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  count(comm_group, comm_keep_pure_tme, name = "n_cells") %>%
  arrange(desc(n_cells))
write_csv(group_counts, file.path(tables_dir, "CellChat输入_纯TME_comm_group计数.csv"))

patient_group_counts <- full_obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  mutate(Pt_number = as.character(Pt_number)) %>%
  count(Pt_number, comm_group, comm_keep_pure_tme, name = "n_cells") %>%
  arrange(Pt_number, desc(n_cells))
write_csv(patient_group_counts, file.path(tables_dir, "CellChat输入_纯TME_patient_x_comm_group计数.csv"))

excluded_counts <- full_obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  filter(!comm_keep_pure_tme) %>%
  count(comm_exclusion_reason_pure_tme, malignant_call_status_immune, anno_ident, name = "n_cells") %>%
  arrange(comm_exclusion_reason_pure_tme, desc(n_cells))
write_csv(excluded_counts, file.path(tables_dir, "CellChat输入_纯TME_排除细胞审计.csv"))

cellchat_obj <- subset(full_obj, subset = comm_keep_pure_tme)
cellchat_obj$comm_group <- factor(cellchat_obj$comm_group)

qs2::qs_save(
  cellchat_obj,
  file.path(outputs_dir, "CellChat输入_四亚型加纯TME.qs2")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(logs_dir, "01b_构建CellChat输入_纯TME_sessionInfo.txt")
)

message("Pure TME CellChat input object written to ", outputs_dir)
