#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(readr)
  library(tibble)
})

module_dir <- "07_细胞通讯"
tables_dir <- file.path(module_dir, "tables")
logs_dir <- file.path(module_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

full_obj_path <- "04_inferCNV_免疫参考验证/outputs/GBM.RNA.qc_doubletfinder.infercnv_immune_reference_calls.qs2"
malignant_obj_path <- "05_恶性细胞分亚群与Neftel对照/outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"

stopifnot(file.exists(full_obj_path), file.exists(malignant_obj_path))

full_obj <- qs2::qs_read(full_obj_path)
mal_obj <- qs2::qs_read(malignant_obj_path)

required_full <- c("anno_ident", "malignant_call_status_immune", "is_malignant_for_downstream_immune", "Pt_number")
required_mal <- c("subtype_label_final")

missing_full <- setdiff(required_full, colnames(full_obj@meta.data))
missing_mal <- setdiff(required_mal, colnames(mal_obj@meta.data))
if (length(missing_full) > 0) stop("Full object missing metadata: ", paste(missing_full, collapse = ", "))
if (length(missing_mal) > 0) stop("Malignant object missing metadata: ", paste(missing_mal, collapse = ", "))

barcode_audit <- tibble(
  metric = c(
    "full_cells",
    "malignant_cells",
    "barcode_overlap",
    "malignant_all_in_full"
  ),
  value = c(
    ncol(full_obj),
    ncol(mal_obj),
    length(intersect(colnames(full_obj), colnames(mal_obj))),
    as.character(all(colnames(mal_obj) %in% colnames(full_obj)))
  )
)
write_csv(barcode_audit, file.path(tables_dir, "Phase0_barcode对齐审计.csv"))

full_md <- full_obj@meta.data %>%
  as_tibble(rownames = "cell")

mal_md <- mal_obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  select(cell, subtype_label_final)

anno_counts <- full_md %>%
  count(anno_ident, malignant_call_status_immune, name = "n_cells") %>%
  arrange(desc(n_cells))
write_csv(anno_counts, file.path(tables_dir, "Phase0_完整对象注释与恶性状态计数.csv"))

subtype_counts <- mal_md %>%
  count(subtype_label_final, name = "n_cells") %>%
  arrange(desc(n_cells))
write_csv(subtype_counts, file.path(tables_dir, "Phase0_当前四亚型计数.csv"))

ambiguous_audit <- full_md %>%
  mutate(
    ambiguity_class = case_when(
      malignant_call_status_immune == "ambiguous_burden_only" ~ "ambiguous_burden_only",
      anno_ident == "Ambiguous" ~ "annotation_Ambiguous",
      TRUE ~ "main_analysis_candidate"
    )
  ) %>%
  count(ambiguity_class, malignant_call_status_immune, anno_ident, name = "n_cells") %>%
  arrange(ambiguity_class, desc(n_cells))
write_csv(ambiguous_audit, file.path(tables_dir, "Phase0_ambiguous细胞审计.csv"))

gene_check <- c("PLAU", "PLAUR", "MMP9", "MMP14", "SERPINE1", "VIM", "CD44")
gene_presence <- tibble(
  gene = gene_check,
  present = gene_check %in% rownames(full_obj)
)
write_csv(gene_presence, file.path(tables_dir, "Phase0_PLAUR相关基因存在性.csv"))

if (requireNamespace("CellChat", quietly = TRUE)) {
  data("CellChatDB.human", package = "CellChat")
  db_interaction <- CellChatDB.human$interaction %>%
    as_tibble()

  plaur_db_hits <- db_interaction %>%
    filter(
      ligand %in% c("PLAU", "PLAUR") |
        receptor %in% c("PLAU", "PLAUR") |
        grepl("\\bPLAU\\b|\\bPLAUR\\b", interaction_name, ignore.case = FALSE) |
        grepl("\\bPLAU\\b|\\bPLAUR\\b", pathway_name, ignore.case = FALSE)
    ) %>%
    arrange(pathway_name, interaction_name)

  write_csv(plaur_db_hits, file.path(tables_dir, "Phase0_CellChatDB_PLAUR相关互作.csv"))

  db_annotation_cols <- intersect(
    c("annotation", "evidence", "interaction_name", "pathway_name", "ligand", "receptor"),
    colnames(plaur_db_hits)
  )
  plaur_db_summary <- plaur_db_hits %>%
    select(all_of(db_annotation_cols))
  write_csv(plaur_db_summary, file.path(tables_dir, "Phase0_CellChatDB_PLAUR互作类别摘要.csv"))
} else {
  write_csv(
    tibble(
      status = "CellChat_not_installed",
      note = "Install/load CellChat before confirming PLAU-PLAUR interaction category."
    ),
    file.path(tables_dir, "Phase0_CellChatDB_PLAUR相关互作.csv")
  )
}

writeLines(
  capture.output(sessionInfo()),
  file.path(logs_dir, "00_Phase0_对象与metadata审计_sessionInfo.txt")
)

message("Phase 0 audit written to ", tables_dir)
