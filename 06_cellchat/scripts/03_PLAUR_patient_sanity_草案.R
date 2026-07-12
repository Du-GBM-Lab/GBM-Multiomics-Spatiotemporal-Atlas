#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(Matrix)
})

module_dir <- "07_细胞通讯"
outputs_dir <- file.path(module_dir, "outputs")
tables_dir <- file.path(module_dir, "tables")
logs_dir <- file.path(module_dir, "logs")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

obj_path <- file.path(outputs_dir, "CellChat输入_四亚型加TME.qs2")
stopifnot(file.exists(obj_path))

obj <- qs2::qs_read(obj_path)
DefaultAssay(obj) <- "RNA"

genes <- intersect(c("PLAU", "PLAUR", "VTN"), rownames(obj))
if (length(genes) == 0) stop("None of PLAU/PLAUR/VTN found in object.")

expr_mat <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(obj, assay = "RNA", layer = "counts")
)
expr <- expr_mat[genes, , drop = FALSE]
if (nrow(expr) == 0 || ncol(expr) == 0) stop("Empty expression matrix.")

md <- obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  select(cell, Pt_number, comm_group)

expr_long <- lapply(genes, function(g) {
  x <- as.numeric(expr[g, ])
  tibble(
    cell = colnames(obj),
    gene = g,
    expr = x,
    detected = x > 0
  )
}) %>%
  bind_rows() %>%
  left_join(md, by = "cell")

gene_group_patient <- expr_long %>%
  group_by(gene, comm_group, Pt_number) %>%
  summarise(
    n_cells = n(),
    pct_detected = mean(detected) * 100,
    mean_expr = mean(expr),
    .groups = "drop"
  ) %>%
  arrange(gene, comm_group, desc(n_cells))
write_csv(gene_group_patient, file.path(tables_dir, "PLAU_PLAUR_patient_x_comm_group表达审计.csv"))

gene_group_summary <- expr_long %>%
  group_by(gene, comm_group) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(Pt_number),
    pct_detected = mean(detected) * 100,
    mean_expr = mean(expr),
    .groups = "drop"
  ) %>%
  arrange(gene, desc(mean_expr), desc(pct_detected))
write_csv(gene_group_summary, file.path(tables_dir, "PLAU_PLAUR_comm_group表达摘要.csv"))

patient_dominance <- gene_group_patient %>%
  group_by(gene, comm_group) %>%
  mutate(total_detected_weight = sum(pct_detected * n_cells, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(detected_weight = pct_detected * n_cells) %>%
  group_by(gene, comm_group) %>%
  summarise(
    n_patients = n_distinct(Pt_number),
    top_patient = Pt_number[which.max(detected_weight)],
    top_patient_detected_weight_fraction = max(detected_weight, na.rm = TRUE) / sum(detected_weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(gene, desc(top_patient_detected_weight_fraction))
write_csv(patient_dominance, file.path(tables_dir, "PLAU_PLAUR_patient主导性审计.csv"))

writeLines(
  capture.output(sessionInfo()),
  file.path(logs_dir, "03_PLAUR_patient_sanity_sessionInfo.txt")
)

message("PLAU/PLAUR patient sanity tables written to ", tables_dir)
