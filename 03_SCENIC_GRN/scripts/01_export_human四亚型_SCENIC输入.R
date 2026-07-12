#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(Seurat)
  library(SeuratObject)
  library(qs2)
  library(Matrix)
})

module_dir <- normalizePath("07_发育时间_TF_ATAC验证", winslash = "/", mustWork = TRUE)
project_root <- dirname(module_dir)
project_dir <- file.path(module_dir, "02_TF_regulon_SCENIC")
data_dir <- file.path(project_dir, "data")
logs_dir <- file.path(project_dir, "logs")
tables_dir <- file.path(module_dir, "tables")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

human_obj_path <- file.path(
  project_root,
  "05_恶性细胞分亚群与Neftel对照/outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"
)
stopifnot(file.exists(human_obj_path))

short_label_map <- c(
  "Proliferative-NPC" = "NPC-P",
  "OPC-Myelination" = "OPC-M",
  "Vascular-niche MES" = "MES-V",
  "MES-Antigen-presenting" = "MES-I",
  "Proliferative-NPC subtype" = "NPC-P",
  "OPC-Myelination subtype" = "OPC-M",
  "Vascular-niche MES subtype" = "MES-V",
  "MES-Antigen-presenting subtype" = "MES-I"
)
subtype_levels <- c("NPC-P", "OPC-M", "MES-V", "MES-I")

obj <- qs2::qs_read(human_obj_path)

required_meta <- c("subtype_k4", "subtype_label_final")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Human object missing metadata: ", paste(missing_meta, collapse = ", "))
}

obj$cell_type <- unname(short_label_map[obj$subtype_label_final])
if (any(is.na(obj$cell_type))) {
  bad <- sort(unique(obj$subtype_label_final[is.na(obj$cell_type)]))
  stop("Unmapped subtype_label_final values: ", paste(bad, collapse = ", "))
}
if (!setequal(unique(obj$cell_type), subtype_levels)) {
  stop("cell_type does not resolve to NPC-P / OPC-M / MES-V / MES-I.")
}
obj$cell_type <- factor(obj$cell_type, levels = subtype_levels)

DefaultAssay(obj) <- "RNA"
counts <- if (packageVersion("SeuratObject") >= "5.0.0") {
  SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
} else {
  Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")
}
counts <- counts[, colnames(obj), drop = FALSE]

Matrix::writeMM(counts, file.path(data_dir, "gbm_counts.mtx"))
write.csv(
  data.frame(gene = rownames(counts), check.names = FALSE),
  file.path(data_dir, "gbm_genes.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  data.frame(cell = colnames(counts), check.names = FALSE),
  file.path(data_dir, "gbm_barcodes.csv"),
  row.names = FALSE,
  quote = TRUE
)

metadata <- obj@meta.data
metadata$cell <- rownames(metadata)
metadata <- metadata[colnames(counts), , drop = FALSE]
metadata$cell_type <- as.character(obj$cell_type)
metadata$subtype_short <- as.character(obj$cell_type)
write.csv(metadata, file.path(data_dir, "gbm_metadata.csv"), row.names = TRUE, quote = TRUE)

reduction_name <- intersect(
  c("umap", "umap_compact_balanced", "umap_baseline", "umap_tight"),
  names(obj@reductions)
)[1]
if (is.na(reduction_name)) {
  stop("No usable UMAP reduction found.")
}
umap <- Embeddings(obj, reduction = reduction_name)
umap <- umap[colnames(counts), , drop = FALSE]
write.csv(umap, file.path(data_dir, "gbm_umap.csv"), row.names = TRUE, quote = TRUE)

audit <- data.frame(
  subtype = names(table(obj$cell_type)),
  n_cells = as.integer(table(obj$cell_type)),
  stringsAsFactors = FALSE
)
write.csv(audit, file.path(tables_dir, "scenic_phase1_human四亚型输入审计.csv"), row.names = FALSE)

log_lines <- c(
  paste0("human_obj_path=", normalizePath(human_obj_path, winslash = "/")),
  paste0("n_cells=", ncol(obj)),
  paste0("n_features=", nrow(obj)),
  paste0("counts_dim=", paste(dim(counts), collapse = " x ")),
  paste0("umap_reduction=", reduction_name),
  paste0("subtype_counts=", paste(audit$subtype, audit$n_cells, sep = ":", collapse = "; "))
)
writeLines(log_lines, file.path(logs_dir, "01_export_SCENIC输入日志.txt"))
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "01_export_SCENIC输入_sessionInfo.txt"))

