suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(qs2)
})

module_dir <- normalizePath("07_发育时间_TF_ATAC验证", winslash = "/", mustWork = TRUE)
project_dir <- file.path(module_dir, "01_mouse_developmental")
source_dir <- "<DATA_ROOT>/项目/分型/分型代码/9.小鼠不同时期亚型变化"
data_dir <- file.path(project_dir, "data")
tables_dir <- file.path(module_dir, "tables")
outputs_dir <- file.path(project_dir, "outputs")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)

objects <- data.frame(
  label = c("tumourigenesis_atlas", "brain_injury_samples",
            "malignant_models", "annotated_injury"),
  path = file.path(source_dir, c(
    "GSE278511_tumourigenesis.atlas.rds",
    "GSE278511_brain.injury.samples.rds",
    "GSE278511_malignant.cells.Sox2CEPPT.and.NestinCPPT.models.rds",
    "注释好的小鼠单细胞injury.rds"
  )),
  stringsAsFactors = FALSE
)

qc_patterns <- c(
  "nCount", "nFeature", "percent", "mito", "mt", "ribo", "doublet",
  "sample", "stage", "cell_type", "annotation", "cluster", "seurat_clusters"
)

classify_stage <- function(x) {
  out <- rep(NA_character_, length(x))
  out[grepl("Preneoplastic", x, ignore.case = TRUE)] <- "Preneoplastic"
  out[grepl("Early_lesion", x, ignore.case = TRUE)] <- "Early lesion"
  out[grepl("Mid_lesion", x, ignore.case = TRUE)] <- "Mid lesion"
  out[grepl("Endpoint", x, ignore.case = TRUE)] <- "Endpoint"
  out
}

safe_dim <- function(expr) {
  out <- tryCatch(expr, error = function(e) NULL)
  if (is.null(out)) return(NA_character_)
  paste(dim(out), collapse = " x ")
}

safe_layer_dim <- function(obj, assay, layer) {
  assay_obj <- obj[[assay]]
  layers <- tryCatch(Layers(assay_obj), error = function(e) character())
  if (!(layer %in% layers)) return(NA_character_)
  safe_dim(LayerData(obj, assay = assay, layer = layer))
}

inspect_object <- function(obj, label, path) {
  meta <- obj[[]]
  assays <- Assays(obj)
  reductions <- Reductions(obj)
  default_assay <- DefaultAssay(obj)

  assay_summary <- do.call(rbind, lapply(assays, function(a) {
    DefaultAssay(obj) <- a
    data.frame(
      object = label,
      assay = a,
      layers = paste(tryCatch(Layers(obj[[a]]), error = function(e) character()), collapse = ";"),
      counts_dim = safe_layer_dim(obj, a, "counts"),
      data_dim = safe_layer_dim(obj, a, "data"),
      scale_data_dim = safe_layer_dim(obj, a, "scale.data"),
      stringsAsFactors = FALSE
    )
  }))
  DefaultAssay(obj) <- default_assay

  reduction_summary <- do.call(rbind, lapply(reductions, function(r) {
    emb <- Embeddings(obj, reduction = r)
    data.frame(
      object = label,
      reduction = r,
      dim = paste(dim(emb), collapse = " x "),
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(reduction_summary)) {
    reduction_summary <- data.frame(object = label, reduction = NA_character_, dim = NA_character_)
  }

  meta_fields <- colnames(meta)
  qc_fields <- meta_fields[grepl(paste(qc_patterns, collapse = "|"), meta_fields, ignore.case = TRUE)]
  field_summary <- data.frame(
    object = label,
    field = meta_fields,
    class = vapply(meta, function(x) paste(class(x), collapse = ";"), character(1)),
    n_unique = vapply(meta, function(x) length(unique(x)), integer(1)),
    n_na = vapply(meta, function(x) sum(is.na(x)), integer(1)),
    example_values = vapply(meta, function(x) {
      vals <- unique(as.character(x))
      vals <- vals[!is.na(vals)]
      paste(utils::head(vals, 8), collapse = " | ")
    }, character(1)),
    is_qc_candidate = meta_fields %in% qc_fields,
    stringsAsFactors = FALSE
  )

  sample_summary <- data.frame()
  if ("sample" %in% meta_fields) {
    sample_summary <- as.data.frame(table(sample = meta$sample), stringsAsFactors = FALSE)
    names(sample_summary)[2] <- "n_cells"
    sample_summary$object <- label
    sample_summary$stage_inferred <- classify_stage(as.character(sample_summary$sample))
    sample_summary <- sample_summary[, c("object", "sample", "stage_inferred", "n_cells")]
  }

  stage_summary <- data.frame()
  if (nrow(sample_summary) > 0) {
    sample_summary$stage_inferred[is.na(sample_summary$stage_inferred)] <- "unmatched"
    stage_summary <- aggregate(n_cells ~ object + stage_inferred, data = sample_summary, sum)
  }

  object_summary <- data.frame(
    object = label,
    source_path = path,
    file_size_gb = round(file.info(path)$size / 1024^3, 3),
    class = paste(class(obj), collapse = ";"),
    n_features = nrow(obj),
    n_cells = ncol(obj),
    n_metadata_fields = ncol(meta),
    default_assay = default_assay,
    assays = paste(assays, collapse = ";"),
    reductions = paste(reductions, collapse = ";"),
    has_sample = "sample" %in% meta_fields,
    qc_candidate_fields = paste(qc_fields, collapse = ";"),
    stringsAsFactors = FALSE
  )

  list(
    object_summary = object_summary,
    assay_summary = assay_summary,
    reduction_summary = reduction_summary,
    field_summary = field_summary,
    sample_summary = sample_summary,
    stage_summary = stage_summary
  )
}

all_results <- list()
for (i in seq_len(nrow(objects))) {
  label <- objects$label[i]
  path <- objects$path[i]
  message("Reading: ", label, " <- ", path)
  obj <- readRDS(path)
  all_results[[label]] <- inspect_object(obj, label, path)
  rm(obj)
  invisible(gc())
}

bind_part <- function(part) {
  do.call(rbind, lapply(all_results, `[[`, part))
}

write.csv(bind_part("object_summary"),
          file.path(tables_dir, "小鼠对象总体审计.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(bind_part("assay_summary"),
          file.path(tables_dir, "小鼠对象assay审计.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(bind_part("reduction_summary"),
          file.path(tables_dir, "小鼠对象reduction审计.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(bind_part("field_summary"),
          file.path(tables_dir, "小鼠对象metadata字段审计.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(bind_part("sample_summary"),
          file.path(tables_dir, "小鼠对象sample阶段审计.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(bind_part("stage_summary"),
          file.path(tables_dir, "小鼠对象stage细胞数审计.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

save_clean_seurat <- function(source_path, out_name) {
  obj <- readRDS(source_path)
  meta <- obj[[]]
  meta$stage_inferred <- if ("sample" %in% colnames(meta)) {
    classify_stage(as.character(meta$sample))
  } else {
    NA_character_
  }
  obj <- AddMetaData(obj, metadata = meta$stage_inferred, col.name = "stage_inferred")

  keep_reductions <- intersect(c("umap", "pca", "harmony"), Reductions(obj))
  clean <- DietSeurat(
    obj,
    assays = intersect("RNA", Assays(obj)),
    dimreducs = keep_reductions,
    graphs = NULL,
    misc = FALSE,
    layers = c("counts", "data")
  )
  clean@commands <- list()
  if ("tools" %in% slotNames(clean)) clean@tools <- list()

  out_path <- file.path(data_dir, out_name)
  qs2::qs_save(clean, out_path)
  rm(obj, clean)
  invisible(gc())
  out_path
}

message("Creating cleaned copies")
cleaned_atlas <- save_clean_seurat(
  objects$path[objects$label == "tumourigenesis_atlas"],
  "GSE278511小鼠发育时间_atlas精简.qs2"
)
cleaned_malignant <- save_clean_seurat(
  objects$path[objects$label == "malignant_models"],
  "GSE278511小鼠恶性模型_精简.qs2"
)
cleaned_annotated <- save_clean_seurat(
  objects$path[objects$label == "annotated_injury"],
  "GSE278511小鼠injury注释对象_精简.qs2"
)

cnv_ref_dir <- file.path(data_dir, "cnv_ref")
dir.create(cnv_ref_dir, recursive = TRUE, showWarnings = FALSE)
gene_order_src <- file.path(
  source_dir,
  "cnv_ref",
  "mouse_gencode.GRCm39.vM32.basic.annotation.by_gene_name.infercnv_positions.txt"
)
gene_order_dst <- file.path(cnv_ref_dir, basename(gene_order_src))
invisible(file.copy(gene_order_src, gene_order_dst, overwrite = TRUE))

manifest <- data.frame(
  item = c(
    "source_dir",
    "source_object_tumourigenesis_atlas",
    "cleaned_object_tumourigenesis_atlas",
    "cleaned_object_malignant_models",
    "cleaned_object_annotated_injury",
    "gene_order_file"
  ),
  path = c(
    source_dir,
    objects$path[objects$label == "tumourigenesis_atlas"],
    cleaned_atlas,
    cleaned_malignant,
    cleaned_annotated,
    gene_order_dst
  ),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(tables_dir, "小鼠数据文件清单.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

message("Done.")
