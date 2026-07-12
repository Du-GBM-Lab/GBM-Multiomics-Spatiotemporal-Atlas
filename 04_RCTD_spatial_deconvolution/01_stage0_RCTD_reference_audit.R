# R9 RCTD preflight: reference candidate audit only.
# This script does not run RCTD or create deconvolution outputs.

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(Matrix)
})

root <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/R9_空间转录组"
dir.create(file.path(root, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(root, "logs"), showWarnings = FALSE, recursive = TRUE)

candidate_full <- c(
  immune_calls = "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/04_inferCNV_免疫参考验证/outputs/GBM.RNA.qc_doubletfinder.infercnv_immune_reference_calls.qs2",
  broad_calls = "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/03_inferCNV_恶性识别/outputs/GBM.RNA.qc_doubletfinder.infercnv_calls.qs2",
  qc_filtered = "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/02_scRNA_QC/outputs/GBM.RNA.qc_doubletfinder.filtered.qs2",
  cellchat_input = "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/07_细胞通讯/outputs/CellChat输入_四亚型加纯TME.qs2"
)
path_malignant <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/05_恶性细胞分亚群与Neftel对照/outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"

audit_one <- function(nm, path) {
  if (!file.exists(path)) {
    return(data.frame(candidate = nm, path = path, exists = FALSE))
  }
  obj <- tryCatch(qs2::qs_read(path), error = function(e) e)
  if (inherits(obj, "error")) {
    return(data.frame(candidate = nm, path = path, exists = TRUE, read_error = conditionMessage(obj)))
  }
  is_seurat <- inherits(obj, "Seurat")
  out <- data.frame(
    candidate = nm,
    path = path,
    exists = TRUE,
    class = paste(class(obj), collapse = "/"),
    is_seurat = is_seurat,
    n_features = if (is_seurat) nrow(obj) else NA_integer_,
    n_cells = if (is_seurat) ncol(obj) else NA_integer_,
    assays = if (is_seurat) paste(Assays(obj), collapse = ";") else NA_character_,
    default_assay = if (is_seurat) DefaultAssay(obj) else NA_character_,
    n_meta_cols = if (is_seurat) ncol(obj@meta.data) else NA_integer_,
    stringsAsFactors = FALSE
  )
  if (is_seurat) {
    meta_cols <- data.frame(candidate = nm, column = colnames(obj@meta.data))
    write.csv(meta_cols, file.path(root, "tables", paste0("stage0_RCTD_reference_meta_columns_", nm, ".csv")), row.names = FALSE)
    key_patterns <- "anno|annot|cell.?type|predicted|infercnv|malig|subtype|Subtype|comm_group|Pt|patient|sample|orig.ident|DF|doublet"
    hits <- grep(key_patterns, colnames(obj@meta.data), value = TRUE, ignore.case = TRUE)
    write.csv(data.frame(candidate = nm, column = hits), file.path(root, "tables", paste0("stage0_RCTD_reference_key_columns_", nm, ".csv")), row.names = FALSE)
  }
  rm(obj); gc()
  out
}

cat("== spacexr installed:", requireNamespace("spacexr", quietly = TRUE), "\n")
if (requireNamespace("spacexr", quietly = TRUE)) {
  ns <- asNamespace("spacexr")
  cat("== spacexr version:", as.character(utils::packageVersion("spacexr")), "\n")
  cat("== classic API create.RCTD/run.RCTD/Reference/SpatialRNA:",
      exists("create.RCTD", ns, inherits = FALSE),
      exists("run.RCTD", ns, inherits = FALSE),
      exists("Reference", ns, inherits = FALSE),
      exists("SpatialRNA", ns, inherits = FALSE), "\n")
  cat("== Bioc API createRctd/runRctd:",
      exists("createRctd", ns, inherits = FALSE),
      exists("runRctd", ns, inherits = FALSE), "\n")
}

candidate_audit <- do.call(rbind, Map(audit_one, names(candidate_full), candidate_full))
write.csv(candidate_audit, file.path(root, "tables", "stage0_RCTD_reference_candidate_objects.csv"), row.names = FALSE)
print(candidate_audit)

cat("\n== Load selected full candidate for label audit: immune_calls\n")
full <- qs2::qs_read(candidate_full[["immune_calls"]])
mal <- qs2::qs_read(path_malignant)
DefaultAssay(full) <- "RNA"
DefaultAssay(mal) <- "RNA"

cat("full dim:", nrow(full), ncol(full), "\n")
cat("mal dim:", nrow(mal), ncol(mal), "\n")

full_meta <- full@meta.data
mal_meta <- mal@meta.data

write.csv(data.frame(column = colnames(full_meta)), file.path(root, "tables", "stage0_RCTD_selected_full_meta_columns.csv"), row.names = FALSE)
write.csv(data.frame(column = colnames(mal_meta)), file.path(root, "tables", "stage0_RCTD_malignant_meta_columns.csv"), row.names = FALSE)

malignant_cols <- intersect(c("subtype_k4", "subtype_label_final", "Pt_number"), colnames(mal_meta))
stopifnot("subtype_k4" %in% malignant_cols)

full_barcodes <- colnames(full)
mal_barcodes <- colnames(mal)
mal_in_full <- intersect(full_barcodes, mal_barcodes)
cat("malignant barcodes in full:", length(mal_in_full), "\n")

candidate_nonmal_cols <- grep("anno|annot|cell.?type|predicted|comm_group|infercnv|malig", colnames(full_meta), value = TRUE, ignore.case = TRUE)
cat("candidate nonmalignant/type columns:", paste(candidate_nonmal_cols, collapse = ", "), "\n")

type_summaries <- list()
for (cc in candidate_nonmal_cols) {
  if (length(unique(full_meta[[cc]])) > 100) next
  tab <- sort(table(full_meta[[cc]], useNA = "ifany"), decreasing = TRUE)
  type_summaries[[cc]] <- data.frame(column = cc, label = names(tab), n = as.integer(tab))
}
type_summaries <- if (length(type_summaries)) do.call(rbind, type_summaries) else data.frame()
write.csv(type_summaries, file.path(root, "tables", "stage0_RCTD_selected_full_candidate_label_counts.csv"), row.names = FALSE)

# Build a provisional RCTD label audit using the most likely non-malignant column.
# Prefer broad human-readable annotation columns. This is only an audit, not a final object.
preferred_nonmal <- intersect(c("annotation", "anno_ident", "predicted.level3", "cell_type", "celltype", "comm_group"), colnames(full_meta))[1]
if (is.na(preferred_nonmal)) preferred_nonmal <- candidate_nonmal_cols[1]
cat("preferred nonmalignant label column for audit:", preferred_nonmal, "\n")

provisional <- as.character(full_meta[[preferred_nonmal]])
names(provisional) <- full_barcodes
subtype <- as.character(mal_meta$subtype_k4)
names(subtype) <- mal_barcodes
provisional[mal_in_full] <- subtype[mal_in_full]

tab <- sort(table(provisional, useNA = "ifany"), decreasing = TRUE)
write.csv(data.frame(RCTD_label = names(tab), n_cells = as.integer(tab)),
          file.path(root, "tables", "stage0_RCTD_provisional_label_counts.csv"),
          row.names = FALSE)

cnt <- tryCatch(GetAssayData(full, assay = "RNA", layer = "counts"),
                error = function(e1) tryCatch(GetAssayData(full, assay = "RNA", slot = "counts"), error = function(e2) NULL))
count_ok <- !is.null(cnt) && inherits(cnt, "dgCMatrix")
cat("RNA counts available:", count_ok, "\n")
if (count_ok) {
  nUMI <- Matrix::colSums(cnt)
  summary_df <- data.frame(
    metric = c("n_cells", "n_genes", "nUMI_min", "nUMI_median", "nUMI_mean", "nUMI_max"),
    value = c(ncol(cnt), nrow(cnt), min(nUMI), median(nUMI), mean(nUMI), max(nUMI))
  )
  write.csv(summary_df, file.path(root, "tables", "stage0_RCTD_reference_counts_summary.csv"), row.names = FALSE)
}

cat("\n[STOP RCTD reference audit]\n")
