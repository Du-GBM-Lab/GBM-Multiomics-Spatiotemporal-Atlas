root <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(Matrix)
})

path_full <- file.path(root, "data/reference/GBM.RNA.qc_doubletfinder.infercnv_immune_reference_calls.qs2")
path_mal <- file.path(root, "data/reference/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2")

full <- qs2::qs_read(path_full)
mal <- qs2::qs_read(path_mal)
DefaultAssay(full) <- "RNA"
DefaultAssay(mal) <- "RNA"

cat("full:", nrow(full), ncol(full), "assays=", paste(Assays(full), collapse = ","), "\n")
cat("mal:", nrow(mal), ncol(mal), "assays=", paste(Assays(mal), collapse = ","), "\n")

ov <- intersect(colnames(mal), colnames(full))
cat("malignant overlap:", length(ov), "/", ncol(mal), "\n")

ct <- as.character(full$anno_ident)
names(ct) <- colnames(full)
ct[ov] <- as.character(mal$subtype_k4[match(ov, colnames(mal))])
tab <- sort(table(ct), decreasing = TRUE)

print(tab)
cat("all_gt_50:", all(tab > 50), "\n")

dir.create(file.path(root, "tables"), showWarnings = FALSE, recursive = TRUE)
write.csv(
  data.frame(celltype = names(tab), n_cells = as.integer(tab)),
  file.path(root, "tables/server_reference_only_celltype_counts.csv"),
  row.names = FALSE
)

cat("[STOP reference only sanity]\n")
