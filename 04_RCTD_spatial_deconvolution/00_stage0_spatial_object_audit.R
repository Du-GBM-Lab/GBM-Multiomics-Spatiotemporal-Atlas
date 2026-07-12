# R9 spatial transcriptomics Stage 0 object audit
# Read-only audit of the legacy merged spatial object. No deconvolution,
# correlation, mapping, or plotting is performed here.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

root <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/R9_空间转录组"
dir.create(file.path(root, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(root, "docs"), showWarnings = FALSE, recursive = TRUE)

path_st <- "<DATA_ROOT>/项目/分型/分型代码/0.对象/5.ST_merge.rds"
stopifnot(file.exists(path_st))

cat("== Input:", path_st, "\n")
cat("== File size GB:", round(file.info(path_st)$size / 1024^3, 3), "\n")

obj <- readRDS(path_st)

cat("== Class:", paste(class(obj), collapse = ", "), "\n")
cat("== Dim:", nrow(obj), "features x", ncol(obj), "spots/cells\n")
cat("== Assays:", paste(Assays(obj), collapse = ", "), "\n")
cat("== DefaultAssay:", DefaultAssay(obj), "\n")
cat("== Reductions:", paste(Reductions(obj), collapse = ", "), "\n")
cat("== Images:", paste(names(obj@images), collapse = ", "), "\n")
cat("== Metadata columns:", ncol(obj@meta.data), "\n")

meta <- obj@meta.data
sample_col <- intersect(c("orig.ident", "sample", "Sample", "slice", "slice_id"), colnames(meta))[1]
if (is.na(sample_col)) sample_col <- NA_character_

if (!is.na(sample_col)) {
  sample_tab <- sort(table(meta[[sample_col]], useNA = "ifany"), decreasing = TRUE)
  sample_df <- data.frame(sample = names(sample_tab), n_spots = as.integer(sample_tab))
  write.csv(sample_df, file.path(root, "tables", "stage0_spatial_sample_spot_counts.csv"), row.names = FALSE)
  cat("\n== Sample column:", sample_col, "\n")
  cat("== N samples/slices:", length(sample_tab), "\n")
  print(sample_tab)
} else {
  cat("\n== Sample column: NOT FOUND\n")
}

assay_rows <- do.call(rbind, lapply(Assays(obj), function(a) {
  layers <- tryCatch(Layers(obj[[a]]), error = function(e) character(0))
  mat <- tryCatch(
    GetAssayData(obj, assay = a, layer = "counts"),
    error = function(e1) tryCatch(
      GetAssayData(obj, assay = a, slot = "counts"),
      error = function(e2) NULL
    )
  )
  data.frame(
    assay = a,
    assay_class = class(obj[[a]])[1],
    layers = if (length(layers)) paste(layers, collapse = ";") else NA_character_,
    n_features = if (is.null(mat)) NA_integer_ else nrow(mat),
    n_cells = if (is.null(mat)) NA_integer_ else ncol(mat),
    counts_class = if (is.null(mat)) NA_character_ else class(mat)[1],
    stringsAsFactors = FALSE
  )
}))
write.csv(assay_rows, file.path(root, "tables", "stage0_spatial_assays.csv"), row.names = FALSE)

spatial_layers <- tryCatch(Layers(obj[["Spatial"]]), error = function(e) character(0))
if (length(spatial_layers)) {
  layer_dims <- do.call(rbind, lapply(spatial_layers, function(ly) {
    mat <- LayerData(obj[["Spatial"]], layer = ly)
    data.frame(layer = ly, n_features = nrow(mat), n_spots = ncol(mat), class = class(mat)[1])
  }))
  write.csv(layer_dims, file.path(root, "tables", "stage0_spatial_layers.csv"), row.names = FALSE)
}

patterns <- list(
  card_cols = "^CARD_",
  old_subtype_score_cols = "^Subtype_Score",
  current_subtype_cols = "subtype|Subtype|MES|NPC|OPC|label",
  location_cols = "^Location$|^PAN$|^MVP$|^LE$|^CT$|niche|Ivy|IVY",
  spatial_coord_cols = "row|col|image|imagerow|imagecol|x|y"
)

col_hits <- lapply(names(patterns), function(nm) {
  hits <- grep(patterns[[nm]], colnames(meta), value = TRUE, ignore.case = TRUE)
  if (!length(hits)) hits <- NA_character_
  data.frame(category = nm, column = hits, stringsAsFactors = FALSE)
})
col_hits <- do.call(rbind, col_hits)
write.csv(col_hits, file.path(root, "tables", "stage0_spatial_metadata_column_hits.csv"), row.names = FALSE)
write.csv(data.frame(column = colnames(meta)), file.path(root, "tables", "stage0_spatial_metadata_columns.csv"), row.names = FALSE)

cat("\n== Metadata column hits:\n")
cat("== All metadata columns:", paste(colnames(meta), collapse = ", "), "\n")
for (nm in names(patterns)) {
  hits <- grep(patterns[[nm]], colnames(meta), value = TRUE, ignore.case = TRUE)
  cat("  ", nm, ":", if (length(hits)) paste(hits, collapse = ", ") else "NONE", "\n")
}

if (any(grepl("^CARD_", colnames(meta)))) {
  card_cols <- grep("^CARD_", colnames(meta), value = TRUE)
  card_summary <- data.frame(
    column = card_cols,
    min = sapply(card_cols, function(x) min(meta[[x]], na.rm = TRUE)),
    median = sapply(card_cols, function(x) median(meta[[x]], na.rm = TRUE)),
    mean = sapply(card_cols, function(x) mean(meta[[x]], na.rm = TRUE)),
    max = sapply(card_cols, function(x) max(meta[[x]], na.rm = TRUE)),
    nonzero_rate = sapply(card_cols, function(x) mean(meta[[x]] > 0, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  write.csv(card_summary, file.path(root, "tables", "stage0_spatial_CARD_column_summary.csv"), row.names = FALSE)
  cat("\n== CARD columns:", length(card_cols), "\n")
  print(card_summary)
}

if ("Location" %in% colnames(meta)) {
  loc_tab <- if (!is.na(sample_col)) {
    as.data.frame.matrix(table(meta[[sample_col]], meta$Location, useNA = "ifany"))
  } else {
    as.data.frame(table(meta$Location, useNA = "ifany"))
  }
  write.csv(loc_tab, file.path(root, "tables", "stage0_spatial_location_table.csv"))
  cat("\n== Location table:\n")
  print(table(meta$Location, useNA = "ifany"))
}

genes_interest <- c("FOSL1", "PLAU", "PLAUR", "POSTN", "CDH2", "ANXA2")
gene_rows <- do.call(rbind, lapply(Assays(obj), function(a) {
  feats <- rownames(obj[[a]])
  data.frame(assay = a, gene = genes_interest, present = genes_interest %in% feats)
}))
write.csv(gene_rows, file.path(root, "tables", "stage0_spatial_gene_presence.csv"), row.names = FALSE)
cat("\n== Gene presence:\n")
print(gene_rows)

image_df <- data.frame(
  image = names(obj@images),
  class = sapply(obj@images, function(x) paste(class(x), collapse = "/")),
  stringsAsFactors = FALSE
)
write.csv(image_df, file.path(root, "tables", "stage0_spatial_images.csv"), row.names = FALSE)
cat("\n== Image classes:\n")
print(image_df)

summary_lines <- c(
  "# R9 Stage 0 spatial object audit",
  "",
  paste0("- input: `", path_st, "`"),
  paste0("- class: ", paste(class(obj), collapse = ", ")),
  paste0("- dimensions: ", nrow(obj), " features x ", ncol(obj), " spots/cells"),
  paste0("- assays: ", paste(Assays(obj), collapse = ", ")),
  paste0("- default assay: ", DefaultAssay(obj)),
  paste0("- reductions: ", paste(Reductions(obj), collapse = ", ")),
  paste0("- images: ", paste(names(obj@images), collapse = ", ")),
  paste0("- sample column: ", ifelse(is.na(sample_col), "NOT FOUND", sample_col)),
  paste0("- n samples/slices: ", ifelse(is.na(sample_col), "NA", length(unique(meta[[sample_col]])))),
  "",
  "Interpretation boundary:",
  "- This audit only verifies the legacy spatial object structure.",
  "- It does not rerun CARD, subtype mapping, spatial correlation, CellChat, Fuzzy Cosine, or pseudotime.",
  "- Old `CARD_*` columns, if present, are legacy labels and must not be treated as current MES-V/MES-I evidence without remapping."
)
writeLines(summary_lines, file.path(root, "docs", "stage0_spatial_object_audit.md"), useBytes = TRUE)

cat("\n[STOP R9 Stage0 object audit]\n")
