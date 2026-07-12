#!/usr/bin/env Rscript
# =============================================================================
# R9 Batch 3 B1 | prepare COMMOT input matrices
#
# Boundary:
#   This only exports per-slice Visium counts and locked R9 scores for the
#   spatial-aware COMMOT run. It does not run communication inference and does
#   not make any biological claim.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
})

base_dir <- getwd()
out_dir <- file.path(base_dir, "tables/R9_batch3_B1_spatial_communication/input")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

spatial_obj_path <- "<DATA_ROOT>/项目/分型/分型代码/0.对象/5.ST_merge.rds"
r5_dir <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/07_细胞通讯/tables"
score_path <- file.path(base_dir, "tables/R9_stage3_signature_scores/R9_stage3_IVY_hypoxia_scores_per_spot.csv")

primary <- data.table(
  tier = "primary",
  interaction_name = "PLAU_PLAUR",
  ligand = "PLAU",
  receptor = "PLAUR",
  pathway_name = "PLAU",
  source_table = file.path(r5_dir, "PLAU_PLAUR_纯TME方向表达支撑表.csv"),
  rule = "pre-declared primary LR"
)

plaur_db <- fread(file.path(r5_dir, "Phase0_CellChatDB_PLAUR相关互作.csv"))
secondary <- unique(plaur_db[interaction_name != "PLAU_PLAUR",
  .(interaction_name, ligand, receptor, pathway_name)
])
if (nrow(secondary) > 0) {
  secondary[, `:=`(
    tier = "secondary",
    source_table = file.path(r5_dir, "Phase0_CellChatDB_PLAUR相关互作.csv"),
    rule = "pre-existing R5 CellChatDB PLAUR-related LR, excluding primary"
  )]
  setcolorder(secondary, names(primary))
} else {
  secondary <- primary[0]
}

all_comm <- fread(file.path(r5_dir, "CellChat纯TME_全部通讯结果.csv"))
exploratory <- unique(all_comm[, .(interaction_name, ligand, receptor, pathway_name)])
exploratory <- exploratory[interaction_name != "PLAU_PLAUR"]
exploratory[, `:=`(
  tier = "exploratory",
  source_table = file.path(r5_dir, "CellChat纯TME_全部通讯结果.csv"),
  rule = "unique LR from R5 pure-TME CellChat result; locked internal-audit/background landscape"
)]
setcolorder(exploratory, names(primary))

lr_all <- unique(rbindlist(list(primary, secondary, exploratory), use.names = TRUE, fill = TRUE))

split_complex <- function(x) {
  unique(unlist(strsplit(x, "_", fixed = TRUE), use.names = FALSE))
}
candidate_genes <- unique(c(lr_all$ligand, unlist(lapply(lr_all$receptor, split_complex), use.names = FALSE)))

message("Reading Spatial object: ", spatial_obj_path)
st <- readRDS(spatial_obj_path)
if (!"Spatial" %in% names(st@assays)) stop("Spatial assay not found.")

spatial_genes <- rownames(st[["Spatial"]])
lr_all[, ligand_in_spatial := ligand %in% spatial_genes]
lr_all[, receptor_all_in_spatial := vapply(receptor, function(r) all(split_complex(r) %in% spatial_genes), logical(1))]
lr_all[, usable_in_spatial := ligand_in_spatial & receptor_all_in_spatial]
lr_usable <- lr_all[usable_in_spatial == TRUE]
if (!any(lr_usable$tier == "primary")) stop("Primary PLAU_PLAUR is not covered in Spatial genes.")

fwrite(lr_all, file.path(out_dir, "B1_LR_database_all_with_spatial_coverage.csv"))
fwrite(lr_usable, file.path(out_dir, "B1_LR_database_usable_for_COMMOT.csv"))

genes_needed <- unique(c(lr_usable$ligand, unlist(lapply(lr_usable$receptor, split_complex), use.names = FALSE)))
genes_needed <- intersect(genes_needed, spatial_genes)

scores <- fread(score_path)
required_score_cols <- c("spot_id", "slice", "image", "x", "y", "MES-V", "MES-lineage", "vascular", "neuron_control")
missing_scores <- setdiff(required_score_cols, names(scores))
if (length(missing_scores)) stop("Missing score columns: ", paste(missing_scores, collapse = ", "))

image_slice_map <- data.table(
  image = names(st@images),
  slice = vapply(names(st@images), function(img) {
    cells <- Cells(st@images[[img]])
    si <- unique(st$orig.ident[cells])
    if (length(si) == 1) as.character(si) else NA_character_
  }, character(1))
)[!is.na(slice)]
fwrite(image_slice_map, file.path(out_dir, "B1_image_slice_map.csv"))

manifest <- list()
ri <- 0L
for (sl in sort(unique(scores$slice))) {
  img <- image_slice_map[slice == sl, image][1]
  if (is.na(img)) stop("No image mapping for slice: ", sl)
  cells <- Cells(st@images[[img]])
  layer_candidates <- c(paste0("counts.", sl), paste0("counts.", gsub("^#", "", sl)), "counts")
  layer_candidates <- layer_candidates[layer_candidates %in% Layers(st[["Spatial"]])]
  if (!length(layer_candidates)) stop("No counts layer for slice: ", sl)
  mat <- LayerData(st[["Spatial"]], layer = layer_candidates[1])
  present_cells <- intersect(cells, colnames(mat))
  mat <- mat[genes_needed, present_cells, drop = FALSE]
  slice_scores <- scores[slice == sl & spot_id %in% present_cells, ..required_score_cols]
  slice_scores <- slice_scores[match(present_cells, spot_id)]
  if (!identical(slice_scores$spot_id, present_cells)) stop("Spot order mismatch for slice: ", sl)

  safe_sl <- gsub("[^A-Za-z0-9]+", "_", gsub("^#", "", sl))
  slice_dir <- file.path(out_dir, safe_sl)
  dir.create(slice_dir, showWarnings = FALSE, recursive = TRUE)
  Matrix::writeMM(t(mat), file.path(slice_dir, "counts_spots_by_genes.mtx"))
  fwrite(data.table(gene = rownames(mat)), file.path(slice_dir, "genes.csv"))
  fwrite(slice_scores, file.path(slice_dir, "obs_scores.csv"))
  ri <- ri + 1L
  manifest[[ri]] <- data.table(
    slice = sl,
    safe_slice = safe_sl,
    slice_dir = slice_dir,
    n_spots = ncol(mat),
    n_genes = nrow(mat),
    counts_layer = layer_candidates[1]
  )
}

manifest <- rbindlist(manifest, use.names = TRUE)
fwrite(manifest, file.path(out_dir, "B1_COMMOT_input_manifest.csv"))

params <- data.table(
  parameter = c("spatial_object", "score_table", "R5_primary_source", "R5_secondary_source", "R5_exploratory_source", "exported_genes_n"),
  value = c(
    spatial_obj_path,
    score_path,
    primary$source_table[1],
    ifelse(nrow(secondary) > 0, secondary$source_table[1], "none"),
    file.path(r5_dir, "CellChat纯TME_全部通讯结果.csv"),
    as.character(length(genes_needed))
  )
)
fwrite(params, file.path(out_dir, "B1_COMMOT_input_parameters.csv"))

message("Done. Output: ", out_dir)
