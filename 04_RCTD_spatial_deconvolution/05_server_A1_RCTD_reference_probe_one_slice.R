suppressPackageStartupMessages({
  library(spacexr)
  library(qs2)
  library(Seurat)
  library(Matrix)
  library(data.table)
})

set.seed(1)

root <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
path_full <- file.path(root, "data/reference/GBM.RNA.qc_doubletfinder.infercnv_immune_reference_calls.qs2")
path_mal  <- file.path(root, "data/reference/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2")
path_st   <- file.path(root, "data/spatial/5.ST_merge.rds")

out_dir <- file.path(root, "outputs")
tab_dir <- file.path(root, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

get_counts <- function(obj, assay, layer = "counts") {
  out <- tryCatch(
    GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) GetAssayData(obj, assay = assay, slot = layer)
  )
  if (!inherits(out, "dgCMatrix")) out <- as(out, "dgCMatrix")
  out
}

image_to_orig <- function(image_name) {
  sub("^X\\.", "#", image_name)
}

extract_xy <- function(coords) {
  if ("cell" %in% colnames(coords)) rownames(coords) <- coords$cell
  if (all(c("x", "y") %in% colnames(coords))) {
    xy <- coords[, c("x", "y"), drop = FALSE]
  } else if (all(c("imagecol", "imagerow") %in% colnames(coords))) {
    xy <- coords[, c("imagecol", "imagerow"), drop = FALSE]
    colnames(xy) <- c("x", "y")
  } else if (all(c("col", "row") %in% colnames(coords))) {
    xy <- coords[, c("col", "row"), drop = FALSE]
    colnames(xy) <- c("x", "y")
  } else {
    stop("Unsupported coordinate columns: ", paste(colnames(coords), collapse = ","))
  }
  xy
}

cat("== spacexr version:", as.character(packageVersion("spacexr")), "\n")
cat("== API branch: CLASSIC create.RCTD/run.RCTD\n")

## 1) Build and save RCTD reference.
full <- qs2::qs_read(path_full)
DefaultAssay(full) <- "RNA"
mal <- qs2::qs_read(path_mal)

stopifnot("anno_ident" %in% colnames(full@meta.data))
stopifnot("subtype_k4" %in% colnames(mal@meta.data))

ref_counts <- get_counts(full, assay = "RNA", layer = "counts")
ct <- as.character(full$anno_ident)
names(ct) <- colnames(full)
ov <- intersect(colnames(mal), colnames(full))
stopifnot(length(ov) == ncol(mal))
ct[ov] <- as.character(mal$subtype_k4[match(ov, colnames(mal))])
ct <- factor(ct)
nUMI_ref <- Matrix::colSums(ref_counts)

celltype_table <- sort(table(ct), decreasing = TRUE)
cat("== reference celltypes:\n")
print(celltype_table)
stopifnot(all(celltype_table > 25))

reference <- Reference(ref_counts, ct, nUMI_ref)
reference_path <- file.path(out_dir, "R9_A1_RCTD_reference.qs2")
qs2::qs_save(reference, reference_path)
cat("== Reference saved:", reference_path, "\n")

## 2) Prepare one Visium slice.
st <- readRDS(path_st)
imgs <- Images(st)
stopifnot(length(imgs) == 18)
probe_image <- imgs[1]
probe_orig <- image_to_orig(probe_image)

cells_by_orig <- colnames(st)[st$orig.ident == probe_orig]
cells_by_image <- Cells(st@images[[probe_image]])
cat("== probe image:", probe_image, "| orig.ident:", probe_orig, "\n")
cat("== spots by orig.ident:", length(cells_by_orig), "| by image:", length(cells_by_image), "\n")
cat("== spot setequal:", setequal(cells_by_orig, cells_by_image),
    "| identical order:", identical(cells_by_orig, cells_by_image), "\n")
stopifnot(setequal(cells_by_orig, cells_by_image))

layer_name <- paste0("counts.", probe_orig)
sp_counts <- LayerData(st[["Spatial"]], layer = layer_name)
if (!inherits(sp_counts, "dgCMatrix")) sp_counts <- as(sp_counts, "dgCMatrix")
sp_counts <- sp_counts[, cells_by_image, drop = FALSE]

coords_raw <- GetTissueCoordinates(st, image = probe_image)
coords <- extract_xy(coords_raw)
coords <- coords[colnames(sp_counts), , drop = FALSE]
stopifnot(!anyNA(coords$x), !anyNA(coords$y), identical(rownames(coords), colnames(sp_counts)))

nUMI_sp <- Matrix::colSums(sp_counts)
cat("== probe slice counts:", nrow(sp_counts), "genes x", ncol(sp_counts), "spots\n")
cat("== coordinate columns used: x,y | nrow:", nrow(coords), "\n")

puck <- SpatialRNA(coords, sp_counts, nUMI_sp)

## 3) Run one-slice RCTD full mode.
t0 <- Sys.time()
RCTD <- create.RCTD(puck, reference, max_cores = 4)
RCTD <- run.RCTD(RCTD, doublet_mode = "full")
dt <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

w <- RCTD@results$weights
w <- as.matrix(w)
cat("== raw weights dim:", paste(dim(w), collapse = "x"), "\n")
stopifnot(nrow(w) == ncol(sp_counts))
rs <- rowSums(w)
if (any(rs <= 0 | is.na(rs))) stop("Found non-positive/NA row sums in weights")
w_norm <- sweep(w, 1, rs, "/")

cat("\n== single-slice full minutes:", round(dt, 2), "\n")
cat("== normalized weights dim:", paste(dim(w_norm), collapse = "x"), "\n")
cat("== celltype columns:\n")
print(colnames(w_norm))
cat("== mean proportions:\n")
print(round(sort(colMeans(w_norm), decreasing = TRUE), 4))
cat("== extrapolated 18-slice serial minutes:", round(dt * 18, 1), "\n")

probe_out <- list(
  slice_image = probe_image,
  slice_orig_ident = probe_orig,
  n_spots = ncol(sp_counts),
  minutes = dt,
  weights_raw = w,
  weights_norm = w_norm,
  celltype_table = celltype_table,
  coordinate_columns = colnames(coords_raw)
)
qs2::qs_save(probe_out, file.path(out_dir, "R9_A1_probe_oneSlice.qs2"))
data.table::fwrite(
  data.table(spot = rownames(w_norm), w_norm),
  file.path(tab_dir, paste0("R9_A1_probe_weights_norm_", gsub("[^A-Za-z0-9]+", "_", probe_orig), ".csv"))
)

summary_tab <- data.table(
  slice_image = probe_image,
  slice_orig_ident = probe_orig,
  n_spots = ncol(sp_counts),
  n_celltypes = ncol(w_norm),
  minutes = dt,
  extrapolated_18_serial_minutes = dt * 18
)
data.table::fwrite(summary_tab, file.path(tab_dir, "R9_A1_probe_summary.csv"))

cat("\n[STOP A1] Reference + one-slice RCTD full probe completed. Do not run all slices before review.\n")
