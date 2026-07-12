#!/usr/bin/env Rscript
# =============================================================================
# R9 batch 1 | A1 SVG with SPARK-X on server
# Purpose:
#   Identify spatial variable genes per Visium slice using SPARK::sparkx().
#
# Boundaries:
#   - SVGs are single-gene spatial pattern results, not proximity/interaction.
#   - Results are landscape/method evidence unless later overlap tests pass.
#   - No pooled-spot significance claim; recurrence is summarized across slices.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(Matrix)
})

set.seed(1)

base_dir <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
st_path <- file.path(base_dir, "data/spatial/5.ST_merge.rds")
out_dir <- file.path(base_dir, "outputs/R9_batch1_unbiased_landscape/A1_SVG")
tab_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A1_SVG")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

n_cores <- as.integer(Sys.getenv("R9_A1_SPARKX_CORES", "4"))
min_detect_frac <- as.numeric(Sys.getenv("R9_A1_MIN_DETECT_FRAC", "0.03"))
max_genes_per_slice <- as.integer(Sys.getenv("R9_A1_MAX_GENES", "8000"))
sparkx_option <- Sys.getenv("R9_A1_SPARKX_OPTION", "mixture")

pkg_audit <- data.table(
  package = c("SPARK", "Seurat", "Matrix", "data.table"),
  available = c(
    requireNamespace("SPARK", quietly = TRUE),
    requireNamespace("Seurat", quietly = TRUE),
    requireNamespace("Matrix", quietly = TRUE),
    requireNamespace("data.table", quietly = TRUE)
  )
)
fwrite(pkg_audit, file.path(tab_dir, "A1_SVG_package_audit.csv"))
if (!all(pkg_audit$available)) {
  stop("Missing required packages for A1 SVG: ",
       paste(pkg_audit$package[!pkg_audit$available], collapse = ", "))
}

st <- readRDS(st_path)
sample_col <- "orig.ident"
if (!sample_col %in% colnames(st@meta.data)) stop("Missing orig.ident in spatial metadata.")

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

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

imgs <- Images(st)
oid <- unique(as.character(st$orig.ident))
key_img <- sub("^X\\.", "", imgs)
key_oid <- sub("^#", "", oid)
map <- data.table(
  image = imgs,
  slice = oid[match(key_img, key_oid)]
)
if (anyNA(map$slice)) stop("Could not map image names to orig.ident.")
fwrite(map, file.path(tab_dir, "A1_SVG_image_slice_map.csv"))

run_one_slice <- function(img, slice) {
  out_csv <- file.path(tab_dir, paste0("A1_SVG_SPARKX_", safe_name(slice), ".csv"))
  if (file.exists(out_csv)) {
    cat("[skip existing]", slice, "\n")
    return(fread(out_csv))
  }

  layer_name <- paste0("counts.", slice)
  if (!layer_name %in% Layers(st[["Spatial"]])) stop("Missing Spatial layer: ", layer_name)
  cells <- Cells(st@images[[img]])
  counts <- LayerData(st[["Spatial"]], layer = layer_name)
  counts <- counts[, cells, drop = FALSE]
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")

  coords <- extract_xy(GetTissueCoordinates(st, image = img))
  coords <- coords[colnames(counts), , drop = FALSE]
  if (!identical(rownames(coords), colnames(counts))) {
    stop("Coordinate/count order mismatch for slice: ", slice)
  }

  detected_frac <- Matrix::rowMeans(counts > 0)
  gene_mean <- Matrix::rowMeans(counts)
  gene_var <- Matrix::rowMeans(counts^2) - gene_mean^2
  keep <- detected_frac >= min_detect_frac & gene_var > 0
  if (sum(keep) > max_genes_per_slice) {
    ord <- order(gene_var[keep], decreasing = TRUE)
    keep_genes <- names(gene_var[keep])[ord[seq_len(max_genes_per_slice)]]
  } else {
    keep_genes <- names(which(keep))
  }
  if (length(keep_genes) < 50) stop("Too few genes after filtering for slice: ", slice)

  count_in <- counts[keep_genes, , drop = FALSE]
  locus_in <- as.matrix(coords[, c("x", "y"), drop = FALSE])

  cat("[SPARK-X start]", slice, "spots=", ncol(count_in), "genes=", nrow(count_in), "\n")
  fit <- SPARK::sparkx(
    count_in = count_in,
    locus_in = locus_in,
    numCores = n_cores,
    option = sparkx_option,
    verbose = TRUE
  )
  res <- as.data.table(fit$res_mtest, keep.rownames = "gene")
  setnames(res, old = names(res), new = make.names(names(res), unique = TRUE))
  res[, slice := slice]
  res[, n_spots := ncol(count_in)]
  res[, n_genes_tested := nrow(count_in)]
  res[, detected_fraction := detected_frac[gene]]
  fwrite(res, out_csv)
  res
}

all_res <- rbindlist(lapply(seq_len(nrow(map)), function(i) run_one_slice(map$image[i], map$slice[i])), fill = TRUE)
fwrite(all_res, file.path(tab_dir, "A1_SVG_SPARKX_all_slices_long.csv"))

padj_col <- intersect(c("adjustedPval", "adjusted.P.value", "combinedPval.adj", "Pval.adj", "BY.adjusted.P.value"), names(all_res))
if (!length(padj_col)) {
  padj_col <- grep("adj|adjust", names(all_res), value = TRUE, ignore.case = TRUE)
}
p_col <- grep("pval|p.value", names(all_res), value = TRUE, ignore.case = TRUE)
summary <- all_res[, .(
  n_slices_tested = uniqueN(slice),
  n_slices_fdr_005 = if (length(padj_col)) sum(get(padj_col[1]) < 0.05, na.rm = TRUE) else NA_integer_,
  min_p = if (length(p_col)) min(get(p_col[1]), na.rm = TRUE) else NA_real_
), by = gene][order(-n_slices_fdr_005, min_p)]
fwrite(summary, file.path(tab_dir, "A1_SVG_SPARKX_consensus_summary.csv"))

params <- data.table(
  parameter = c("method", "function", "n_cores", "min_detect_frac", "max_genes_per_slice", "sparkx_option", "statistical_unit", "boundary"),
  value = c("SPARK-X", "SPARK::sparkx", n_cores, min_detect_frac, max_genes_per_slice, sparkx_option,
            "per-slice SVG; cross-slice recurrence summarized after per-slice tests",
            "SVG landscape only; no proximity/contact/interaction claim")
)
fwrite(params, file.path(tab_dir, "A1_SVG_parameters.csv"))

cat("\n[STOP A1 SVG] SPARK-X per-slice source tables written. Review before overlap/figure decisions.\n")
