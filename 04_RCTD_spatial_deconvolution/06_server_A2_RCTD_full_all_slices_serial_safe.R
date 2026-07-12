suppressPackageStartupMessages({
  library(spacexr)
  library(qs2)
  library(Seurat)
  library(Matrix)
  library(data.table)
})

set.seed(1)

root <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
path_st  <- file.path(root, "data/spatial/5.ST_merge.rds")
ref_path <- file.path(root, "outputs/R9_A1_RCTD_reference.qs2")
out_dir  <- file.path(root, "outputs")
tab_dir  <- file.path(root, "tables")
slice_dir <- file.path(out_dir, "R9_A2_per_slice")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(slice_dir, showWarnings = FALSE, recursive = TRUE)

n_parallel <- 1
max_cores_each <- 1

image_to_orig <- function(image_name) {
  sub("^X\\.", "#", image_name)
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
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

get_spatial_layer <- function(st, orig_ident) {
  layer_name <- paste0("counts.", orig_ident)
  layers <- Layers(st[["Spatial"]])
  if (!layer_name %in% layers) {
    stop("Missing Spatial raw count layer: ", layer_name,
         " | available layers: ", paste(layers, collapse = ","))
  }
  out <- LayerData(st[["Spatial"]], layer = layer_name)
  if (!inherits(out, "dgCMatrix")) out <- as(out, "dgCMatrix")
  out
}

cat("== spacexr version:", as.character(packageVersion("spacexr")), "\n")
cat("== A2 mode: serial-safe | n_parallel =", n_parallel,
    "| max_cores_each =", max_cores_each, "\n")

reference <- qs2::qs_read(ref_path)
st <- readRDS(path_st)

imgs <- Images(st)
oid <- unique(as.character(st$orig.ident))
key_img <- sub("^X\\.", "", imgs)
key_oid <- sub("^#", "", oid)
map <- data.frame(
  image = imgs,
  image_key = key_img,
  orig.ident = oid[match(key_img, key_oid)],
  orig_key = key_oid[match(key_img, key_oid)],
  n_spots_meta = as.integer(table(factor(st$orig.ident, levels = oid))[match(key_img, key_oid)]),
  stringsAsFactors = FALSE
)
data.table::fwrite(map, file.path(tab_dir, "R9_A2_image_orig_ident_mapping.csv"))
cat("== image <-> orig.ident mapping:\n")
print(map)
stopifnot(nrow(map) == 18, !any(is.na(map$orig.ident)), identical(map$image_key, map$orig_key))

run_one <- function(i) {
  img <- map$image[i]
  oi <- map$orig.ident[i]
  out_path <- file.path(slice_dir, paste0("R9_A2_weights_", safe_name(oi), ".qs2"))
  if (file.exists(out_path)) {
    cat("\n== [skip existing]", oi, "\n")
    old <- qs2::qs_read(out_path)
    return(old$summary)
  }

  cat("\n== [start]", i, "/", nrow(map), img, oi, "\n")
  t0 <- Sys.time()
  ans <- tryCatch({
    cells_by_orig <- colnames(st)[st$orig.ident == oi]
    cells_by_image <- Cells(st@images[[img]])
    if (!setequal(cells_by_orig, cells_by_image)) {
      stop("Spot mismatch between orig.ident and image for ", oi)
    }
    cells <- cells_by_image

    sp_counts <- get_spatial_layer(st, oi)
    sp_counts <- sp_counts[, cells, drop = FALSE]

    co_raw <- GetTissueCoordinates(st, image = img)
    co <- extract_xy(co_raw)
    co <- co[colnames(sp_counts), , drop = FALSE]
    if (!identical(rownames(co), colnames(sp_counts))) {
      stop("Coordinate/count spot order mismatch for ", oi)
    }

    nUMI <- Matrix::colSums(sp_counts)
    keep <- nUMI > 0
    puck <- SpatialRNA(co[keep, , drop = FALSE], sp_counts[, keep, drop = FALSE], nUMI[keep])

    rc <- create.RCTD(puck, reference, max_cores = max_cores_each)
    rc <- run.RCTD(rc, doublet_mode = "full")
    w <- as.matrix(rc@results$weights)
    if (nrow(w) != sum(keep)) {
      stop("Unexpected weights rows: ", nrow(w), " != kept spots ", sum(keep))
    }
    rs <- rowSums(w)
    if (any(rs <= 0 | is.na(rs))) stop("Found non-positive/NA weight row sums")
    w_norm <- sweep(w, 1, rs, "/")

    df <- data.frame(
      spot_id = rownames(w_norm),
      slice = oi,
      image = img,
      x = co[rownames(w_norm), "x"],
      y = co[rownames(w_norm), "y"],
      w_norm,
      check.names = FALSE
    )
    minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    summary <- data.frame(
      ok = TRUE,
      image = img,
      slice = oi,
      n_spot_in = length(cells),
      n_spot_nUMI_positive = sum(keep),
      n_spot_out = nrow(w_norm),
      n_reject = length(cells) - nrow(w_norm),
      minutes = minutes,
      message = "",
      stringsAsFactors = FALSE
    )
    qs2::qs_save(list(weights = df, summary = summary), out_path)
    cat("== [done]", oi, "| spots out:", nrow(w_norm),
        "| minutes:", round(minutes, 2), "\n")
    summary
  }, error = function(e) {
    minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    summary <- data.frame(
      ok = FALSE,
      image = img,
      slice = oi,
      n_spot_in = NA_integer_,
      n_spot_nUMI_positive = NA_integer_,
      n_spot_out = NA_integer_,
      n_reject = NA_integer_,
      minutes = minutes,
      message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
    qs2::qs_save(list(weights = NULL, summary = summary), out_path)
    cat("== [fail]", oi, "|", conditionMessage(e), "\n")
    summary
  })
  ans
}

summaries <- lapply(seq_len(nrow(map)), run_one)
summary_tab <- data.table::rbindlist(summaries, fill = TRUE)
data.table::fwrite(summary_tab, file.path(tab_dir, "R9_A2_perslice_spot_counts.csv"))

ok_slices <- summary_tab[ok == TRUE, slice]
weight_list <- lapply(ok_slices, function(oi) {
  qs2::qs_read(file.path(slice_dir, paste0("R9_A2_weights_", safe_name(oi), ".qs2")))$weights
})
big <- data.table::rbindlist(weight_list, fill = TRUE)
qs2::qs_save(big, file.path(out_dir, "R9_A2_RCTD_weights_allslices_long.qs2"))
data.table::fwrite(big, file.path(tab_dir, "R9_A2_RCTD_weights_allslices_long.csv"))

ct_cols <- setdiff(colnames(big), c("spot_id", "slice", "image", "x", "y"))
mean_prop <- data.table(
  celltype = ct_cols,
  mean_prop = as.numeric(colMeans(as.data.frame(big[, ..ct_cols])))
)
setorder(mean_prop, -mean_prop)
data.table::fwrite(mean_prop, file.path(tab_dir, "R9_A2_global_mean_proportions.csv"))

cat("\n== A2 success:", sum(summary_tab$ok), "/", nrow(summary_tab), "\n")
if (any(!summary_tab$ok)) {
  cat("== failed slices:\n")
  print(summary_tab[ok == FALSE, .(slice, message)])
}
cat("== total output spots:", nrow(big), "\n")
cat("== per-slice spot counts:\n")
print(summary_tab)
cat("== global mean proportions:\n")
print(mean_prop)

cat("\n[STOP A2] Full-slice RCTD serial-safe run completed or stopped with recorded per-slice failures. Do not run Step B before review.\n")
