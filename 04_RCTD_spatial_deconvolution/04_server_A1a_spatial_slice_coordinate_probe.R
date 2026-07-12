root <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
suppressPackageStartupMessages({
  library(Seurat)
})

path_st <- file.path(root, "data/spatial/5.ST_merge.rds")
st <- readRDS(path_st)

cat("== ST object:", nrow(st), "features x", ncol(st), "spots\n")
cat("== assays:", paste(Assays(st), collapse = ","), "\n")
cat("== default assay:", DefaultAssay(st), "\n")
cat("== meta columns:", paste(colnames(st@meta.data), collapse = ","), "\n")

imgs <- Images(st)
cat("== images:", length(imgs), "\n")
print(imgs)

cat("\n== orig.ident table:\n")
print(sort(table(st$orig.ident), decreasing = TRUE))

strip_img <- sub("^X\\.", "#", imgs)
map <- data.frame(
  image = imgs,
  stripped_to_orig_ident = strip_img,
  n_by_orig_ident = as.integer(table(factor(st$orig.ident, levels = strip_img))),
  stringsAsFactors = FALSE
)
cat("\n== image -> orig.ident mapping by stripped image name:\n")
print(map)
write.csv(map, file.path(root, "tables/A1a_image_orig_ident_mapping.csv"), row.names = FALSE)

probe_image <- imgs[1]
probe_orig <- strip_img[1]
cat("\n== probe image:", probe_image, "| orig:", probe_orig, "\n")
cells_by_orig <- colnames(st)[st$orig.ident == probe_orig]
cells_by_image <- Cells(st@images[[probe_image]])
cat("cells by orig.ident:", length(cells_by_orig), "\n")
cat("cells by image:", length(cells_by_image), "\n")
cat("setequal:", setequal(cells_by_orig, cells_by_image), "\n")
cat("identical order:", identical(cells_by_orig, cells_by_image), "\n")

coords <- GetTissueCoordinates(st, image = probe_image)
cat("\n== GetTissueCoordinates columns:\n")
print(colnames(coords))
cat("== coordinates head:\n")
print(head(coords))
cat("== coordinate rownames/cells setequal:", setequal(rownames(coords), cells_by_image), "\n")
cat("== coordinate nrow:", nrow(coords), "\n")

write.csv(
  data.frame(column = colnames(coords)),
  file.path(root, "tables/A1a_coordinate_columns.csv"),
  row.names = FALSE
)

cat("\n[STOP A1a slice coordinate probe]\n")
