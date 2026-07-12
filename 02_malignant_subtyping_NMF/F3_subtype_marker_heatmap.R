# F3_subtype_marker_heatmap.R
# Purpose: Draw final malignant subtype pseudobulk marker Z-score heatmap.
#
# Input:
#   outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2
#
# Output:
#   figures/F3_subtype_marker_heatmap.pdf
#   figures/source_data/F3_marker_zscore_matrix.csv
#   figures/source_data/F3_top_markers_table.csv
#   figures/source_data/F3_marker_heatmap_sanity_checks.csv

.libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(Matrix)
  library(grid)
})

set.seed(42)

base_dir <- "05_恶性细胞分亚群与Neftel对照"
fig_dir <- file.path(base_dir, "figures")
source_dir <- file.path(fig_dir, "source_data")
dir.create(source_dir, showWarnings = FALSE, recursive = TRUE)

obj_path <- file.path(base_dir, "outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2")
out_pdf <- file.path(fig_dir, "F3_subtype_marker_heatmap.pdf")
zscore_csv <- file.path(source_dir, "F3_marker_zscore_matrix.csv")
marker_csv <- file.path(source_dir, "F3_top_markers_table.csv")
sanity_csv <- file.path(source_dir, "F3_marker_heatmap_sanity_checks.csv")

subtype_levels <- c(
  "Subtype1",
  "Subtype2",
  "Subtype3",
  "Subtype4"
)

subtype_colors <- c(
  "Subtype1" = "#0072B5",
  "Subtype2" = "#E18727",
  "Subtype3" = "#20854E",
  "Subtype4" = "#BC3C29"
)

obj <- qs2::qs_read(obj_path)
DefaultAssay(obj) <- "RNA"
stopifnot("subtype_k4" %in% colnames(obj@meta.data))
obj@meta.data$subtype_k4 <- factor(obj@meta.data$subtype_k4, levels = subtype_levels)
Idents(obj) <- "subtype_k4"

cell_counts <- table(obj@meta.data$subtype_k4)
cat("Subtype cell counts:\n")
print(cell_counts)
stopifnot(setequal(names(cell_counts), subtype_levels))
stopifnot(sum(cell_counts) == 28213)
stopifnot(!anyNA(obj@meta.data$subtype_k4))

cat("Running downsampled FindAllMarkers...\n")
all_markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  max.cells.per.ident = 1000,
  random.seed = 42,
  verbose = FALSE
)

stopifnot(nrow(all_markers) > 0)
stopifnot(all(c("gene", "cluster", "p_val", "avg_log2FC") %in% colnames(all_markers)))
all_markers$BH_q <- p.adjust(all_markers$p_val, method = "BH")
all_markers$combined_score <- all_markers$avg_log2FC * (-log10(all_markers$BH_q + 1e-300))

top_markers <- all_markers %>%
  filter(BH_q < 0.05, avg_log2FC > 0.5) %>%
  mutate(cluster = factor(cluster, levels = subtype_levels)) %>%
  group_by(cluster) %>%
  arrange(desc(combined_score), .by_group = TRUE) %>%
  slice_head(n = 15) %>%
  ungroup() %>%
  group_by(gene) %>%
  slice_max(combined_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, desc(combined_score))

marker_genes <- top_markers$gene
marker_to_subtype <- setNames(as.character(top_markers$cluster), top_markers$gene)

cat(sprintf("Total markers after dedupe: %d\n", length(marker_genes)))
cat("Markers per subtype:\n")
print(table(marker_to_subtype))

stopifnot(length(marker_genes) >= 55, length(marker_genes) <= 60)
stopifnot(all(table(marker_to_subtype) >= 10))
stopifnot(!anyDuplicated(marker_genes))
stopifnot(all(marker_genes %in% rownames(obj)))

data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")[marker_genes, ]
subtype_vec <- as.character(obj@meta.data$subtype_k4)

pseudobulk <- sapply(subtype_levels, function(s) {
  cells_s <- which(subtype_vec == s)
  Matrix::rowMeans(data_mat[, cells_s, drop = FALSE])
})
colnames(pseudobulk) <- subtype_levels
rownames(pseudobulk) <- marker_genes

pseudobulk_z <- t(scale(t(pseudobulk)))
stopifnot(ncol(pseudobulk_z) == 4)
stopifnot(!any(is.na(pseudobulk_z)))
pseudobulk_z[pseudobulk_z > 2.5] <- 2.5
pseudobulk_z[pseudobulk_z < -2.5] <- -2.5

diag_checks <- lapply(subtype_levels, function(s) {
  s_markers <- names(marker_to_subtype)[marker_to_subtype == s]
  s_block <- pseudobulk_z[s_markers, , drop = FALSE]
  argmax_col <- apply(s_block, 1, which.max)
  s_col_idx <- which(colnames(pseudobulk_z) == s)
  pct_correct <- mean(argmax_col == s_col_idx) * 100
  data.frame(
    subtype = s,
    n_markers = length(s_markers),
    pct_markers_peak_in_own_subtype = pct_correct,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows()

cat("Diagonal pattern sanity:\n")
print(diag_checks)
stopifnot(all(diag_checks$pct_markers_peak_in_own_subtype >= 80))

row_anno_vec <- factor(marker_to_subtype[rownames(pseudobulk_z)], levels = subtype_levels)
col_anno_vec <- factor(colnames(pseudobulk_z), levels = subtype_levels)

row_anno <- rowAnnotation(
  Subtype = row_anno_vec,
  col = list(Subtype = subtype_colors),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  width = unit(2.5, "mm")
)

col_anno <- HeatmapAnnotation(
  Subtype = col_anno_vec,
  col = list(Subtype = subtype_colors),
  show_legend = TRUE,
  show_annotation_name = FALSE,
  height = unit(3, "mm"),
  annotation_legend_param = list(
    Subtype = list(
      title = "Subtype",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      grid_height = unit(3, "mm"),
      grid_width = unit(3, "mm"),
      nrow = 4
    )
  )
)

heat_col <- colorRamp2(c(-2.5, 0, 2.5), c("#4393C3", "white", "#D6604D"))

ht <- Heatmap(
  pseudobulk_z,
  name = "Z-score",
  col = heat_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_anno_vec,
  column_split = col_anno_vec,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 0),
  row_names_side = "right",
  row_title = NULL,
  column_title = NULL,
  row_gap = unit(0.8, "mm"),
  column_gap = unit(0.8, "mm"),
  top_annotation = col_anno,
  left_annotation = row_anno,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.2),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7),
    legend_height = unit(2.8, "cm"),
    grid_width = unit(3, "mm")
  )
)

pdf(out_pdf, width = 4.5, height = 7.5, useDingbats = FALSE)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(3, 3, 3, 5), "mm")
)
dev.off()

write.csv(pseudobulk_z, zscore_csv)
write.csv(
  top_markers %>% select(gene, cluster, avg_log2FC, BH_q, combined_score),
  marker_csv,
  row.names = FALSE
)
write.csv(diag_checks, sanity_csv, row.names = FALSE)

cat("Output PDF:", out_pdf, "\n")
cat("Z-score source data:", zscore_csv, "\n")
cat("Marker source data:", marker_csv, "\n")
cat("Sanity source data:", sanity_csv, "\n")
cat("Sanity checks passed.\n")
