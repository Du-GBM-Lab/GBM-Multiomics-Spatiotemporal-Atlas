# F1_fine_cluster_correlation_k4.R
# Purpose: Draw fine-cluster signature correlation heatmap supporting the
# coarse-grained k=4 malignant subtype framework.
# Input:
#   tables/02c_sweep_theta5/cluster_correlation_matrix.csv
#   tables/03c_cluster_to_subtype_k4_theta5.csv
# Output:
#   figures/F1_fine_cluster_correlation_k4.pdf
#   figures/source_data/F1_correlation_matrix.csv
#   figures/source_data/F1_cluster_subtype_mapping.csv
#
# Interpretation boundary:
#   This figure shows the fine-cluster correlation structure and the selected
#   k=4 coarse partition. Do not describe gap/silhouette as determining k=4.

.libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

base_dir <- "05_恶性细胞分亚群与Neftel对照"
fig_dir <- file.path(base_dir, "figures")
source_dir <- file.path(fig_dir, "source_data")
dir.create(source_dir, showWarnings = FALSE, recursive = TRUE)

cor_path <- file.path(base_dir, "tables/02c_sweep_theta5/cluster_correlation_matrix.csv")
map_path <- file.path(base_dir, "tables/03c_cluster_to_subtype_k4_theta5.csv")
out_pdf <- file.path(fig_dir, "F1_fine_cluster_correlation_k4.pdf")

subtype_colors <- c(
  "Subtype1" = "#0072B5",
  "Subtype2" = "#E18727",
  "Subtype3" = "#20854E",
  "Subtype4" = "#BC3C29"
)

cor_color <- circlize::colorRamp2(
  c(0, 0.20, 0.35, 0.50, 0.75, 1),
  c("#2166AC", "#92C5DE", "white", "#F4A582", "#D6604D", "#B2182B")
)

adjusted_rand_index <- function(x, y) {
  x <- as.factor(x)
  y <- as.factor(y)
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  n <- sum(tab)
  sum_comb <- sum(choose2(tab))
  row_comb <- sum(choose2(rowSums(tab)))
  col_comb <- sum(choose2(colSums(tab)))
  total_comb <- choose2(n)
  expected <- row_comb * col_comb / total_comb
  max_index <- (row_comb + col_comb) / 2
  if (isTRUE(all.equal(max_index, expected))) return(0)
  (sum_comb - expected) / (max_index - expected)
}

cor_mat <- as.matrix(read.csv(cor_path, row.names = 1, check.names = FALSE))
storage.mode(cor_mat) <- "numeric"
mapping <- read.csv(map_path, check.names = FALSE, stringsAsFactors = FALSE)

cluster_col <- if ("fine_cluster_id" %in% colnames(mapping)) "fine_cluster_id" else "cluster"
subtype_col <- if ("subtype_k4_final" %in% colnames(mapping)) "subtype_k4_final" else "subtype"

mapping[[cluster_col]] <- as.character(mapping[[cluster_col]])
mapping[[subtype_col]] <- as.character(mapping[[subtype_col]])

stopifnot(nrow(cor_mat) == 26, ncol(cor_mat) == 26)
stopifnot(setequal(rownames(cor_mat), colnames(cor_mat)))
stopifnot(setequal(rownames(cor_mat), mapping[[cluster_col]]))
stopifnot(all(mapping[[subtype_col]] %in% names(subtype_colors)))

cor_mat <- cor_mat[rownames(cor_mat), rownames(cor_mat), drop = FALSE]
cluster_to_subtype <- setNames(mapping[[subtype_col]], mapping[[cluster_col]])
cluster_subtype_vec <- cluster_to_subtype[rownames(cor_mat)]
cluster_subtype_vec <- factor(cluster_subtype_vec, levels = names(subtype_colors))
stopifnot(!any(is.na(cluster_subtype_vec)))

methods_to_try <- c("average", "ward.D2", "complete", "ward.D")
match_results <- vapply(methods_to_try, function(m) {
  hc <- hclust(as.dist(1 - cor_mat), method = m)
  cuts <- cutree(hc, k = 4)
  adjusted_rand_index(cuts, cluster_subtype_vec)
}, numeric(1))

cat("ARI by hclust method:\n")
print(round(match_results, 4))

best_method <- names(which.max(match_results))
ari_best <- max(match_results)
cat(sprintf("Best method: %s (ARI = %.4f)\n", best_method, ari_best))

if (isTRUE(all.equal(ari_best, 1, tolerance = 1e-8))) {
  hc_use <- hclust(as.dist(1 - cor_mat), method = best_method)
  row_split_use <- 4
  split_strategy <- sprintf("automatic_cutree_k4_%s", best_method)
} else {
  hc_use <- hclust(as.dist(1 - cor_mat), method = "average")
  row_split_use <- cluster_subtype_vec
  split_strategy <- "explicit_subtype_split_with_average_dendrogram"
}

cat(sprintf("Used hclust method: %s\n", hc_use$method))
cat(sprintf("Row/column split strategy: %s\n", split_strategy))
cat("Fine cluster x subtype mapping:\n")
print(data.frame(
  fine_cluster = rownames(cor_mat),
  subtype = as.character(cluster_subtype_vec),
  row.names = NULL
))

row_anno <- rowAnnotation(
  Subtype = cluster_subtype_vec,
  col = list(Subtype = subtype_colors),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  width = unit(2.5, "mm")
)

col_anno <- HeatmapAnnotation(
  Subtype = cluster_subtype_vec,
  col = list(Subtype = subtype_colors),
  show_legend = TRUE,
  show_annotation_name = FALSE,
  height = unit(2.5, "mm"),
  annotation_legend_param = list(
    Subtype = list(
      title = "Subtype",
      title_gp = gpar(fontsize = 9.5, fontface = "bold"),
      labels_gp = gpar(fontsize = 8.5),
      grid_height = unit(3, "mm"),
      grid_width = unit(3, "mm")
    )
  )
)

ht <- Heatmap(
  cor_mat,
  name = "Cluster\ncorrelation",
  col = cor_color,
  cluster_rows = hc_use,
  cluster_columns = hc_use,
  row_split = row_split_use,
  column_split = row_split_use,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  row_dend_width = unit(8, "mm"),
  column_dend_height = unit(8, "mm"),
  row_names_gp = gpar(fontsize = 8.5),
  column_names_gp = gpar(fontsize = 8.5),
  row_names_side = "right",
  column_names_side = "bottom",
  row_title = NULL,
  column_title = NULL,
  top_annotation = col_anno,
  left_annotation = row_anno,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.3),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 9.5, fontface = "bold"),
    labels_gp = gpar(fontsize = 8.5),
    legend_height = unit(2.8, "cm"),
    grid_width = unit(3, "mm")
  )
)

pdf(out_pdf, width = 6.2, height = 5.8, useDingbats = FALSE)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(3, 3, 3, 5), "mm")
)
dev.off()

write.csv(cor_mat, file.path(source_dir, "F1_correlation_matrix.csv"))
write.csv(
  data.frame(
    fine_cluster = rownames(cor_mat),
    subtype = as.character(cluster_subtype_vec),
    stringsAsFactors = FALSE
  ),
  file.path(source_dir, "F1_cluster_subtype_mapping.csv"),
  row.names = FALSE
)

cat("Output PDF:", out_pdf, "\n")
cat("Source data directory:", source_dir, "\n")
cat("Sanity checks passed.\n")
