suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir) {
  setwd(file.path(getwd(), step_dir))
}

score_path <- "tables/10c_metaprogram_score_by_subtype.csv"
label_path <- "tables/10b_metaprogram_label_final.csv"
out_pdf <- "figures/F5_nmf_metaprogram_heatmap.pdf"
out_source <- "figures/source_data/F5_nmf_metaprogram_heatmap_source.csv"
out_z <- "figures/source_data/F5_nmf_metaprogram_heatmap_zscore_matrix.csv"
out_sanity <- "figures/source_data/F5_nmf_metaprogram_heatmap_sanity_checks.csv"

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/source_data", showWarnings = FALSE, recursive = TRUE)

subtype_levels <- paste0("Subtype", 1:4)
subtype_colors <- c(
  "Subtype1" = "#0072B5",
  "Subtype2" = "#E18727",
  "Subtype3" = "#20854E",
  "Subtype4" = "#BC3C29"
)

mp_levels <- paste0("MP0", 1:6)
mp_label_map <- c(
  "MP01" = "MP01: IGFBP-signaling / Stress-response",
  "MP02" = "MP02: ECM-organization",
  "MP03" = "MP03: Myelination",
  "MP04" = "MP04: MHC-II Antigen-presentation",
  "MP05" = "MP05: Cell-cycle",
  "MP06" = "MP06: Angiogenesis-signaling"
)
mp_family <- c(
  "MP01" = "Stress",
  "MP02" = "ECM",
  "MP03" = "Myelination",
  "MP04" = "Antigen-presentation",
  "MP05" = "Cell-cycle",
  "MP06" = "Angiogenesis"
)

score_tbl <- readr::read_csv(score_path, show_col_types = FALSE)
label_tbl <- readr::read_csv(label_path, show_col_types = FALSE)

score_tbl <- score_tbl %>%
  mutate(
    metaprogram = ifelse(metaprogram == "MP05_cycling_supplement", "MP05", metaprogram),
    subtype_k4 = factor(subtype_k4, levels = subtype_levels)
  ) %>%
  filter(metaprogram %in% mp_levels)

stopifnot(setequal(unique(score_tbl$metaprogram), mp_levels))
stopifnot(setequal(as.character(unique(score_tbl$subtype_k4)), subtype_levels))

mean_mat <- score_tbl %>%
  select(subtype_k4, metaprogram, mean) %>%
  tidyr::pivot_wider(names_from = subtype_k4, values_from = mean) %>%
  tibble::column_to_rownames("metaprogram") %>%
  as.matrix()
mean_mat <- mean_mat[mp_levels, subtype_levels, drop = FALSE]
stopifnot(!any(is.na(mean_mat)))

z_mat <- t(scale(t(mean_mat)))
z_mat[is.na(z_mat)] <- 0
z_mat[z_mat > 2.5] <- 2.5
z_mat[z_mat < -2.5] <- -2.5
rownames(z_mat) <- mp_label_map[rownames(z_mat)]

row_family <- factor(mp_family[mp_levels], levels = unique(mp_family[mp_levels]))
column_subtype <- factor(colnames(z_mat), levels = subtype_levels)

col_anno <- HeatmapAnnotation(
  Subtype = column_subtype,
  col = list(Subtype = subtype_colors),
  show_annotation_name = FALSE,
  height = unit(3, "mm"),
  annotation_legend_param = list(
    Subtype = list(
      title = "Subtype",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      grid_height = unit(3, "mm"),
      grid_width = unit(3, "mm")
    )
  )
)

heat_col <- circlize::colorRamp2(c(-2.5, 0, 2.5), c("#4393C3", "white", "#D6604D"))

ht <- Heatmap(
  z_mat,
  name = "Z-score",
  col = heat_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_family,
  column_split = column_subtype,
  top_annotation = col_anno,
  row_title = NULL,
  column_title = NULL,
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 8),
  row_names_side = "right",
  row_gap = unit(0.8, "mm"),
  column_gap = unit(0.8, "mm"),
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.25),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", z_mat[i, j]), x, y, gp = gpar(fontsize = 6.5))
  },
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7),
    legend_height = unit(2.7, "cm"),
    grid_width = unit(3, "mm")
  )
)

pdf(out_pdf, width = 6.8, height = 4.8, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(3, 3, 3, 3), "mm"))
dev.off()

source_tbl <- as.data.frame(as.table(mean_mat), stringsAsFactors = FALSE) %>%
  rename(metaprogram = Var1, subtype_k4 = Var2, mean_score = Freq) %>%
  mutate(
    metaprogram_label = mp_label_map[metaprogram],
    program_family = mp_family[metaprogram],
    subtype_k4 = factor(subtype_k4, levels = subtype_levels),
    z_score = mapply(
      function(mp, st) z_mat[mp_label_map[mp], as.character(st)],
      metaprogram,
      subtype_k4
    )
  ) %>%
  arrange(factor(metaprogram, levels = mp_levels), subtype_k4)

readr::write_csv(source_tbl, out_source)
write.csv(z_mat, out_z)

sanity_tbl <- tibble::tibble(
  metric = c("n_metaprograms", "n_subtypes", "missing_mean_values", "labels_from_file"),
  value = c(
    nrow(mean_mat),
    ncol(mean_mat),
    sum(is.na(mean_mat)),
    paste(label_tbl$metaprogram_id, collapse = ";")
  )
)
readr::write_csv(sanity_tbl, out_sanity)

cat("F5 NMF metaprogram heatmap completed.\n")
cat("Mean metaprogram score matrix:\n")
print(round(mean_mat, 3))
cat("Output PDF:", out_pdf, "\n")
cat("Source:", out_source, "\n")
