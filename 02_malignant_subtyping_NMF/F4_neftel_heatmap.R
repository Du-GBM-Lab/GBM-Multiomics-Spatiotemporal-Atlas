suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir) {
  setwd(file.path(getwd(), step_dir))
}

obj_path <- "outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"
out_pdf <- "figures/F4_neftel_heatmap.pdf"
out_source <- "figures/source_data/F4_neftel_heatmap_source.csv"
out_z <- "figures/source_data/F4_neftel_heatmap_zscore_matrix.csv"
out_sanity <- "figures/source_data/F4_neftel_heatmap_sanity_checks.csv"

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/source_data", showWarnings = FALSE, recursive = TRUE)

subtype_levels <- paste0("Subtype", 1:4)
subtype_colors <- c(
  "Subtype1" = "#0072B5",
  "Subtype2" = "#E18727",
  "Subtype3" = "#20854E",
  "Subtype4" = "#BC3C29"
)

module_map <- tibble::tribble(
  ~module, ~metadata_col, ~family,
  "MES1", "AMS_MES1", "MES",
  "MES2", "AMS_MES2", "MES",
  "AC",   "AMS_AC",   "AC",
  "OPC",  "AMS_OPC",  "OPC",
  "NPC1", "AMS_NPC1", "NPC",
  "NPC2", "AMS_NPC2", "NPC",
  "G1S",  "AMS_G1S",  "Cycling",
  "G2M",  "AMS_G2M",  "Cycling"
)
module_levels <- module_map$module
family_levels <- c("MES", "AC", "OPC", "NPC", "Cycling")

obj <- qs2::qs_read(obj_path)
md <- obj@meta.data

required_cols <- c("subtype_k4", module_map$metadata_col)
missing_cols <- setdiff(required_cols, colnames(md))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

md$subtype_k4 <- factor(as.character(md$subtype_k4), levels = subtype_levels)
stopifnot(!any(is.na(md$subtype_k4)))
stopifnot(nrow(md) == 28213)

mean_tbl <- md %>%
  select(subtype_k4, all_of(module_map$metadata_col)) %>%
  group_by(subtype_k4) %>%
  summarise(
    across(all_of(module_map$metadata_col), ~ mean(.x, na.rm = TRUE)),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

mean_mat <- mean_tbl %>%
  select(subtype_k4, all_of(module_map$metadata_col)) %>%
  tibble::column_to_rownames("subtype_k4") %>%
  as.matrix() %>%
  t()

rownames(mean_mat) <- module_map$module[match(rownames(mean_mat), module_map$metadata_col)]
mean_mat <- mean_mat[module_levels, subtype_levels, drop = FALSE]
stopifnot(!any(is.na(mean_mat)))

z_mat <- t(scale(t(mean_mat)))
z_mat[is.na(z_mat)] <- 0
z_mat[z_mat > 2.5] <- 2.5
z_mat[z_mat < -2.5] <- -2.5

row_family <- factor(module_map$family[match(rownames(z_mat), module_map$module)], levels = family_levels)
column_subtype <- factor(colnames(z_mat), levels = subtype_levels)

row_anno <- rowAnnotation(
  State = row_family,
  col = list(State = c(
    "MES" = "#8C510A",
    "AC" = "#5E3C99",
    "OPC" = "#01665E",
    "NPC" = "#4D4D4D",
    "Cycling" = "#B2182B"
  )),
  show_annotation_name = FALSE,
  width = unit(2.5, "mm"),
  annotation_legend_param = list(
    State = list(
      title = "State",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      grid_height = unit(3, "mm"),
      grid_width = unit(3, "mm")
    )
  )
)

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

heat_col <- circlize::colorRamp2(
  c(-2.5, 0, 2.5),
  c("#4393C3", "white", "#D6604D")
)

ht <- Heatmap(
  z_mat,
  name = "Z-score",
  col = heat_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = row_family,
  column_split = column_subtype,
  row_title = NULL,
  column_title = NULL,
  row_gap = unit(0.9, "mm"),
  column_gap = unit(0.8, "mm"),
  top_annotation = col_anno,
  left_annotation = row_anno,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.25),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
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

pdf(out_pdf, width = 4.8, height = 4.6, useDingbats = FALSE)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(3, 3, 3, 3), "mm")
)
dev.off()

source_tbl <- as.data.frame(as.table(mean_mat), stringsAsFactors = FALSE) %>%
  rename(module = Var1, subtype_k4 = Var2, mean_AMS = Freq) %>%
  left_join(module_map %>% select(module, family), by = "module") %>%
  mutate(
    module = factor(module, levels = module_levels),
    family = factor(family, levels = family_levels),
    subtype_k4 = factor(subtype_k4, levels = subtype_levels),
    z_score = mapply(
      function(m, s) z_mat[as.character(m), as.character(s)],
      module,
      subtype_k4
    )
  ) %>%
  arrange(module, subtype_k4)

readr::write_csv(source_tbl, out_source)
write.csv(z_mat, out_z)

sanity_tbl <- tibble::tibble(
  metric = c(
    "n_cells",
    "n_subtypes",
    "n_modules",
    "subtype_cell_counts",
    "missing_score_values",
    "score_columns_used"
  ),
  value = c(
    nrow(md),
    length(unique(md$subtype_k4)),
    nrow(module_map),
    paste(paste(names(table(md$subtype_k4)), as.integer(table(md$subtype_k4)), sep = "="), collapse = ";"),
    sum(is.na(md[, module_map$metadata_col])),
    paste(module_map$metadata_col, collapse = ";")
  )
)
readr::write_csv(sanity_tbl, out_sanity)

cat("F4 Neftel heatmap completed.\n")
cat("Subtype counts:\n")
print(table(md$subtype_k4))
cat("Mean AMS matrix:\n")
print(round(mean_mat, 3))
cat("Output PDF:", out_pdf, "\n")
cat("Source:", out_source, "\n")
