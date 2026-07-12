suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(ggrastr)
})

project_root <- "<DATA_ROOT>/项目/分型/修稿杠生信"
object_path <- file.path(
  project_root,
  "重新分析/02_scRNA_QC/outputs/GBM.RNA.qc_doubletfinder.filtered.qs2"
)
figure_dir <- file.path(project_root, "图片表格/02_恶性筛选（没有过一遍）/正文图")
source_dir <- file.path(project_root, "图片表格/02_恶性筛选（没有过一遍）/源数据/正图_全细胞注释UMAP")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(figure_dir, "正图_全细胞注释UMAP.pdf")
out_png <- file.path(source_dir, "正图_全细胞注释UMAP_preview.png")
out_labeled_pdf <- file.path(source_dir, "正图_全细胞注释UMAP_labeled_internal.pdf")
out_labeled_png <- file.path(source_dir, "正图_全细胞注释UMAP_labeled_internal.png")
out_coords <- file.path(source_dir, "allcell_annotation_umap_coordinates.csv")
out_counts <- file.path(source_dir, "allcell_annotation_counts.csv")
out_session <- file.path(source_dir, "sessionInfo.txt")

celltype_palette <- c(
  "T cells" = "#1F77B4",
  "Macrophages" = "#D62728",
  "Radial glial" = "#2CA02C",
  "Microglial" = "#9467BD",
  "Mural cells" = "#FF7F0E",
  "OPCs" = "#17BECF",
  "Oligodendrocytes" = "#8C6D31",
  "Ambiguous" = "#BDBDBD",
  "Astrocytes" = "#4E79A7",
  "NK cells" = "#E15759",
  "Monocytes" = "#59A14F",
  "Endothelial" = "#B07AA1",
  "Neurons" = "#EDC948",
  "cDCs" = "#76B7B2",
  "B cells" = "#7B61B3",
  "Ependymal cells" = "#9C755F",
  "pDCs" = "#34495E"
)

obj <- qs_read(object_path)
stopifnot("umap" %in% names(obj@reductions))
stopifnot("anno_ident" %in% colnames(obj@meta.data))

emb <- Embeddings(obj, "umap")
meta <- obj@meta.data
umap_df <- data.frame(
  cell_barcode = rownames(emb),
  UMAP_1 = emb[, 1],
  UMAP_2 = emb[, 2],
  annotation = as.character(meta[rownames(emb), "anno_ident"]),
  sample = as.character(meta[rownames(emb), "Pt_number"]),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

counts_df <- umap_df |>
  count(annotation, name = "n_cells") |>
  arrange(desc(n_cells)) |>
  mutate(percent = n_cells / sum(n_cells) * 100)

annotation_levels <- counts_df$annotation
umap_df$annotation <- factor(umap_df$annotation, levels = annotation_levels)

missing_palette <- setdiff(annotation_levels, names(celltype_palette))
if (length(missing_palette) > 0) {
  stop("Missing palette values for: ", paste(missing_palette, collapse = ", "))
}

write.csv(umap_df, out_coords, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(counts_df, out_counts, row.names = FALSE, fileEncoding = "UTF-8")

label_df <- umap_df |>
  group_by(annotation) |>
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

plot_df <- umap_df |>
  left_join(counts_df[, c("annotation", "n_cells")], by = "annotation") |>
  arrange(desc(n_cells))

base_umap <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = annotation)) +
  rasterise(
    geom_point(size = 0.045, alpha = 0.82, stroke = 0),
    dpi = 600
  ) +
  scale_color_manual(values = celltype_palette[annotation_levels], drop = FALSE) +
  guides(color = guide_legend(
    title = "Cell type",
    override.aes = list(size = 2.2, alpha = 1),
    ncol = 1
  )) +
  coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 8) +
  theme(
    axis.line = element_line(linewidth = 0.3, color = "#303030"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 8.5, color = "#202020"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 7.5),
    legend.text = element_text(size = 6.3),
    legend.key.height = unit(0.24, "cm"),
    legend.key.width = unit(0.24, "cm"),
    legend.spacing.y = unit(0.02, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(2, 2, 2, 2)
  )

p <- base_umap

p_labeled <- base_umap +
  geom_text_repel(
    data = label_df,
    aes(x = UMAP_1, y = UMAP_2, label = annotation),
    size = 2.35,
    color = "#111111",
    min.segment.length = Inf,
    segment.size = 0,
    max.overlaps = Inf,
    seed = 20260612,
    inherit.aes = FALSE
  ) +
  theme(legend.position = "none")

ggsave(out_pdf, p, width = 5.2, height = 4.6, units = "in", device = cairo_pdf, bg = "white")
ggsave(out_png, p, width = 5.2, height = 4.6, units = "in", dpi = 300, bg = "white")
ggsave(out_labeled_pdf, p_labeled, width = 4.8, height = 4.8, units = "in", device = cairo_pdf, bg = "white")
ggsave(out_labeled_png, p_labeled, width = 4.8, height = 4.8, units = "in", dpi = 300, bg = "white")

capture.output(sessionInfo(), file = out_session)

message("Saved: ", out_pdf)
message("Saved source data: ", source_dir)
