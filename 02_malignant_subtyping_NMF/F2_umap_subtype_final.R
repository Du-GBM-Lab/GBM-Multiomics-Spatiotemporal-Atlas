# F2_umap_subtype_final.R
# Purpose: Draw final-labeled malignant subtype UMAP.
#
# Input:
#   outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2
#     - final subtype labels
#   outputs/GBM.malignant.subtyped.umap_candidates.qs2
#     - malignant-only display UMAP reduction `umap_closer4`
#
# Output:
#   figures/F2_umap_subtype_final.pdf
#   figures/source_data/F2_umap_coordinates.csv

.libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(ggplot2)
  library(dplyr)
  library(grid)
  library(scales)
})

base_dir <- "05_恶性细胞分亚群与Neftel对照"
fig_dir <- file.path(base_dir, "figures")
source_dir <- file.path(fig_dir, "source_data")
dir.create(source_dir, showWarnings = FALSE, recursive = TRUE)

label_obj_path <- file.path(base_dir, "outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2")
coord_obj_path <- file.path(base_dir, "outputs/GBM.malignant.subtyped.umap_candidates.qs2")
out_pdf <- file.path(fig_dir, "F2_umap_subtype_final.pdf")
source_csv <- file.path(source_dir, "F2_umap_coordinates.csv")

subtype_colors <- c(
  "Proliferative-NPC subtype" = "#0072B5",
  "OPC-Myelination subtype" = "#E18727",
  "MES-ECM subtype" = "#20854E",
  "MES-Antigen-presenting subtype" = "#BC3C29"
)

label_obj <- qs2::qs_read(label_obj_path)
coord_obj <- qs2::qs_read(coord_obj_path)
stopifnot("subtype_label_final" %in% colnames(label_obj@meta.data))
stopifnot("subtype_k4" %in% colnames(label_obj@meta.data))
stopifnot("umap_closer4" %in% names(coord_obj@reductions))
stopifnot(setequal(colnames(label_obj), colnames(coord_obj)))

umap_df <- as.data.frame(Embeddings(coord_obj, "umap_closer4"))
stopifnot(ncol(umap_df) >= 2)
umap_df <- umap_df[, 1:2, drop = FALSE]
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
label_md <- label_obj@meta.data[rownames(umap_df), , drop = FALSE]
stopifnot(identical(rownames(label_md), rownames(umap_df)))
umap_df$subtype <- factor(label_md$subtype_label_final, levels = names(subtype_colors))
umap_df$subtype_id <- factor(label_md$subtype_k4, levels = paste0("Subtype", 1:4))
umap_df$cell_id <- rownames(umap_df)
umap_df$reduction <- "umap_closer4"

stopifnot(nrow(umap_df) == 28213)
stopifnot(!any(is.na(umap_df$subtype)))
stopifnot(setequal(as.character(unique(umap_df$subtype)), names(subtype_colors)))
stopifnot(!anyNA(umap_df$UMAP_1), !anyNA(umap_df$UMAP_2))

cat("Subtype cell counts:\n")
print(table(umap_df$subtype))
cat("Coordinate reduction: umap_closer4\n")
cat("Label object:", label_obj_path, "\n")
cat("Coordinate object:", coord_obj_path, "\n")

label_pos <- data.frame(
  subtype = factor(names(subtype_colors), levels = names(subtype_colors)),
  label = c(
    "Proliferative-NPC\nsubtype",
    "OPC-Myelination\nsubtype",
    "MES-ECM\nsubtype",
    "MES-Antigen-presenting\nsubtype"
  ),
  UMAP_1 = c(2.5, -14.0, 0.0, 4.5),
  UMAP_2 = c(-1.0, -1.0, -10.0, 9.5),
  stringsAsFactors = FALSE
)

p <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = subtype)) +
  geom_point(size = 0.25, alpha = 0.62, stroke = 0, shape = 16) +
  geom_label(
    data = label_pos,
    aes(UMAP_1, UMAP_2, label = label, fill = subtype),
    color = "white",
    fontface = "bold",
    size = 3.0,
    linewidth = 0.18,
    label.padding = unit(0.18, "lines"),
    label.r = unit(0.08, "lines"),
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = subtype_colors, name = NULL, drop = FALSE) +
  scale_fill_manual(values = alpha(subtype_colors, 0.88), guide = "none") +
  guides(color = "none") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  coord_fixed() +
  theme_classic(base_size = 8) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_line(linewidth = 0.3, color = "black"),
    axis.title = element_text(size = 8, color = "black"),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(3, 3, 3, 3, "pt")
  )

ggsave(out_pdf, p, width = 4.8, height = 4.8, device = cairo_pdf)

write.csv(umap_df, source_csv, row.names = FALSE)

cat("Output PDF:", out_pdf, "\n")
cat("Source data:", source_csv, "\n")
cat("Sanity checks passed.\n")
