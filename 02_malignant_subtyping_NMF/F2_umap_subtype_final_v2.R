# F2_umap_subtype_final_v2.R
# Purpose: Polished final-labeled malignant subtype UMAP.
# Style: bottom legend with cell counts, corner axis arrows, and in-plot subtype
# labels with text-halo background. No rectangular label boxes.
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
  "Subtype1" = "#0072B5",
  "Subtype2" = "#E18727",
  "Subtype3" = "#20854E",
  "Subtype4" = "#BC3C29"
)

label_obj <- qs2::qs_read(label_obj_path)
coord_obj <- qs2::qs_read(coord_obj_path)

stopifnot("subtype_k4" %in% colnames(label_obj@meta.data))
stopifnot("umap_closer4" %in% names(coord_obj@reductions))
stopifnot(setequal(colnames(label_obj), colnames(coord_obj)))

umap_df <- as.data.frame(Embeddings(coord_obj, "umap_closer4"))
stopifnot(ncol(umap_df) >= 2)
umap_df <- umap_df[, 1:2, drop = FALSE]
colnames(umap_df) <- c("UMAP_1", "UMAP_2")

label_md <- label_obj@meta.data[rownames(umap_df), , drop = FALSE]
stopifnot(identical(rownames(label_md), rownames(umap_df)))

umap_df$subtype <- factor(label_md$subtype_k4, levels = paste0("Subtype", 1:4))
umap_df$cell_id <- rownames(umap_df)
umap_df$reduction <- "umap_closer4"

stopifnot(nrow(umap_df) == 28213)
stopifnot(!any(is.na(umap_df$subtype)))
stopifnot(setequal(as.character(unique(umap_df$subtype)), names(subtype_colors)))
stopifnot(!anyNA(umap_df$UMAP_1), !anyNA(umap_df$UMAP_2))

subtype_size <- table(umap_df$subtype)
umap_df <- umap_df %>%
  mutate(subtype_n = as.numeric(subtype_size[as.character(subtype)])) %>%
  arrange(desc(subtype_n))

legend_labels <- sprintf(
  "%s  (n=%s)",
  names(subtype_size),
  format(as.numeric(subtype_size), big.mark = ",")
)
names(legend_labels) <- names(subtype_size)

x_range <- range(umap_df$UMAP_1)
y_range <- range(umap_df$UMAP_2)
arrow_x_start <- x_range[1] - diff(x_range) * 0.02
arrow_x_end <- x_range[1] + diff(x_range) * 0.15
arrow_y_start <- y_range[1] - diff(y_range) * 0.02
arrow_y_end <- y_range[1] + diff(y_range) * 0.15
arrow_pad <- diff(y_range) * 0.04

label_pos <- data.frame(
  subtype = factor(names(subtype_colors), levels = names(subtype_colors)),
  label = names(subtype_colors),
  UMAP_1 = c(2.5, -14.0, 0.0, 4.5),
  UMAP_2 = c(-1.0, -1.0, -10.0, 9.5),
  stringsAsFactors = FALSE
)

halo_dx <- diff(x_range) * 0.004
halo_dy <- diff(y_range) * 0.004
halo_offsets <- expand.grid(dx = c(-halo_dx, 0, halo_dx), dy = c(-halo_dy, 0, halo_dy))
halo_offsets <- halo_offsets[!(halo_offsets$dx == 0 & halo_offsets$dy == 0), ]
label_halo <- merge(label_pos, halo_offsets)
label_halo$UMAP_1 <- label_halo$UMAP_1 + label_halo$dx
label_halo$UMAP_2 <- label_halo$UMAP_2 + label_halo$dy

cat("Subtype cell counts:\n")
print(subtype_size)
cat("Coordinate reduction: umap_closer4\n")
cat("Label object:", label_obj_path, "\n")
cat("Coordinate object:", coord_obj_path, "\n")

p <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color = subtype), size = 0.4, alpha = 0.85, stroke = 0, shape = 16) +
  annotate(
    "segment",
    x = arrow_x_start,
    xend = arrow_x_end,
    y = arrow_y_start - arrow_pad,
    yend = arrow_y_start - arrow_pad,
    arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
    linewidth = 0.4,
    color = "black"
  ) +
  annotate(
    "segment",
    x = arrow_x_start - arrow_pad,
    xend = arrow_x_start - arrow_pad,
    y = arrow_y_start,
    yend = arrow_y_end,
    arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
    linewidth = 0.4,
    color = "black"
  ) +
  annotate(
    "text",
    x = (arrow_x_start + arrow_x_end) / 2,
    y = arrow_y_start - arrow_pad * 2.2,
    label = "UMAP 1",
    size = 2.5,
    hjust = 0.5,
    color = "black"
  ) +
  annotate(
    "text",
    x = arrow_x_start - arrow_pad * 2.2,
    y = (arrow_y_start + arrow_y_end) / 2,
    label = "UMAP 2",
    size = 2.5,
    angle = 90,
    hjust = 0.5,
    color = "black"
  ) +
  geom_text(
    data = label_halo,
    aes(UMAP_1, UMAP_2, label = label),
    color = "white",
    alpha = 0.9,
    size = 3.5,
    fontface = "bold",
    lineheight = 0.9,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_text(
    data = label_pos,
    aes(UMAP_1, UMAP_2, label = label, color = subtype),
    size = 3.5,
    fontface = "bold",
    lineheight = 0.9,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_color_manual(values = subtype_colors, labels = legend_labels, name = NULL, drop = FALSE) +
  guides(color = guide_legend(
    override.aes = list(size = 3, alpha = 1),
    nrow = 2,
    byrow = TRUE
  )) +
  coord_fixed() +
  theme_void(base_size = 8) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7, color = "black"),
    legend.key.size = unit(3, "mm"),
    legend.spacing.x = unit(3, "mm"),
    legend.margin = margin(t = 8, b = 0),
    plot.margin = margin(8, 8, 8, 8),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(out_pdf, p, width = 5, height = 5.5, device = cairo_pdf)

write.csv(
  umap_df %>%
    select(cell_id, UMAP_1, UMAP_2, subtype, reduction),
  source_csv,
  row.names = FALSE
)

cat("Output PDF:", out_pdf, "\n")
cat("Source data:", source_csv, "\n")
cat("Sanity checks passed.\n")
