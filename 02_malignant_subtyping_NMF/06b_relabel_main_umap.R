# 05_恶性细胞分亚群与Neftel对照/06b_relabel_main_umap.R
# Relabel final closer4 UMAP with manually placed subtype labels.
# This fixes label placement only; it does not change cells, clusters, or metadata.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(Seurat)
  library(qs2)
  library(ggplot2)
  library(dplyr)
})

set.seed(42)

proj <- "05_恶性细胞分亚群与Neftel对照"
in_obj <- file.path(proj, "outputs", "GBM.malignant.subtyped.named.qs2")
fig_dir <- file.path(proj, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

named_levels <- c("NPC-Cycling", "OPC-like", "MES-Perivascular", "MES-Inflammatory")
named_palette <- c(
  "NPC-Cycling" = "#00468B",
  "OPC-like" = "#ED0000",
  "MES-Perivascular" = "#42B540",
  "MES-Inflammatory" = "#0099B4"
)

obj <- qs2::qs_read(in_obj)
stopifnot("umap_closer4" %in% Reductions(obj))
stopifnot("subtype_named" %in% colnames(obj@meta.data))

emb <- as.data.frame(Embeddings(obj, "umap_closer4"))
colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
emb$subtype <- factor(obj$subtype_named, levels = named_levels)
emb <- emb[sample(nrow(emb)), ]

label_pos <- data.frame(
  subtype = factor(named_levels, levels = named_levels),
  x = c(2.5, -14.0, 0.0, 4.5),
  y = c(-1.0, -1.0, -10.0, 9.5),
  stringsAsFactors = FALSE
)

p <- ggplot(emb, aes(UMAP_1, UMAP_2, color = subtype)) +
  geom_point(size = 0.30, alpha = 0.62, stroke = 0, shape = 16) +
  geom_text(
    data = label_pos,
    aes(x = x, y = y, label = subtype, color = subtype),
    size = 5,
    fontface = "bold",
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = named_palette, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    title = "Malignant subtypes (n=28,213 cells, 24 patients)",
    x = "UMAP 1",
    y = "UMAP 2",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.size = unit(0.42, "cm"),
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.4)
  )

ggsave(
  file.path(fig_dir, "06b_final_main_UMAP_relabeled.pdf"),
  p,
  width = 7.5,
  height = 5.5,
  device = cairo_pdf
)

ggsave(
  file.path(fig_dir, "06b_final_main_UMAP_relabeled.png"),
  p,
  width = 7.5,
  height = 5.5,
  dpi = 300
)

cat("Done.\n")
cat("PDF:", file.path(fig_dir, "06b_final_main_UMAP_relabeled.pdf"), "\n")
cat("PNG:", file.path(fig_dir, "06b_final_main_UMAP_relabeled.png"), "\n")
