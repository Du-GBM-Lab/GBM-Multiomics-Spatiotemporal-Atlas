# 05_恶性细胞分亚群与Neftel对照/03e_closer_soft_umap_grid.R
# Display-only UMAP grid: closer and softer variants based on theta=5 + k=4.
# This script does not change Harmony, clustering, subtype assignment, or k.

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

in_path <- "05_恶性细胞分亚群与Neftel对照/outputs/GBM.malignant.subtyped.theta5.k4.qs2"
fig_dir <- "05_恶性细胞分亚群与Neftel对照/figures"
out_dir <- "05_恶性细胞分亚群与Neftel对照/outputs"
tab_dir <- "05_恶性细胞分亚群与Neftel对照/tables"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# Larger n.neighbors strengthens global attraction; moderate min.dist softens edges;
# spread below 1 keeps the layout from looking like separated islands.
CLOSER_GRID <- list(
  list(n.neighbors = 90,  min.dist = 0.18, spread = 0.70, label = "closer1_nn90_md0.18_sp0.7"),
  list(n.neighbors = 120, min.dist = 0.20, spread = 0.70, label = "closer2_nn120_md0.20_sp0.7"),
  list(n.neighbors = 150, min.dist = 0.22, spread = 0.75, label = "closer3_nn150_md0.22_sp0.75"),
  list(n.neighbors = 120, min.dist = 0.30, spread = 0.80, label = "closer4_nn120_md0.30_sp0.8"),
  list(n.neighbors = 180, min.dist = 0.25, spread = 0.70, label = "closer5_nn180_md0.25_sp0.7"),
  list(n.neighbors = 200, min.dist = 0.35, spread = 0.85, label = "closer6_nn200_md0.35_sp0.85")
)

DOWNSTREAM_DIMS <- 1:30
LANCET_K4 <- c(
  "Subtype1" = "#00468B",
  "Subtype2" = "#ED0000",
  "Subtype3" = "#42B540",
  "Subtype4" = "#0099B4"
)

PT_SIZE <- 0.30
PT_ALPHA <- 0.62
LABEL_SIZE <- 5

theme_publication <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.4),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(hjust = 0, margin = margin(b = 4)),
      legend.key = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      plot.margin = margin(6, 6, 6, 6)
    )
}

plot_umap_clean <- function(obj, reduction, title = "", show_label = TRUE) {
  emb <- as.data.frame(Embeddings(obj, reduction))
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  emb$subtype_k4 <- factor(obj$subtype_k4, levels = names(LANCET_K4))
  emb <- emb[sample(nrow(emb)), ]

  centroids <- emb |>
    group_by(subtype_k4) |>
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

  p <- ggplot(emb, aes(UMAP_1, UMAP_2, color = subtype_k4)) +
    geom_point(size = PT_SIZE, alpha = PT_ALPHA, stroke = 0, shape = 16) +
    scale_color_manual(values = LANCET_K4, name = NULL) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) +
    theme_publication()

  if (show_label) {
    p <- p + ggrepel::geom_text_repel(
      data = centroids,
      aes(UMAP_1, UMAP_2, label = subtype_k4),
      color = "black",
      size = LABEL_SIZE,
      fontface = "bold",
      bg.color = "white",
      bg.r = 0.15,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.color = NA,
      inherit.aes = FALSE
    )
  }

  p
}

malig <- qs2::qs_read(in_path)
DefaultAssay(malig) <- "RNA"
stopifnot("harmony" %in% Reductions(malig))
stopifnot("subtype_k4" %in% colnames(malig@meta.data))

cat("Loaded:", ncol(malig), "cells\n")

for (cfg in CLOSER_GRID) {
  rd_name <- paste0("umap_", cfg$label)
  cat("Computing", rd_name, "\n")
  set.seed(42)
  malig <- RunUMAP(
    malig,
    reduction = "harmony",
    dims = DOWNSTREAM_DIMS,
    n.neighbors = cfg$n.neighbors,
    min.dist = cfg$min.dist,
    spread = cfg$spread,
    reduction.name = rd_name,
    reduction.key = paste0(gsub("[^A-Za-z0-9]", "", rd_name), "_"),
    verbose = FALSE
  )
}

plots <- lapply(CLOSER_GRID, function(cfg) {
  rd_name <- paste0("umap_", cfg$label)
  plot_umap_clean(
    malig,
    rd_name,
    title = paste0(
      cfg$label,
      "\nNN=", cfg$n.neighbors,
      " | min.dist=", cfg$min.dist,
      " | spread=", cfg$spread
    ),
    show_label = TRUE
  )
})

ggsave(
  file.path(fig_dir, "03e_closer_soft_umap_grid_Lancet.pdf"),
  wrap_plots(plots, ncol = 2),
  width = 13,
  height = 16,
  device = cairo_pdf
)

for (cfg in CLOSER_GRID) {
  rd_name <- paste0("umap_", cfg$label)
  p <- plot_umap_clean(
    malig,
    rd_name,
    title = paste0("Malignant subtypes | ", cfg$label),
    show_label = TRUE
  )
  ggsave(
    file.path(fig_dir, paste0("03e_final_", cfg$label, "_Lancet.pdf")),
    p,
    width = 6.5,
    height = 5.5,
    device = cairo_pdf
  )
}

param_table <- bind_rows(lapply(CLOSER_GRID, as.data.frame))
write.csv(param_table, file.path(tab_dir, "03e_closer_soft_umap_parameters.csv"), row.names = FALSE)

qs2::qs_save(malig, file.path(out_dir, "GBM.malignant.subtyped.theta5.k4.closer_soft_umap_candidates.qs2"))

cat("\nDone.\n")
cat("Review: figures/03e_closer_soft_umap_grid_Lancet.pdf\n")
cat("Then pick one: 03e_final_closer*_Lancet.pdf\n")
