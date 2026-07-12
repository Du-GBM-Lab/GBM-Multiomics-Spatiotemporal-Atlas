# 05_恶性细胞分亚群与Neftel对照/03d_soft_umap_parameter_grid.R
# display-only: softer / closer UMAP candidates for theta=5 + k=4
# Goal: bring separated islands closer and make edges less sharp.

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

SOFT_GRID <- list(
  list(n.neighbors = 20, min.dist = 0.03, spread = 0.50, label = "soft1_nn20_md0.03_sp0.5"),
  list(n.neighbors = 30, min.dist = 0.05, spread = 0.60, label = "soft2_nn30_md0.05_sp0.6"),
  list(n.neighbors = 40, min.dist = 0.08, spread = 0.70, label = "soft3_nn40_md0.08_sp0.7"),
  list(n.neighbors = 50, min.dist = 0.10, spread = 0.80, label = "soft4_nn50_md0.10_sp0.8"),
  list(n.neighbors = 75, min.dist = 0.15, spread = 1.00, label = "soft5_nn75_md0.15_sp1.0")
)

DOWNSTREAM_DIMS <- 1:30
LANCET_K4 <- c("Subtype1" = "#00468B", "Subtype2" = "#ED0000", "Subtype3" = "#42B540", "Subtype4" = "#0099B4")

PT_SIZE <- 0.32
PT_ALPHA <- 0.68
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

  centroids <- emb %>%
    group_by(subtype_k4) %>%
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

for (cfg in SOFT_GRID) {
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

plots <- lapply(SOFT_GRID, function(cfg) {
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
  file.path(fig_dir, "03d_soft_umap_grid_Lancet.pdf"),
  wrap_plots(plots, ncol = 2),
  width = 13,
  height = 14,
  device = cairo_pdf
)

for (cfg in SOFT_GRID) {
  rd_name <- paste0("umap_", cfg$label)
  p <- plot_umap_clean(
    malig,
    rd_name,
    title = paste0("Malignant subtypes | ", cfg$label),
    show_label = TRUE
  )
  ggsave(
    file.path(fig_dir, paste0("03d_final_", cfg$label, "_Lancet.pdf")),
    p,
    width = 6.5,
    height = 5.5,
    device = cairo_pdf
  )
}

# sanity overlays on the middle candidate.
fixed_rd <- "umap_soft3_nn40_md0.08_sp0.7"
emb <- as.data.frame(Embeddings(malig, fixed_rd))
colnames(emb) <- c("UMAP_1", "UMAP_2")
emb$pt_hl <- factor(
  ifelse(malig$Pt_number == "Pt8", "Pt8", ifelse(malig$Pt_number == "Pt18", "Pt18", "Other")),
  levels = c("Other", "Pt18", "Pt8")
)
emb$Phase <- factor(malig$Phase, levels = c("G1", "S", "G2M"))

p_pt <- ggplot(emb[sample(nrow(emb)), ], aes(UMAP_1, UMAP_2, color = pt_hl)) +
  geom_point(size = PT_SIZE, alpha = PT_ALPHA, stroke = 0) +
  scale_color_manual(values = c("Other" = "grey85", "Pt18" = "#3C5488", "Pt8" = "#E64B35"), name = NULL) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Pt8 / Pt18 | soft3") +
  theme_publication()

p_phase <- ggplot(emb[sample(nrow(emb)), ], aes(UMAP_1, UMAP_2, color = Phase)) +
  geom_point(size = PT_SIZE, alpha = PT_ALPHA, stroke = 0) +
  scale_color_manual(values = c("G1" = "grey75", "S" = "#3C5488", "G2M" = "#E64B35"), name = NULL) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Cell cycle | soft3") +
  theme_publication()

ggsave(
  file.path(fig_dir, "03d_diagnostic_soft3.pdf"),
  p_pt | p_phase,
  width = 13,
  height = 5.5,
  device = cairo_pdf
)

param_table <- bind_rows(lapply(SOFT_GRID, as.data.frame))
write.csv(param_table, file.path(tab_dir, "03d_soft_umap_parameters.csv"), row.names = FALSE)

qs2::qs_save(malig, file.path(out_dir, "GBM.malignant.subtyped.theta5.k4.soft_umap_candidates.qs2"))

cat("\nDone.\n")
cat("Review: figures/03d_soft_umap_grid_Lancet.pdf\n")
cat("Then pick one: 03d_final_soft*_Lancet.pdf\n")
