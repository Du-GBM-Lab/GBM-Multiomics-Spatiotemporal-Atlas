# 05_恶性细胞分亚群与Neftel对照/03b_umap_publication_figures.R
# 纯展示层: 重投 UMAP(多参数候选) + publication palettes + 重画
# 不动 PCA / Harmony / cluster / subtype 归属
# k=3 和 k=4 都画

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(42)

# ---- 路径 ----
in_path <- "05_恶性细胞分亚群与Neftel对照/outputs/GBM.malignant.reclustered.qs2"
evidence <- "05_恶性细胞分亚群与Neftel对照/outputs/03_subtype_evidence_bundle.rds"
out_dir <- "05_恶性细胞分亚群与Neftel对照/outputs"
fig_dir <- "05_恶性细胞分亚群与Neftel对照/figures"

# ============================================================
# 参数: publication display layer
# ============================================================
UMAP_GRID <- list(
  list(n.neighbors = 30, min.dist = 0.3, spread = 1.0, label = "nn30_md0.3"),
  list(n.neighbors = 30, min.dist = 0.5, spread = 1.0, label = "nn30_md0.5"),
  list(n.neighbors = 50, min.dist = 0.5, spread = 1.0, label = "nn50_md0.5"),
  list(n.neighbors = 50, min.dist = 0.8, spread = 1.0, label = "nn50_md0.8")
)
UMAP_DIMS <- 1:30

PALETTE_NPG <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")
PALETTE_JCO <- c("#0073C2", "#EFC000", "#868686", "#CD534C", "#7AA6DC", "#003C67")
PALETTE_LANCET <- c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91")

K_CANDIDATES <- c(3, 4)

PT_SIZE <- 0.4
PT_ALPHA <- 0.85
LABEL_SIZE <- 5
FIG_W <- 6.5
FIG_H <- 5.5

# ============================================================
# 加载
# ============================================================
malig <- qs2::qs_read(in_path)
ev <- readRDS(evidence)
DefaultAssay(malig) <- "RNA"
cat("Loaded:", ncol(malig), "cells\n")

cluster_col_ref <- paste0("harmony_res.", ev$ref_resolution)
for (k in K_CANDIDATES) {
  sub_map <- ev$k_assignments %>% filter(k == !!k)
  malig[[paste0("subtype_k", k)]] <- factor(
    sub_map$subtype[match(as.character(malig[[cluster_col_ref]][, 1]), sub_map$cluster)],
    levels = paste0("Subtype", seq_len(k))
  )
}

# ============================================================
# 重投 UMAP: 只基于已有 harmony reduction
# ============================================================
for (cfg in UMAP_GRID) {
  rd_name <- paste0("umap_", cfg$label)
  if (!(rd_name %in% Reductions(malig))) {
    cat("  Computing UMAP:", cfg$label, "\n")
    set.seed(42)
    malig <- RunUMAP(
      malig,
      reduction = "harmony",
      dims = UMAP_DIMS,
      n.neighbors = cfg$n.neighbors,
      min.dist = cfg$min.dist,
      spread = cfg$spread,
      reduction.name = rd_name,
      reduction.key = paste0(gsub("[^A-Za-z0-9]", "", rd_name), "_"),
      verbose = FALSE
    )
  }
}

# ============================================================
# publication theme + helper
# ============================================================
theme_publication <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = "") +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.4, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      axis.text = element_text(size = base_size - 2, color = "black"),
      axis.title = element_text(size = base_size, color = "black"),
      plot.title = element_text(size = base_size + 1, hjust = 0, face = "plain", margin = margin(b = 4)),
      legend.position = "right",
      legend.key = element_blank(),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      legend.key.size = unit(0.4, "cm"),
      plot.margin = margin(6, 6, 6, 6)
    )
}

plot_umap_clean <- function(obj, reduction, group_col, palette,
                            title = "", show_label = TRUE,
                            pt_size = PT_SIZE, pt_alpha = PT_ALPHA) {
  emb <- as.data.frame(Embeddings(obj, reduction))
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  emb$grp <- obj[[group_col]][, 1]
  emb <- emb[sample(nrow(emb)), ]

  centroids <- emb %>%
    group_by(grp) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

  n_grp <- length(unique(emb$grp))
  cols <- palette[seq_len(n_grp)]
  names(cols) <- levels(emb$grp)

  p <- ggplot(emb, aes(UMAP_1, UMAP_2, color = grp)) +
    geom_point(size = pt_size, alpha = pt_alpha, stroke = 0, shape = 16) +
    scale_color_manual(values = cols, name = NULL) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) +
    theme_publication()

  if (show_label) {
    p <- p + ggrepel::geom_text_repel(
      data = centroids,
      aes(UMAP_1, UMAP_2, label = grp),
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

# ============================================================
# Step A: UMAP 参数对比
# ============================================================
palettes_list <- list(NPG = PALETTE_NPG, JCO = PALETTE_JCO, Lancet = PALETTE_LANCET)

for (k in K_CANDIDATES) {
  for (pal_name in names(palettes_list)) {
    pal <- palettes_list[[pal_name]]
    plots <- lapply(UMAP_GRID, function(cfg) {
      rd_name <- paste0("umap_", cfg$label)
      plot_umap_clean(
        malig,
        rd_name,
        paste0("subtype_k", k),
        pal,
        title = paste0("nn=", cfg$n.neighbors, " | min.dist=", cfg$min.dist),
        show_label = TRUE
      )
    })
    fname <- file.path(fig_dir, sprintf("03b_umap_grid_k%d_%s.pdf", k, pal_name))
    ggsave(
      fname,
      wrap_plots(plots, ncol = 2) +
        plot_annotation(
          title = sprintf("UMAP grid | k=%d | palette: %s", k, pal_name),
          theme = theme(plot.title = element_text(face = "bold"))
        ),
      width = 13,
      height = 11,
      device = cairo_pdf
    )
  }
}

# ============================================================
# Step B: 配色对比
# ============================================================
fixed_rd <- "umap_nn50_md0.5"
for (k in K_CANDIDATES) {
  plots <- lapply(names(palettes_list), function(pal_name) {
    plot_umap_clean(
      malig,
      fixed_rd,
      paste0("subtype_k", k),
      palettes_list[[pal_name]],
      title = pal_name,
      show_label = TRUE
    )
  })
  fname <- file.path(fig_dir, sprintf("03b_palette_compare_k%d_%s.pdf", k, fixed_rd))
  ggsave(
    fname,
    wrap_plots(plots, ncol = 3) +
      plot_annotation(
        title = sprintf("Palette comparison | k=%d | %s", k, fixed_rd),
        theme = theme(plot.title = element_text(face = "bold"))
      ),
    width = 16,
    height = 5.5,
    device = cairo_pdf
  )
}

# ============================================================
# Step C: 单图清版
# ============================================================
for (k in K_CANDIDATES) {
  for (pal_name in names(palettes_list)) {
    p <- plot_umap_clean(
      malig,
      fixed_rd,
      paste0("subtype_k", k),
      palettes_list[[pal_name]],
      title = sprintf("Malignant subtypes (k = %d)", k),
      show_label = TRUE
    )
    fname <- file.path(fig_dir, sprintf("03b_final_k%d_%s.pdf", k, pal_name))
    ggsave(fname, p, width = FIG_W, height = FIG_H, device = cairo_pdf)
  }
}

# ============================================================
# Step D: 诊断叠加
# ============================================================
patient_palette_24 <- colorRampPalette(
  c(
    "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
    "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
    "#0073C2", "#EFC000", "#868686", "#CD534C", "#7AA6DC",
    "#003C67", "#8F7700", "#3B3B3B", "#A73030", "#4A6990",
    "#00468B", "#ED0000", "#42B540", "#0099B4"
  )
)(24)

phase_palette <- c("G1" = "#BBBBBB", "S" = "#3C5488", "G2M" = "#E64B35")

p_pt <- plot_umap_clean(
  malig,
  fixed_rd,
  "Pt_number",
  patient_palette_24,
  title = "By Patient",
  show_label = FALSE,
  pt_size = 0.3,
  pt_alpha = 0.7
) +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 2)))

p_ph <- plot_umap_clean(
  malig,
  fixed_rd,
  "Phase",
  phase_palette,
  title = "By Cell-Cycle Phase",
  show_label = FALSE,
  pt_size = 0.3,
  pt_alpha = 0.7
)

ggsave(file.path(fig_dir, "03b_diagnostic_overlay.pdf"), p_pt | p_ph, width = 13, height = 5.5, device = cairo_pdf)

pt_highlight <- ifelse(
  malig$Pt_number == "Pt8",
  "Pt8",
  ifelse(malig$Pt_number == "Pt18", "Pt18", "Other")
)
malig$pt_highlight <- factor(pt_highlight, levels = c("Other", "Pt18", "Pt8"))
highlight_palette <- c("Other" = "grey85", "Pt18" = "#3C5488", "Pt8" = "#E64B35")

p_hl <- plot_umap_clean(
  malig,
  fixed_rd,
  "pt_highlight",
  highlight_palette,
  title = "Pt8 / Pt18 highlight",
  show_label = FALSE,
  pt_size = 0.35,
  pt_alpha = 0.8
)
ggsave(file.path(fig_dir, "03b_Pt8_Pt18_on_new_umap.pdf"), p_hl, width = FIG_W, height = FIG_H, device = cairo_pdf)

# ============================================================
# 保存
# ============================================================
qs2::qs_save(malig, file.path(out_dir, "GBM.malignant.reclustered.qs2"))

cat("\nDone.\n")
cat("\n看图顺序:\n")
cat("  1. 03b_umap_grid_k{3,4}_{NPG,JCO,Lancet}.pdf\n")
cat("  2. 03b_palette_compare_k{3,4}_umap_nn50_md0.5.pdf\n")
cat("  3. 03b_final_k{3,4}_{NPG,JCO,Lancet}.pdf\n")
cat("  4. 03b_diagnostic_overlay.pdf\n")
cat("  5. 03b_Pt8_Pt18_on_new_umap.pdf\n")
