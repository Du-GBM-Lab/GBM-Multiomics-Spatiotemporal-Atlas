# 05_恶性细胞分亚群与Neftel对照/03c_finalize_theta5_k4.R
# 收尾: theta=5 + k=4 + ultra1 UMAP + Lancet 配色
# 输出最终候选对象 + sensitivity 对比表

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(cluster)
})

set.seed(42)

# ---- 路径 ----
in_path <- "05_恶性细胞分亚群与Neftel对照/outputs/GBM.malignant.reclustered.theta5.qs2"
out_dir <- "05_恶性细胞分亚群与Neftel对照/outputs"
tab_dir <- "05_恶性细胞分亚群与Neftel对照/tables"
fig_dir <- "05_恶性细胞分亚群与Neftel对照/figures"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 参数 ----
TARGET_K <- 4
DOWNSTREAM_DIMS <- 30
DEFAULT_RES <- 0.8

ULTRA_TIGHT_GRID <- list(
  list(n.neighbors = 8, min.dist = 0.001, spread = 0.3, label = "ultra1")
)

PALETTES <- list(
  Lancet = c("#00468B", "#ED0000", "#42B540", "#0099B4")
)

TOP_N_MARKER <- 50
COR_METHOD <- "spearman"
HCLUST_METHOD <- "complete"
MIN_LOGFC <- 0.25
MIN_PCT <- 0.25

# ============================================================
# Helper
# ============================================================
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

plot_umap_clean <- function(emb_df, group_col, palette,
                            title = "", show_label = TRUE,
                            pt_size = 0.4, pt_alpha = 0.85) {
  emb <- emb_df
  emb$grp <- emb[[group_col]]
  emb <- emb[sample(nrow(emb)), ]

  centroids <- emb %>%
    group_by(grp) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

  groups <- if (is.factor(emb$grp)) levels(emb$grp) else as.character(unique(emb$grp))
  if (!is.null(names(palette)) && all(groups %in% names(palette))) {
    cols <- palette[groups]
  } else {
    cols <- palette[seq_along(groups)]
    names(cols) <- groups
  }

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
      size = 5,
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
# 1. 加载 theta=5 对象
# ============================================================
malig <- qs2::qs_read(in_path)
DefaultAssay(malig) <- "RNA"
cat("Loaded theta=5 object:", ncol(malig), "cells\n")

# ============================================================
# 2. cluster signature -> cor -> cutree(k=4), 写回 subtype_k4
# ============================================================
cluster_col_ref <- paste0("harmony_res.", DEFAULT_RES)
marker_features <- VariableFeatures(malig)
if (length(marker_features) == 0) {
  marker_features <- rownames(malig)
}

Idents(malig) <- cluster_col_ref
cl_sizes <- table(Idents(malig))
small_cl <- names(cl_sizes)[cl_sizes < 10]
if (length(small_cl) > 0) {
  cat("Skipping small clusters:", paste(small_cl, collapse = ", "), "\n")
  obj_for_marker <- subset(malig, idents = setdiff(names(cl_sizes), small_cl))
} else {
  obj_for_marker <- malig
}

cat("Computing markers for cluster signatures on ", length(marker_features), " features...\n", sep = "")
markers <- FindAllMarkers(
  obj_for_marker,
  features = marker_features,
  only.pos = TRUE,
  min.pct = MIN_PCT,
  logfc.threshold = MIN_LOGFC,
  verbose = FALSE
) %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = TOP_N_MARKER, with_ties = FALSE) %>%
  ungroup()

union_genes <- unique(markers$gene)
sig_mat <- AverageExpression(
  obj_for_marker,
  assays = "RNA",
  features = union_genes,
  slot = "data",
  verbose = FALSE
)$RNA
colnames(sig_mat) <- sub("^g(?=\\d)", "", colnames(sig_mat), perl = TRUE)
sig_mat <- as.matrix(sig_mat)

cor_mat <- cor(sig_mat, method = COR_METHOD)
hc_ref <- hclust(as.dist(1 - cor_mat), method = HCLUST_METHOD)

ass <- cutree(hc_ref, k = TARGET_K)
sub_map <- data.frame(cluster = names(ass), subtype = paste0("Subtype", ass))

malig$subtype_k4 <- factor(
  sub_map$subtype[match(as.character(malig[[cluster_col_ref]][, 1]), sub_map$cluster)],
  levels = paste0("Subtype", seq_len(TARGET_K))
)

write.csv(sub_map, file.path(tab_dir, "03c_cluster_to_subtype_k4_theta5.csv"), row.names = FALSE)
write.csv(markers, file.path(tab_dir, "03c_top50_markers_cluster_theta5_res0.8.csv"), row.names = FALSE)

subtype_count <- malig@meta.data %>% count(subtype_k4, name = "n_cells")
write.csv(subtype_count, file.path(tab_dir, "03c_subtype_k4_cell_counts.csv"), row.names = FALSE)
cat("\nSubtype k=4 cell counts:\n")
print(subtype_count)

# ============================================================
# 3. 重投 3 套 ultra-tight UMAP
# ============================================================
cat("\nRunning ultra-tight UMAPs...\n")
for (cfg in ULTRA_TIGHT_GRID) {
  rd_name <- paste0("umap_", cfg$label)
  set.seed(42)
  malig <- RunUMAP(
    malig,
    reduction = "harmony",
    dims = 1:DOWNSTREAM_DIMS,
    n.neighbors = cfg$n.neighbors,
    min.dist = cfg$min.dist,
    spread = cfg$spread,
    reduction.name = rd_name,
    reduction.key = paste0(gsub("[^A-Za-z0-9]", "", rd_name), "_"),
    verbose = FALSE
  )
  cat("  ", cfg$label, " (nn=", cfg$n.neighbors, ", md=", cfg$min.dist, ", spread=", cfg$spread, ") done\n", sep = "")
}

# ============================================================
# 4. 出图: final Lancet UMAP
# ============================================================
for (cfg in ULTRA_TIGHT_GRID) {
  rd_name <- paste0("umap_", cfg$label)
  emb <- as.data.frame(Embeddings(malig, rd_name))
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  emb$subtype_k4 <- malig$subtype_k4

  for (pn in names(PALETTES)) {
    p <- plot_umap_clean(
      emb,
      "subtype_k4",
      PALETTES[[pn]],
      title = sprintf("Malignant subtypes (k=4) | %s | %s", cfg$label, pn),
      show_label = TRUE
    )
    ggsave(file.path(fig_dir, sprintf("03c_final_%s_%s.pdf", cfg$label, pn)), p, width = 6.5, height = 5.5, device = cairo_pdf)
  }
}

# ============================================================
# 5. Patient + Phase on ultra1 sanity check
# ============================================================
emb_ultra1 <- as.data.frame(Embeddings(malig, "umap_ultra1"))
colnames(emb_ultra1) <- c("UMAP_1", "UMAP_2")

emb_ultra1$pt_hl <- factor(
  ifelse(malig$Pt_number == "Pt8", "Pt8", ifelse(malig$Pt_number == "Pt18", "Pt18", "Other")),
  levels = c("Other", "Pt18", "Pt8")
)

p_pt <- plot_umap_clean(
  emb_ultra1,
  "pt_hl",
  c("Other" = "grey85", "Pt18" = "#3C5488", "Pt8" = "#E64B35"),
  title = "Pt8 / Pt18 on ultra1 UMAP",
  show_label = FALSE
)

emb_ultra1$Phase <- factor(malig$Phase, levels = c("G1", "S", "G2M"))
p_ph <- plot_umap_clean(
  emb_ultra1,
  "Phase",
  c("G1" = "grey75", "S" = "#3C5488", "G2M" = "#E64B35"),
  title = "Cell cycle phase",
  show_label = FALSE
)

ggsave(file.path(fig_dir, "03c_diagnostic_ultra1.pdf"), p_pt | p_ph, width = 13, height = 5.5, device = cairo_pdf)

# ============================================================
# 6. Subtype x patient / cycling 诊断表
# ============================================================
sub_patient <- malig@meta.data %>%
  count(subtype_k4, Pt_number) %>%
  group_by(subtype_k4) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()
write.csv(sub_patient, file.path(tab_dir, "03c_subtype_k4_patient_composition.csv"), row.names = FALSE)

top_pat <- sub_patient %>%
  group_by(subtype_k4) %>%
  slice_max(order_by = pct, n = 1, with_ties = FALSE) %>%
  ungroup()
cat("\nTop patient per subtype:\n")
print(top_pat)

sub_phase <- malig@meta.data %>%
  count(subtype_k4, Phase) %>%
  group_by(subtype_k4) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()
write.csv(sub_phase, file.path(tab_dir, "03c_subtype_k4_phase_composition.csv"), row.names = FALSE)

# ============================================================
# 7. theta 3/5/8 sensitivity 对比表
# ============================================================
sensitivity_summary <- data.frame()
theta_files <- list(
  "3" = "GBM.malignant.reclustered.theta3.qs2",
  "5" = "GBM.malignant.reclustered.theta5.qs2",
  "8" = "GBM.malignant.reclustered.theta8.qs2"
)

for (tname in names(theta_files)) {
  fpath <- file.path(out_dir, theta_files[[tname]])
  if (!file.exists(fpath)) {
    cat("Missing:", fpath, "- skipping\n")
    next
  }
  tmp <- qs2::qs_read(fpath)
  n_cl <- length(unique(tmp[[paste0("harmony_res.", DEFAULT_RES)]][, 1]))

  sk4 <- if ("subtype_k4" %in% colnames(tmp@meta.data)) tmp$subtype_k4 else NULL
  if (!is.null(sk4)) {
    n_per_sub <- table(sk4)
    biggest_sub <- names(which.max(n_per_sub))
    cells_in_biggest <- which(sk4 == biggest_sub)
    pt_in_biggest <- table(tmp$Pt_number[cells_in_biggest])
    pt_in_biggest_pct <- prop.table(pt_in_biggest) * 100
    pt8_pct <- ifelse("Pt8" %in% names(pt_in_biggest_pct), round(pt_in_biggest_pct["Pt8"], 1), 0)
    pt18_pct <- ifelse("Pt18" %in% names(pt_in_biggest_pct), round(pt_in_biggest_pct["Pt18"], 1), 0)
  } else {
    biggest_sub <- NA
    pt8_pct <- NA
    pt18_pct <- NA
  }

  sensitivity_summary <- rbind(
    sensitivity_summary,
    data.frame(
      theta = tname,
      n_clusters_at_default_res = n_cl,
      biggest_subtype_k4 = biggest_sub,
      Pt8_in_biggest_pct = pt8_pct,
      Pt18_in_biggest_pct = pt18_pct
    )
  )
  rm(tmp)
  gc()
}

write.csv(sensitivity_summary, file.path(tab_dir, "03c_theta_sensitivity_summary.csv"), row.names = FALSE)
cat("\nTheta sensitivity summary:\n")
print(sensitivity_summary)

# ============================================================
# 8. 保存最终候选对象
# ============================================================
qs2::qs_save(malig, file.path(out_dir, "GBM.malignant.subtyped.theta5.k4.qs2"))

cat("\n========== Done ==========\n")
cat("Final object:", file.path(out_dir, "GBM.malignant.subtyped.theta5.k4.qs2"), "\n")
cat("\nFigures to review:\n")
cat("  03c_final_ultra1_Lancet.pdf\n")
cat("  03c_diagnostic_ultra1.pdf\n")
