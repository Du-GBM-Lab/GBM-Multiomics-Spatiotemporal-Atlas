# 05_恶性细胞分亚群与Neftel对照/02_recluster_malignant_cells.R
# v2: 对齐旧代码 Harmony / UMAP / HVG 参数
# 1. 备份 v1 对象
# 2. 用旧代码参数重跑 Harmony + UMAP + cluster
# 3. no-integration 分支保留作 sensitivity

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(qs2)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# ---- 路径 ----
in_path <- "05_恶性细胞分亚群与Neftel对照/outputs/GBM.malignant.high_confidence.qs2"
out_dir <- "05_恶性细胞分亚群与Neftel对照/outputs"
tab_dir <- "05_恶性细胞分亚群与Neftel对照/tables"
fig_dir <- "05_恶性细胞分亚群与Neftel对照/figures"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 备份 v1 对象 ----
v1_targets <- c(
  "GBM.malignant.reclustered.qs2",
  "GBM.malignant.subtyped.qs2"
)

for (f in v1_targets) {
  src <- file.path(out_dir, f)
  dst <- file.path(out_dir, sub("\\.qs2$", ".v1_default_harmony.qs2", f))
  if (file.exists(src) && !file.exists(dst)) {
    file.rename(src, dst)
    cat("Backed up:", f, "->", basename(dst), "\n")
  }
}

# ---- 参数: 对齐旧代码 ----
N_HVG <- 2000
HARMONY_DIMS <- 50
DOWNSTREAM_DIMS <- 30
HARMONY_VAR <- "Pt_number"
HARMONY_THETA <- 3
HARMONY_LAMBDA <- 1
HARMONY_SIGMA <- 0.1
HARMONY_MAXITER <- 20
UMAP_NN <- 30
UMAP_MINDIST <- 1
CLUSTER_RES <- c(0.4, 0.6, 0.8, 1.0, 1.2)
DEFAULT_RES <- 0.8

# ---- 加载 ----
malig <- qs2::qs_read(in_path)
if (inherits(malig[["RNA"]], "Assay5")) {
  malig <- JoinLayers(malig, assay = "RNA")
}
DefaultAssay(malig) <- "RNA"

# 清掉旧 reduction / cluster, 避免继承 upstream 对象结果
for (rd in Reductions(malig)) {
  malig[[rd]] <- NULL
}

old_cluster_cols <- grep(
  "^(seurat_clusters|RNA_snn_res|integrated_snn_res|SCT_snn_res|noint_res|harmony_res|cluster_noint_default|cluster_harmony_default|subtype_raw)",
  colnames(malig@meta.data),
  value = TRUE
)
if (length(old_cluster_cols) > 0) {
  malig@meta.data[, old_cluster_cols] <- NULL
}

cat("Input:", ncol(malig), "malignant cells x", nrow(malig), "genes\n")
cat("Patients:", nlevels(droplevels(factor(malig$Pt_number))), "\n\n")

# ============================================================
# 共用预处理
# ============================================================
malig <- NormalizeData(malig, verbose = FALSE)
malig <- FindVariableFeatures(
  malig,
  selection.method = "vst",
  nfeatures = N_HVG,
  verbose = FALSE
)
malig <- ScaleData(malig, verbose = FALSE)
malig <- RunPCA(malig, npcs = HARMONY_DIMS, verbose = FALSE)

# ============================================================
# Branch A: no-integration (sensitivity)
# ============================================================
cat("=== Branch A: no-integration ===\n")
obj_noint <- malig
obj_noint <- RunUMAP(
  obj_noint,
  dims = 1:DOWNSTREAM_DIMS,
  reduction = "pca",
  n.neighbors = UMAP_NN,
  min.dist = UMAP_MINDIST,
  reduction.name = "umap_noint",
  verbose = FALSE
)
obj_noint <- FindNeighbors(
  obj_noint,
  dims = 1:DOWNSTREAM_DIMS,
  reduction = "pca",
  graph.name = c("nn_noint", "snn_noint"),
  verbose = FALSE
)
for (res in CLUSTER_RES) {
  obj_noint <- FindClusters(
    obj_noint,
    graph.name = "snn_noint",
    resolution = res,
    cluster.name = paste0("noint_res.", res),
    verbose = FALSE
  )
}
obj_noint$cluster_noint_default <- obj_noint[[paste0("noint_res.", DEFAULT_RES)]][, 1]

# ============================================================
# Branch B: Harmony
# ============================================================
cat("=== Branch B: Harmony (theta=3, lambda=1, sigma=0.1, max.iter=20) ===\n")
obj_harm <- malig
obj_harm <- RunHarmony(
  obj_harm,
  group.by.vars = HARMONY_VAR,
  reduction.use = "pca",
  reduction.save = "harmony",
  dims.use = 1:HARMONY_DIMS,
  theta = HARMONY_THETA,
  lambda = HARMONY_LAMBDA,
  sigma = HARMONY_SIGMA,
  max_iter = HARMONY_MAXITER,
  plot_convergence = FALSE,
  verbose = FALSE
)
obj_harm <- RunUMAP(
  obj_harm,
  dims = 1:DOWNSTREAM_DIMS,
  reduction = "harmony",
  n.neighbors = UMAP_NN,
  min.dist = UMAP_MINDIST,
  reduction.name = "umap_harmony",
  verbose = FALSE
)
obj_harm <- FindNeighbors(
  obj_harm,
  dims = 1:DOWNSTREAM_DIMS,
  reduction = "harmony",
  graph.name = c("nn_harmony", "snn_harmony"),
  verbose = FALSE
)
for (res in CLUSTER_RES) {
  obj_harm <- FindClusters(
    obj_harm,
    graph.name = "snn_harmony",
    resolution = res,
    cluster.name = paste0("harmony_res.", res),
    verbose = FALSE
  )
}
obj_harm$cluster_harmony_default <- obj_harm[[paste0("harmony_res.", DEFAULT_RES)]][, 1]

# ============================================================
# 合并两版结果
# ============================================================
malig[["pca"]] <- obj_noint[["pca"]]
malig[["umap_noint"]] <- obj_noint[["umap_noint"]]
malig[["harmony"]] <- obj_harm[["harmony"]]
malig[["umap_harmony"]] <- obj_harm[["umap_harmony"]]

noint_cols <- c(paste0("noint_res.", CLUSTER_RES), "cluster_noint_default")
harmony_cols <- c(paste0("harmony_res.", CLUSTER_RES), "cluster_harmony_default")
malig@meta.data[, noint_cols] <- obj_noint@meta.data[, noint_cols]
malig@meta.data[, harmony_cols] <- obj_harm@meta.data[, harmony_cols]

# ============================================================
# Cell cycle scoring
# ============================================================
cc_genes <- Seurat::cc.genes.updated.2019
malig <- CellCycleScoring(
  malig,
  s.features = cc_genes$s.genes,
  g2m.features = cc_genes$g2m.genes,
  set.ident = FALSE
)

# ============================================================
# 诊断: patient x cluster 占比 + 单 patient 主导 cluster
# ============================================================
diag_dir <- file.path(tab_dir, "02_diagnostics_v2")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

for (branch in c("noint", "harmony")) {
  col <- paste0("cluster_", branch, "_default")
  tab <- table(malig$Pt_number, malig[[col]][, 1])
  pct_by_cluster <- prop.table(tab, margin = 2) * 100
  write.csv(
    round(pct_by_cluster, 2),
    file.path(diag_dir, paste0("02_patient_x_cluster_pct_", branch, "_res", DEFAULT_RES, ".csv"))
  )

  dominated <- which(pct_by_cluster > 70, arr.ind = TRUE)
  cat("[", branch, "] Single-patient-dominated clusters (>70%): ", nrow(dominated), "\n", sep = "")
  if (nrow(dominated) > 0) {
    for (i in seq_len(nrow(dominated))) {
      cat(
        "  Cluster", colnames(pct_by_cluster)[dominated[i, "col"]],
        "<-", rownames(pct_by_cluster)[dominated[i, "row"]],
        "(", round(pct_by_cluster[dominated[i, "row"], dominated[i, "col"]], 1), "%)\n"
      )
    }
  }
}

cluster_count_summary <- data.frame(
  branch = c("no-integration", "harmony"),
  n_clusters_at_default = c(
    length(unique(malig$cluster_noint_default)),
    length(unique(malig$cluster_harmony_default))
  )
)
write.csv(
  cluster_count_summary,
  file.path(diag_dir, "02_cluster_count_summary.csv"),
  row.names = FALSE
)
cat("\nCluster counts at res", DEFAULT_RES, ":\n")
print(cluster_count_summary)

# ============================================================
# 可视化
# ============================================================
plot_umap <- function(obj, reduction, group_col, title) {
  DimPlot(
    obj,
    reduction = reduction,
    group.by = group_col,
    label = TRUE,
    label.size = 3,
    raster = TRUE
  ) +
    ggtitle(title) +
    theme(legend.position = "right")
}

p1 <- plot_umap(
  malig,
  "umap_noint",
  "cluster_noint_default",
  paste0("No-integration | res=", DEFAULT_RES)
)
p2 <- plot_umap(
  malig,
  "umap_harmony",
  "cluster_harmony_default",
  paste0("Harmony (theta=3) | res=", DEFAULT_RES)
)
p3 <- plot_umap(malig, "umap_noint", "Pt_number", "No-integration | by Patient")
p4 <- plot_umap(malig, "umap_harmony", "Pt_number", "Harmony (theta=3) | by Patient")

ggsave(
  file.path(fig_dir, "02_umap_cluster_comparison_v2.pdf"),
  (p1 | p2) / (p3 | p4),
  width = 16,
  height = 12,
  device = cairo_pdf
)

# Pt8 / Pt18 高亮: 检查前一版 Subtype2 的 patient residual 风险
malig$is_Pt8_Pt18 <- ifelse(
  malig$Pt_number %in% c("Pt8", "Pt18"),
  as.character(malig$Pt_number),
  "Other"
)
p_pt_a <- DimPlot(
  malig,
  reduction = "umap_noint",
  group.by = "is_Pt8_Pt18",
  cols = c("Pt8" = "#B40426", "Pt18" = "#3B4CC0", "Other" = "grey85"),
  order = c("Pt8", "Pt18"),
  raster = TRUE
) +
  ggtitle("No-integration | Pt8 / Pt18 highlight")
p_pt_b <- DimPlot(
  malig,
  reduction = "umap_harmony",
  group.by = "is_Pt8_Pt18",
  cols = c("Pt8" = "#B40426", "Pt18" = "#3B4CC0", "Other" = "grey85"),
  order = c("Pt8", "Pt18"),
  raster = TRUE
) +
  ggtitle("Harmony (theta=3) | Pt8 / Pt18 highlight")
ggsave(
  file.path(fig_dir, "02_Pt8_Pt18_highlight_v2.pdf"),
  p_pt_a | p_pt_b,
  width = 14,
  height = 6,
  device = cairo_pdf
)

# Pt1 高亮
malig$is_Pt1 <- ifelse(malig$Pt_number == "Pt1", "Pt1 (n=33)", "Other")
p_pt1_a <- DimPlot(
  malig,
  reduction = "umap_noint",
  group.by = "is_Pt1",
  cols = c("Pt1 (n=33)" = "#B40426", "Other" = "grey85"),
  order = "Pt1 (n=33)",
  raster = TRUE
) +
  ggtitle("No-integration | Pt1")
p_pt1_b <- DimPlot(
  malig,
  reduction = "umap_harmony",
  group.by = "is_Pt1",
  cols = c("Pt1 (n=33)" = "#B40426", "Other" = "grey85"),
  order = "Pt1 (n=33)",
  raster = TRUE
) +
  ggtitle("Harmony (theta=3) | Pt1")
ggsave(
  file.path(fig_dir, "02_Pt1_highlight_v2.pdf"),
  p_pt1_a | p_pt1_b,
  width = 14,
  height = 6,
  device = cairo_pdf
)

# Cell cycle phase
p_cc_a <- DimPlot(malig, reduction = "umap_noint", group.by = "Phase", raster = TRUE) +
  ggtitle("No-integration | Phase")
p_cc_b <- DimPlot(malig, reduction = "umap_harmony", group.by = "Phase", raster = TRUE) +
  ggtitle("Harmony (theta=3) | Phase")
ggsave(
  file.path(fig_dir, "02_cellcycle_phase_v2.pdf"),
  p_cc_a | p_cc_b,
  width = 14,
  height = 6,
  device = cairo_pdf
)

# ============================================================
# 保存
# ============================================================
qs2::qs_save(malig, file.path(out_dir, "GBM.malignant.reclustered.qs2"))

cat("\nDone.\n")
cat("  Object:", file.path(out_dir, "GBM.malignant.reclustered.qs2"), "\n")
cat(
  "  Harmony params: theta=", HARMONY_THETA,
  " lambda=", HARMONY_LAMBDA,
  " sigma=", HARMONY_SIGMA,
  " max.iter=", HARMONY_MAXITER, "\n",
  sep = ""
)
cat("  UMAP params: n.neighbors=", UMAP_NN, " min.dist=", UMAP_MINDIST, "\n", sep = "")
cat("  HVG: ", N_HVG, "\n", sep = "")
cat("  PCA dims: 1:", HARMONY_DIMS, " | downstream dims: 1:", DOWNSTREAM_DIMS, "\n", sep = "")
cat("\nNext: 03 (corrected version) - manual k selection from gap curve.\n")
