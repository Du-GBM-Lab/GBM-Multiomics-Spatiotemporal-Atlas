# 05_恶性细胞分亚群与Neftel对照/03_cluster_top_genes_correlation_subtype.R
# v2: 对齐旧代码方法 + 多维证据共同决策
# 输出 6 类证据, k 由用户看图后定, 03 不写 subtype 回 metadata

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(cluster)
  library(dendextend)
  library(patchwork)
})

set.seed(42)

# ---- 路径 ----
in_path <- "05_恶性细胞分亚群与Neftel对照/outputs/GBM.malignant.reclustered.qs2"
out_dir <- "05_恶性细胞分亚群与Neftel对照/outputs"
tab_dir <- "05_恶性细胞分亚群与Neftel对照/tables"
fig_dir <- "05_恶性细胞分亚群与Neftel对照/figures"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 参数 ----
RESOLUTIONS <- c(0.4, 0.6, 0.8, 1.0, 1.2)
BRANCH <- "harmony"
REF_RES <- "0.8"
TOP_N_MARKER <- 50
GAP_K_MAX <- 10
GAP_B <- 50
COR_METHOD <- "spearman"
HCLUST_METHOD <- "complete"
MIN_LOGFC <- 0.25
MIN_PCT <- 0.25
K_CANDIDATES <- c(3, 4, 5, 6)
PATIENT_DOMINANCE <- 70
USE_HVG_ONLY <- TRUE

# ---- 加载 ----
malig <- qs2::qs_read(in_path)
DefaultAssay(malig) <- "RNA"
cat("Loaded:", ncol(malig), "cells\n\n")

marker_features <- if (USE_HVG_ONLY) {
  VariableFeatures(malig)
} else {
  rownames(malig)
}
if (length(marker_features) == 0) {
  stop("No marker features available. Check VariableFeatures or RNA assay.", call. = FALSE)
}
cat("Marker search features:", length(marker_features), "\n\n")

# ============================================================
# 1. 多 resolution top markers + cluster signature
# ============================================================
get_cluster_signatures <- function(obj, cluster_col, top_n = TOP_N_MARKER) {
  Idents(obj) <- cluster_col
  cl_sizes <- table(Idents(obj))
  small_cl <- names(cl_sizes)[cl_sizes < 10]
  if (length(small_cl) > 0) {
    cat("  Skipping small clusters (<10 cells):", paste(small_cl, collapse = ", "), "\n")
    obj <- subset(obj, idents = setdiff(names(cl_sizes), small_cl))
  }

  markers <- FindAllMarkers(
    obj,
    features = marker_features,
    only.pos = TRUE,
    min.pct = MIN_PCT,
    logfc.threshold = MIN_LOGFC,
    verbose = FALSE
  )

  markers <- markers %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
    ungroup()

  union_genes <- unique(markers$gene)
  if (length(union_genes) == 0) {
    stop("No significant marker genes found for ", cluster_col, call. = FALSE)
  }

  avg_expr <- AverageExpression(
    obj,
    assays = "RNA",
    features = union_genes,
    slot = "data",
    verbose = FALSE
  )$RNA
  colnames(avg_expr) <- sub("^g(?=\\d)", "", colnames(avg_expr), perl = TRUE)

  list(markers = markers, signature_mat = as.matrix(avg_expr))
}

multi_res_results <- list()
for (res in RESOLUTIONS) {
  cluster_col <- paste0(BRANCH, "_res.", res)
  if (!cluster_col %in% colnames(malig@meta.data)) {
    stop("Missing cluster column: ", cluster_col, call. = FALSE)
  }
  cat("=== Resolution", res, "(", length(unique(malig[[cluster_col]][, 1])), "clusters ) ===\n")
  multi_res_results[[as.character(res)]] <- get_cluster_signatures(malig, cluster_col)
}

all_markers_combined <- bind_rows(
  lapply(names(multi_res_results), function(r) {
    multi_res_results[[r]]$markers %>% mutate(resolution = r)
  })
)
write.csv(
  all_markers_combined,
  file.path(tab_dir, "03_top_markers_all_resolutions.csv"),
  row.names = FALSE
)

# ============================================================
# 2. cluster correlation heatmap
# ============================================================
sig_ref <- multi_res_results[[REF_RES]]$signature_mat
cor_mat <- cor(sig_ref, method = COR_METHOD)
write.csv(cor_mat, file.path(tab_dir, paste0("03_cluster_correlation_res", REF_RES, ".csv")))

pheatmap_obj <- pheatmap(
  cor_mat,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = HCLUST_METHOD,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  fontsize = 9,
  border_color = "grey60",
  main = paste0("Cluster correlation | Harmony res=", REF_RES, " | spearman | ", HCLUST_METHOD),
  filename = file.path(fig_dir, paste0("03_cluster_correlation_heatmap_res", REF_RES, ".pdf")),
  width = 9,
  height = 8
)

for (res in setdiff(names(multi_res_results), REF_RES)) {
  sig <- multi_res_results[[res]]$signature_mat
  cm <- cor(sig, method = COR_METHOD)
  write.csv(cm, file.path(tab_dir, paste0("03_cluster_correlation_res", res, ".csv")))
  pheatmap(
    cm,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = HCLUST_METHOD,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    fontsize = 8,
    border_color = "grey60",
    main = paste0("Cluster correlation | res=", res),
    filename = file.path(fig_dir, paste0("03_cluster_correlation_heatmap_res", res, ".pdf")),
    width = 9,
    height = 8,
    silent = TRUE
  )
}

# ============================================================
# 3. gap statistic: hclust + kmeans 两版
# ============================================================
n_clusters_ref <- ncol(sig_ref)
k_max <- min(GAP_K_MAX, n_clusters_ref - 1)

set.seed(42)
gap_hclust <- clusGap(
  t(sig_ref),
  FUN = function(x, k) {
    d <- as.dist(1 - cor(t(x), method = COR_METHOD))
    hc <- hclust(d, method = HCLUST_METHOD)
    list(cluster = cutree(hc, k = k))
  },
  K.max = k_max,
  B = GAP_B
)

set.seed(42)
gap_kmeans <- clusGap(
  cor_mat,
  FUN = kmeans,
  K.max = k_max,
  B = GAP_B
)

gap_combined <- bind_rows(
  as.data.frame(gap_hclust$Tab) %>%
    mutate(k = seq_len(nrow(.)), method = "hclust_complete"),
  as.data.frame(gap_kmeans$Tab) %>%
    mutate(k = seq_len(nrow(.)), method = "kmeans_on_cor")
)
write.csv(
  gap_combined,
  file.path(tab_dir, "03_gap_statistic_hclust_and_kmeans.csv"),
  row.names = FALSE
)

p_gap <- ggplot(gap_combined, aes(x = k, y = gap, color = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.2) +
  geom_vline(xintercept = K_CANDIDATES, linetype = "dotted", color = "grey40", alpha = 0.6) +
  scale_x_continuous(breaks = 1:k_max) +
  scale_color_manual(values = c("hclust_complete" = "#B40426", "kmeans_on_cor" = "#3B4CC0")) +
  labs(
    x = "Number of subtypes (k)",
    y = "Gap statistic",
    title = paste0("Gap statistic | Harmony res=", REF_RES, " | n cluster=", n_clusters_ref),
    subtitle = "Dotted lines: candidate k = 3, 4, 5, 6 | Read for elbow / plateau"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "03_gap_curve.pdf"), p_gap, width = 8, height = 5.5, device = cairo_pdf)

# ============================================================
# 4. silhouette across resolutions
# ============================================================
cross_check <- data.frame()
for (res in names(multi_res_results)) {
  sig <- multi_res_results[[res]]$signature_mat
  d <- as.dist(1 - cor(sig, method = COR_METHOD))
  hc <- hclust(d, method = HCLUST_METHOD)
  n_cl <- ncol(sig)
  for (k_val in 2:min(GAP_K_MAX, n_cl - 1)) {
    grp <- cutree(hc, k = k_val)
    if (length(unique(grp)) < 2) next
    sil <- silhouette(grp, d)
    cross_check <- rbind(
      cross_check,
      data.frame(resolution = res, k = k_val, avg_sil = mean(sil[, 3]))
    )
  }
}
write.csv(cross_check, file.path(tab_dir, "03_silhouette_across_resolutions.csv"), row.names = FALSE)

p_sil <- ggplot(cross_check, aes(x = k, y = avg_sil, color = resolution)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_vline(xintercept = K_CANDIDATES, linetype = "dotted", color = "grey40", alpha = 0.6) +
  scale_x_continuous(breaks = 1:GAP_K_MAX) +
  labs(
    x = "Number of subtypes (k)",
    y = "Average silhouette",
    title = "Silhouette across resolutions (Harmony, hclust complete)",
    subtitle = "Cross-check: which k is stable across resolutions"
  ) +
  theme_bw(base_size = 11)
ggsave(
  file.path(fig_dir, "03_silhouette_across_resolutions.pdf"),
  p_sil,
  width = 8,
  height = 5,
  device = cairo_pdf
)

# ============================================================
# 5. 多 k 候选下的 cluster -> subtype 映射, 叠 patient + cycling 注释
# ============================================================
cluster_col_ref <- paste0(BRANCH, "_res.", REF_RES)

patient_comp <- malig@meta.data %>%
  count(.data[[cluster_col_ref]], Pt_number) %>%
  group_by(.data[[cluster_col_ref]]) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  slice_max(order_by = pct, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(cluster = all_of(cluster_col_ref), top_patient = Pt_number, top_patient_n = n) %>%
  select(cluster, top_patient, top_patient_n, top_patient_pct = pct) %>%
  mutate(patient_dominated = top_patient_pct >= PATIENT_DOMINANCE)

cycling_comp <- malig@meta.data %>%
  count(.data[[cluster_col_ref]], Phase) %>%
  group_by(.data[[cluster_col_ref]]) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  filter(Phase != "G1") %>%
  summarise(cycling_pct = sum(pct), .groups = "drop") %>%
  rename(cluster = all_of(cluster_col_ref))

cluster_size <- malig@meta.data %>%
  count(.data[[cluster_col_ref]], name = "n_cells") %>%
  rename(cluster = all_of(cluster_col_ref))

hc_ref <- hclust(as.dist(1 - cor_mat), method = HCLUST_METHOD)
k_assignments <- lapply(K_CANDIDATES, function(k_val) {
  ass <- cutree(hc_ref, k = k_val)
  data.frame(cluster = names(ass), subtype = paste0("Subtype", ass), k = k_val)
})
k_assignments_df <- bind_rows(k_assignments)

cluster_audit <- cluster_size %>%
  mutate(cluster = as.character(cluster)) %>%
  left_join(patient_comp %>% mutate(cluster = as.character(cluster)), by = "cluster") %>%
  left_join(cycling_comp %>% mutate(cluster = as.character(cluster)), by = "cluster") %>%
  mutate(cycling_pct = ifelse(is.na(cycling_pct), 0, cycling_pct))

for (k_val in K_CANDIDATES) {
  sub_map <- k_assignments_df %>%
    filter(k == k_val) %>%
    select(-k)
  cluster_audit[[paste0("subtype_k", k_val)]] <- sub_map$subtype[
    match(cluster_audit$cluster, sub_map$cluster)
  ]
}

write.csv(
  cluster_audit,
  file.path(tab_dir, "03_cluster_audit_with_k_candidates.csv"),
  row.names = FALSE
)

# ============================================================
# 6. 每个候选 k: subtype 级别的 patient + cycling 主导风险
# ============================================================
subtype_risk <- list()
for (k_val in K_CANDIDATES) {
  sub_map <- k_assignments_df %>% filter(k == k_val)
  tmp_subtype <- sub_map$subtype[
    match(as.character(malig[[cluster_col_ref]][, 1]), sub_map$cluster)
  ]
  tmp_md <- malig@meta.data
  tmp_md$tmp_subtype <- tmp_subtype

  sub_pat <- tmp_md %>%
    filter(!is.na(tmp_subtype)) %>%
    count(tmp_subtype, Pt_number) %>%
    group_by(tmp_subtype) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    slice_max(order_by = pct, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(subtype = tmp_subtype, top_patient = Pt_number, top_patient_pct = pct, top_patient_n = n)

  sub_size <- tmp_md %>%
    filter(!is.na(tmp_subtype)) %>%
    count(tmp_subtype, name = "n_cells") %>%
    rename(subtype = tmp_subtype)

  sub_cycling <- tmp_md %>%
    filter(!is.na(tmp_subtype)) %>%
    count(tmp_subtype, Phase) %>%
    group_by(tmp_subtype) %>%
    mutate(pct = 100 * n / sum(n)) %>%
    filter(Phase != "G1") %>%
    summarise(cycling_pct = sum(pct), .groups = "drop") %>%
    rename(subtype = tmp_subtype)

  subtype_risk[[as.character(k_val)]] <- sub_size %>%
    left_join(sub_pat, by = "subtype") %>%
    left_join(sub_cycling, by = "subtype") %>%
    mutate(cycling_pct = ifelse(is.na(cycling_pct), 0, cycling_pct), k = k_val)
}

subtype_risk_df <- bind_rows(subtype_risk)
write.csv(
  subtype_risk_df,
  file.path(tab_dir, "03_subtype_level_risk_per_k.csv"),
  row.names = FALSE
)

cat("\n===== Subtype-level risk per k =====\n")
print(subtype_risk_df)

# ============================================================
# 7. Annotated dendrogram
# ============================================================
plot_dendro_with_annot <- function(k_show) {
  dend <- as.dendrogram(hc_ref)
  dend <- color_branches(dend, k = k_show, col = scales::hue_pal()(k_show))

  pdf(file.path(fig_dir, paste0("03_dendrogram_k", k_show, "_annotated.pdf")), width = 10, height = 6)
  par(mar = c(8, 4, 3, 1))
  plot(
    dend,
    main = paste0("Cluster dendrogram | k=", k_show, " | patient-dominated clusters listed in audit table"),
    ylab = "1 - correlation"
  )
  dev.off()
}

for (k_val in K_CANDIDATES) {
  plot_dendro_with_annot(k_val)
}

# ============================================================
# 8. UMAP overlay: 每个 k 候选画一张
# ============================================================
malig_tmp <- malig
for (k_val in K_CANDIDATES) {
  sub_map <- k_assignments_df %>% filter(k == k_val)
  malig_tmp[[paste0("subtype_k", k_val)]] <- sub_map$subtype[
    match(as.character(malig_tmp[[cluster_col_ref]][, 1]), sub_map$cluster)
  ]
}

umap_plots <- lapply(K_CANDIDATES, function(k_val) {
  DimPlot(
    malig_tmp,
    reduction = "umap_harmony",
    group.by = paste0("subtype_k", k_val),
    label = TRUE,
    raster = TRUE
  ) +
    ggtitle(paste0("k = ", k_val)) +
    theme(legend.position = "right")
})
ggsave(
  file.path(fig_dir, "03_umap_subtype_candidates.pdf"),
  patchwork::wrap_plots(umap_plots, ncol = 2),
  width = 14,
  height = 11,
  device = cairo_pdf
)

# ============================================================
# 9. 不写回 metadata, 只保存证据
# ============================================================
saveRDS(
  list(
    branch = BRANCH,
    ref_resolution = REF_RES,
    n_clusters_at_ref = n_clusters_ref,
    hc_ref = hc_ref,
    cor_mat = cor_mat,
    k_assignments = k_assignments_df,
    cluster_audit = cluster_audit,
    subtype_risk = subtype_risk_df,
    multi_res_marker_summary = sapply(multi_res_results, function(x) ncol(x$signature_mat)),
    marker_feature_scope = ifelse(USE_HVG_ONLY, "VariableFeatures", "all genes"),
    marker_feature_n = length(marker_features)
  ),
  file.path(out_dir, "03_subtype_evidence_bundle.rds")
)

cat("\nDone.\n")
cat("  Evidence bundle:", file.path(out_dir, "03_subtype_evidence_bundle.rds"), "\n")
cat("  Did NOT write subtype to metadata. Wait for 04+05 + user decision.\n")
cat("\nReview these before deciding k:\n")
cat("  fig: 03_cluster_correlation_heatmap_res0.8.pdf\n")
cat("  fig: 03_gap_curve.pdf\n")
cat("  fig: 03_silhouette_across_resolutions.pdf\n")
cat("  fig: 03_dendrogram_k{3,4,5,6}_annotated.pdf\n")
cat("  fig: 03_umap_subtype_candidates.pdf\n")
cat("  tab: 03_cluster_audit_with_k_candidates.csv\n")
cat("  tab: 03_subtype_level_risk_per_k.csv\n")
