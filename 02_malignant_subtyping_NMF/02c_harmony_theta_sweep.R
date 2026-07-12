# 05_恶性细胞分亚群与Neftel对照/02c_harmony_theta_sweep.R
# 用户决策: 对齐参数下视觉不接受, 扫 theta 看 cluster/subtype 稳健性
# v2(theta=3) 已存为 main analysis backup, 这里跑 theta=5 + theta=8

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(qs2)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(cluster)
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

# 备份现 v2 对象(theta=3)
v2_main <- file.path(out_dir, "GBM.malignant.reclustered.qs2")
v2_keep <- file.path(out_dir, "GBM.malignant.reclustered.theta3.qs2")
if (file.exists(v2_main) && !file.exists(v2_keep)) {
  file.copy(v2_main, v2_keep)
  cat("Backed up v2 (theta=3) to:", basename(v2_keep), "\n")
}

# ---- 共用参数 ----
N_HVG <- 2000
HARMONY_DIMS <- 50
DOWNSTREAM_DIMS <- 30
HARMONY_VAR <- "Pt_number"
HARMONY_LAMBDA <- 1
HARMONY_SIGMA <- 0.1
HARMONY_MAXITER <- 20
CLUSTER_RES <- c(0.4, 0.6, 0.8, 1.0, 1.2)
DEFAULT_RES <- 0.8

UMAP_GRID <- list(
  list(n.neighbors = 30, min.dist = 0.3, spread = 1.0, label = "baseline"),
  list(n.neighbors = 15, min.dist = 0.05, spread = 0.5, label = "tight"),
  list(n.neighbors = 10, min.dist = 0.01, spread = 0.3, label = "extreme"),
  list(n.neighbors = 50, min.dist = 0.1, spread = 0.7, label = "compact_balanced")
)

K_CANDIDATES <- c(3, 4, 5)

TOP_N_MARKER <- 50
COR_METHOD <- "spearman"
HCLUST_METHOD <- "complete"
GAP_K_MAX <- 10
GAP_B <- 50
MIN_LOGFC <- 0.25
MIN_PCT <- 0.25

PALETTE_NPG <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")
THETA_VALUES <- c(5, 8)

# ============================================================
# helper
# ============================================================
theme_publication <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.4, color = "black"),
      axis.ticks = element_line(linewidth = 0.4, color = "black"),
      axis.text = element_text(size = base_size - 2, color = "black"),
      axis.title = element_text(size = base_size, color = "black"),
      plot.title = element_text(size = base_size + 1, hjust = 0, face = "plain", margin = margin(b = 4)),
      legend.key = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      plot.margin = margin(6, 6, 6, 6)
    )
}

plot_umap_clean <- function(emb_df, group_col, palette, title = "", show_label = TRUE) {
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
    geom_point(size = 0.4, alpha = 0.85, stroke = 0, shape = 16) +
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

get_cluster_signatures <- function(obj, cluster_col, marker_features, top_n = TOP_N_MARKER) {
  Idents(obj) <- cluster_col
  cl_sizes <- table(Idents(obj))
  small_cl <- names(cl_sizes)[cl_sizes < 10]
  if (length(small_cl) > 0) {
    obj <- subset(obj, idents = setdiff(names(cl_sizes), small_cl))
  }

  markers <- FindAllMarkers(
    obj,
    features = marker_features,
    only.pos = TRUE,
    min.pct = MIN_PCT,
    logfc.threshold = MIN_LOGFC,
    verbose = FALSE
  ) %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
    ungroup()

  union_genes <- unique(markers$gene)
  if (length(union_genes) == 0) {
    stop("No marker genes found for ", cluster_col, call. = FALSE)
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

# ============================================================
# 主循环
# ============================================================
summary_rows <- list()

for (theta_val in THETA_VALUES) {
  tag <- paste0("theta", theta_val)
  cat("\n========================================\n")
  cat("== Running theta =", theta_val, "\n")
  cat("========================================\n")

  # ----- 加载 + 预处理 -----
  malig <- qs2::qs_read(in_path)
  if (inherits(malig[["RNA"]], "Assay5")) {
    malig <- JoinLayers(malig, assay = "RNA")
  }
  DefaultAssay(malig) <- "RNA"

  malig <- NormalizeData(malig, verbose = FALSE)
  malig <- FindVariableFeatures(malig, selection.method = "vst", nfeatures = N_HVG, verbose = FALSE)
  marker_features <- VariableFeatures(malig)
  malig <- ScaleData(malig, verbose = FALSE)
  malig <- RunPCA(malig, npcs = HARMONY_DIMS, verbose = FALSE)

  # ----- Harmony -----
  cat("  Harmony theta =", theta_val, "\n")
  malig <- RunHarmony(
    malig,
    group.by.vars = HARMONY_VAR,
    reduction.use = "pca",
    reduction.save = "harmony",
    dims.use = 1:HARMONY_DIMS,
    theta = theta_val,
    lambda = HARMONY_LAMBDA,
    sigma = HARMONY_SIGMA,
    max_iter = HARMONY_MAXITER,
    plot_convergence = FALSE,
    verbose = FALSE
  )

  # ----- 多套 UMAP -----
  for (cfg in UMAP_GRID) {
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
  }

  # ----- clustering -----
  malig <- FindNeighbors(
    malig,
    reduction = "harmony",
    dims = 1:DOWNSTREAM_DIMS,
    graph.name = c("nn_h", "snn_h"),
    verbose = FALSE
  )
  for (res in CLUSTER_RES) {
    malig <- FindClusters(
      malig,
      graph.name = "snn_h",
      resolution = res,
      cluster.name = paste0("harmony_res.", res),
      verbose = FALSE
    )
  }
  malig$cluster_default <- malig[[paste0("harmony_res.", DEFAULT_RES)]][, 1]

  # ----- cell cycle -----
  cc_genes <- Seurat::cc.genes.updated.2019
  malig <- CellCycleScoring(
    malig,
    s.features = cc_genes$s.genes,
    g2m.features = cc_genes$g2m.genes,
    set.ident = FALSE
  )

  # ----- cluster signature / gap / silhouette -----
  cluster_col_ref <- paste0("harmony_res.", DEFAULT_RES)
  sig_pack <- get_cluster_signatures(malig, cluster_col_ref, marker_features = marker_features)
  sig_mat <- sig_pack$signature_mat
  cor_mat <- cor(sig_mat, method = COR_METHOD)
  n_clusters_ref <- ncol(sig_mat)
  k_max <- min(GAP_K_MAX, n_clusters_ref - 1)

  set.seed(42)
  gap_hclust <- clusGap(
    t(sig_mat),
    FUN = function(x, k) {
      d <- as.dist(1 - cor(t(x), method = COR_METHOD))
      hc <- hclust(d, method = HCLUST_METHOD)
      list(cluster = cutree(hc, k = k))
    },
    K.max = k_max,
    B = GAP_B
  )
  gap_tab <- as.data.frame(gap_hclust$Tab) %>% mutate(k = seq_len(nrow(.)))

  d_ref <- as.dist(1 - cor_mat)
  hc_ref <- hclust(d_ref, method = HCLUST_METHOD)
  sil_tab <- data.frame()
  for (k_val in 2:k_max) {
    grp <- cutree(hc_ref, k = k_val)
    sil <- silhouette(grp, d_ref)
    sil_tab <- rbind(sil_tab, data.frame(k = k_val, avg_sil = mean(sil[, 3])))
  }

  k_cands <- intersect(K_CANDIDATES, 2:k_max)
  k_assignments <- bind_rows(lapply(k_cands, function(k_val) {
    ass <- cutree(hc_ref, k = k_val)
    data.frame(cluster = names(ass), subtype = paste0("Subtype", ass), k = k_val)
  }))

  for (k_val in k_cands) {
    sub_map <- k_assignments %>% filter(k == k_val)
    malig[[paste0("subtype_k", k_val)]] <- factor(
      sub_map$subtype[match(as.character(malig[[cluster_col_ref]][, 1]), sub_map$cluster)],
      levels = paste0("Subtype", seq_len(k_val))
    )
  }

  # ----- subtype risk -----
  patient_comp_per_subtype <- bind_rows(lapply(k_cands, function(k_val) {
    col <- paste0("subtype_k", k_val)
    malig@meta.data %>%
      filter(!is.na(.data[[col]])) %>%
      count(.data[[col]], Pt_number) %>%
      group_by(.data[[col]]) %>%
      mutate(pct = 100 * n / sum(n)) %>%
      slice_max(order_by = pct, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(k = k_val) %>%
      rename(subtype = all_of(col), top_patient_pct = pct, top_patient_n = n)
  }))

  cycling_per_subtype <- bind_rows(lapply(k_cands, function(k_val) {
    col <- paste0("subtype_k", k_val)
    malig@meta.data %>%
      filter(!is.na(.data[[col]])) %>%
      count(.data[[col]], Phase) %>%
      group_by(.data[[col]]) %>%
      mutate(pct = 100 * n / sum(n)) %>%
      filter(Phase != "G1") %>%
      summarise(cycling_pct = sum(pct), .groups = "drop") %>%
      mutate(k = k_val) %>%
      rename(subtype = all_of(col))
  }))

  subtype_risk <- patient_comp_per_subtype %>%
    select(k, subtype, Pt_number, top_patient_pct, top_patient_n) %>%
    left_join(cycling_per_subtype, by = c("k", "subtype")) %>%
    mutate(cycling_pct = ifelse(is.na(cycling_pct), 0, cycling_pct))

  # ----- 输出表 -----
  sweep_tab_dir <- file.path(tab_dir, paste0("02c_sweep_", tag))
  dir.create(sweep_tab_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(gap_tab, file.path(sweep_tab_dir, "gap_curve.csv"), row.names = FALSE)
  write.csv(sil_tab, file.path(sweep_tab_dir, "silhouette.csv"), row.names = FALSE)
  write.csv(patient_comp_per_subtype, file.path(sweep_tab_dir, "patient_dominance_per_subtype_per_k.csv"), row.names = FALSE)
  write.csv(cycling_per_subtype, file.path(sweep_tab_dir, "cycling_per_subtype_per_k.csv"), row.names = FALSE)
  write.csv(subtype_risk, file.path(sweep_tab_dir, "subtype_risk_per_k.csv"), row.names = FALSE)
  write.csv(sig_pack$markers, file.path(sweep_tab_dir, "top50_markers_default_resolution.csv"), row.names = FALSE)
  write.csv(cor_mat, file.path(sweep_tab_dir, "cluster_correlation_matrix.csv"))
  write.csv(k_assignments, file.path(sweep_tab_dir, "cluster_to_subtype_k_candidates.csv"), row.names = FALSE)
  write.csv(
    data.frame(
      theta = theta_val,
      resolution = CLUSTER_RES,
      n_clusters = sapply(CLUSTER_RES, function(r) length(unique(malig[[paste0("harmony_res.", r)]][, 1])))
    ),
    file.path(sweep_tab_dir, "cluster_count_per_resolution.csv"),
    row.names = FALSE
  )

  # ----- 图: gap + silhouette -----
  p_gap <- ggplot(gap_tab, aes(x = k, y = gap)) +
    geom_line(linewidth = 0.8, color = "#B40426") +
    geom_point(size = 2, color = "#B40426") +
    geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.2, color = "#B40426") +
    geom_vline(xintercept = K_CANDIDATES, linetype = "dotted", color = "grey40", alpha = 0.6) +
    scale_x_continuous(breaks = 1:k_max) +
    labs(
      x = "Number of subtypes (k)",
      y = "Gap statistic",
      title = paste0("Gap statistic | theta=", theta_val, " | n cluster=", n_clusters_ref)
    ) +
    theme_publication()

  p_sil <- ggplot(sil_tab, aes(x = k, y = avg_sil)) +
    geom_line(linewidth = 0.8, color = "#3C5488") +
    geom_point(size = 2, color = "#3C5488") +
    geom_vline(xintercept = K_CANDIDATES, linetype = "dotted", color = "grey40", alpha = 0.6) +
    scale_x_continuous(breaks = 1:k_max) +
    labs(x = "Number of subtypes (k)", y = "Average silhouette", title = paste0("Silhouette | theta=", theta_val)) +
    theme_publication()

  ggsave(file.path(fig_dir, paste0("02c_gap_silhouette_", tag, ".pdf")), p_gap | p_sil, width = 12, height = 5, device = cairo_pdf)

  # ----- 图: cor heatmap -----
  pheatmap(
    cor_mat,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = HCLUST_METHOD,
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    fontsize = 9,
    border_color = "grey60",
    main = paste0("Cluster correlation | theta=", theta_val, " | res=", DEFAULT_RES, " | n cluster=", n_clusters_ref),
    filename = file.path(fig_dir, paste0("02c_correlation_heatmap_", tag, ".pdf")),
    width = 9,
    height = 8
  )

  # ----- 图: UMAP grid -----
  for (k_val in k_cands) {
    plots <- lapply(UMAP_GRID, function(cfg) {
      rd_name <- paste0("umap_", cfg$label)
      emb <- as.data.frame(Embeddings(malig, rd_name))
      colnames(emb) <- c("UMAP_1", "UMAP_2")
      emb[[paste0("subtype_k", k_val)]] <- malig[[paste0("subtype_k", k_val)]][, 1]
      emb[[paste0("subtype_k", k_val)]] <- factor(emb[[paste0("subtype_k", k_val)]], levels = paste0("Subtype", seq_len(k_val)))
      plot_umap_clean(
        emb,
        paste0("subtype_k", k_val),
        PALETTE_NPG,
        title = paste0(cfg$label, " (nn=", cfg$n.neighbors, " md=", cfg$min.dist, ")"),
        show_label = TRUE
      )
    })
    ggsave(
      file.path(fig_dir, sprintf("02c_umap_grid_%s_k%d.pdf", tag, k_val)),
      wrap_plots(plots, ncol = 2) +
        plot_annotation(
          title = sprintf("UMAP grid | theta=%d | k=%d", theta_val, k_val),
          theme = theme(plot.title = element_text(face = "bold"))
        ),
      width = 13,
      height = 11,
      device = cairo_pdf
    )
  }

  # ----- Pt8/Pt18 + Phase on extreme UMAP -----
  emb_ext <- as.data.frame(Embeddings(malig, "umap_extreme"))
  colnames(emb_ext) <- c("UMAP_1", "UMAP_2")
  emb_ext$pt_hl <- factor(
    ifelse(malig$Pt_number == "Pt8", "Pt8", ifelse(malig$Pt_number == "Pt18", "Pt18", "Other")),
    levels = c("Other", "Pt18", "Pt8")
  )

  p_hl <- plot_umap_clean(
    emb_ext,
    "pt_hl",
    c("Other" = "grey85", "Pt18" = "#3C5488", "Pt8" = "#E64B35"),
    title = paste0("Pt8/Pt18 | theta=", theta_val, " | extreme UMAP"),
    show_label = FALSE
  )
  ggsave(file.path(fig_dir, paste0("02c_Pt8_Pt18_extreme_", tag, ".pdf")), p_hl, width = 6.5, height = 5.5, device = cairo_pdf)

  emb_ph <- emb_ext
  emb_ph$Phase <- factor(malig$Phase, levels = c("G1", "S", "G2M"))
  p_ph <- plot_umap_clean(
    emb_ph,
    "Phase",
    c("G1" = "grey75", "S" = "#3C5488", "G2M" = "#E64B35"),
    title = paste0("Cell cycle | theta=", theta_val),
    show_label = FALSE
  )
  ggsave(file.path(fig_dir, paste0("02c_phase_extreme_", tag, ".pdf")), p_ph, width = 6.5, height = 5.5, device = cairo_pdf)

  # ----- 保存对象 -----
  qs2::qs_save(malig, file.path(out_dir, paste0("GBM.malignant.reclustered.", tag, ".qs2")))

  summary_rows[[tag]] <- data.frame(
    theta = theta_val,
    n_clusters_res0.8 = n_clusters_ref,
    gap_best_k = gap_tab$k[which.max(gap_tab$gap)],
    silhouette_best_k = sil_tab$k[which.max(sil_tab$avg_sil)],
    max_top_patient_pct = max(subtype_risk$top_patient_pct, na.rm = TRUE),
    max_cycling_pct = max(subtype_risk$cycling_pct, na.rm = TRUE)
  )

  cat("  theta =", theta_val, "done.\n")
  cat("    n clusters at res", DEFAULT_RES, ":", n_clusters_ref, "\n")
  cat("    subtype risk:\n")
  print(subtype_risk)

  rm(malig)
  gc()
}

summary_df <- bind_rows(summary_rows)
write.csv(summary_df, file.path(tab_dir, "02c_theta_sweep_summary.csv"), row.names = FALSE)

cat("\n\n========================================\n")
cat("Sweep complete.\n")
cat("========================================\n")
print(summary_df)
cat("\n对象:\n")
cat("  GBM.malignant.reclustered.theta3.qs2  (v2 backup)\n")
cat("  GBM.malignant.reclustered.theta5.qs2\n")
cat("  GBM.malignant.reclustered.theta8.qs2\n")
