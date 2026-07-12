# 05_恶性细胞分亚群与Neftel对照/05c_umap_tryout_dims15_phate.R
# Display-only embedding candidates for malignant subtype overview.
# Adds closer4, Harmony dims=1:15 UMAP, and optional PHATE reductions.
# Does not change cluster, subtype assignment, Harmony, or metadata labels.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(Seurat)
  library(qs2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

set.seed(42)

proj <- "05_恶性细胞分亚群与Neftel对照"
in_obj <- file.path(proj, "outputs", "GBM.malignant.subtyped.neftel_scored.submodule_labeled.qs2")
fig_dir <- file.path(proj, "figures")
tab_dir <- file.path(proj, "tables")
out_obj <- file.path(proj, "outputs", "GBM.malignant.subtyped.umap_candidates.qs2")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

LANCET_K4 <- c(
  "Subtype1" = "#00468B",
  "Subtype2" = "#ED0000",
  "Subtype3" = "#42B540",
  "Subtype4" = "#0099B4"
)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

theme_embed_clean <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0, margin = margin(b = 4)),
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.key.size = unit(0.42, "cm"),
      plot.margin = margin(6, 6, 6, 6)
    )
}

get_embedding_df <- function(obj, reduction) {
  emb <- as.data.frame(Embeddings(obj, reduction))
  colnames(emb)[1:2] <- c("dim1", "dim2")
  emb$subtype_k4 <- factor(obj$subtype_k4, levels = names(LANCET_K4))
  emb$cell <- rownames(emb)
  emb
}

add_labels <- function(p, emb, label_size = 5) {
  centers <- emb |>
    group_by(subtype_k4) |>
    summarise(
      dim1 = median(dim1),
      dim2 = median(dim2),
      label = sub("Subtype", "S", first(as.character(subtype_k4))),
      .groups = "drop"
    )

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p + ggrepel::geom_label_repel(
      data = centers,
      aes(dim1, dim2, label = label),
      inherit.aes = FALSE,
      size = label_size,
      fontface = "bold",
      color = "black",
      fill = "white",
      alpha = 0.92,
      label.size = 0.22,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = NA
    )
  } else {
    p + geom_label(
      data = centers,
      aes(dim1, dim2, label = label),
      inherit.aes = FALSE,
      size = label_size,
      fontface = "bold",
      color = "black",
      fill = "white",
      alpha = 0.92,
      label.size = 0.22
    )
  }
}

plot_overview <- function(obj, reduction, title) {
  emb <- get_embedding_df(obj, reduction)
  emb <- emb[sample(nrow(emb)), ]

  p <- ggplot(emb, aes(dim1, dim2, color = subtype_k4)) +
    stat_ellipse(
      aes(fill = subtype_k4),
      geom = "polygon",
      type = "norm",
      level = 0.72,
      alpha = 0.08,
      linewidth = 0,
      show.legend = FALSE
    ) +
    stat_ellipse(
      type = "norm",
      level = 0.72,
      linewidth = 0.55,
      alpha = 0.95,
      show.legend = FALSE
    ) +
    geom_point(size = 0.30, alpha = 0.60, stroke = 0, shape = 16) +
    scale_color_manual(values = LANCET_K4, drop = FALSE) +
    scale_fill_manual(values = LANCET_K4, drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = title) +
    theme_embed_clean()

  add_labels(p, emb)
}

plot_highlight_panels <- function(obj, reduction, title) {
  emb <- get_embedding_df(obj, reduction)
  panels <- lapply(names(LANCET_K4), function(st) {
    bg <- emb
    fg <- emb[emb$subtype_k4 == st, ]
    ggplot() +
      geom_point(data = bg, aes(dim1, dim2), color = "grey84", size = 0.20, alpha = 0.30, stroke = 0) +
      geom_point(data = fg, aes(dim1, dim2), color = LANCET_K4[[st]], size = 0.34, alpha = 0.76, stroke = 0) +
      labs(title = sub("Subtype", "S", st)) +
      theme_embed_clean(base_size = 11) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
  })
  wrap_plots(panels, ncol = 2) + plot_annotation(title = title)
}

plot_compare_grid <- function(obj, reductions) {
  plots <- lapply(names(reductions), function(label) {
    plot_overview(obj, reductions[[label]], label) +
      theme(legend.position = "none")
  })
  wrap_plots(plots, ncol = length(plots)) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right")
}

msg("Loading object:", in_obj)
obj <- qs2::qs_read(in_obj)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

stopifnot("harmony" %in% Reductions(obj))
stopifnot("subtype_k4" %in% colnames(obj@meta.data))
obj$subtype_k4 <- factor(as.character(obj$subtype_k4), levels = names(LANCET_K4))

msg("Computing display UMAP: closer4, Harmony dims=1:30.")
obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 120,
  min.dist = 0.30,
  spread = 0.80,
  reduction.name = "umap_closer4",
  reduction.key = "UMAPcloser4_",
  seed.use = 42,
  verbose = FALSE
)

msg("Computing display UMAP: d15, Harmony dims=1:15.")
obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = 1:15,
  n.neighbors = 120,
  min.dist = 0.30,
  spread = 0.80,
  reduction.name = "umap_d15",
  reduction.key = "UMAPd15_",
  seed.use = 42,
  verbose = FALSE
)

phate_ok <- requireNamespace("phateR", quietly = TRUE)
if (!phate_ok) {
  msg("phateR is not installed; PHATE skipped.")
} else {
  msg("Computing PHATE on Harmony dims=1:30.")
  harmony_emb <- Embeddings(obj, "harmony")[, 1:30]
  phate_res <- phateR::phate(
    harmony_emb,
    knn = 30,
    t = "auto",
    seed = 42,
    verbose = FALSE
  )
  phate_mat <- as.matrix(phate_res$embedding)
  colnames(phate_mat) <- c("PHATE_1", "PHATE_2")
  obj[["phate"]] <- CreateDimReducObject(
    embeddings = phate_mat,
    key = "PHATE_",
    assay = DefaultAssay(obj)
  )
}

msg("Saving overview and highlight figures.")
p_closer4 <- plot_overview(obj, "umap_closer4", "UMAP closer4 | Harmony dims=1:30")
p_d15 <- plot_overview(obj, "umap_d15", "UMAP d15 | Harmony dims=1:15")

ggsave(file.path(fig_dir, "05c_umap_closer4_overview_outlined.pdf"), p_closer4, width = 6.8, height = 5.6, device = cairo_pdf)
ggsave(file.path(fig_dir, "05c_umap_d15_overview_outlined.pdf"), p_d15, width = 6.8, height = 5.6, device = cairo_pdf)

ggsave(
  file.path(fig_dir, "05c_umap_closer4_highlight_panels.pdf"),
  plot_highlight_panels(obj, "umap_closer4", "Subtype highlights | closer4"),
  width = 7.2,
  height = 6.6,
  device = cairo_pdf
)
ggsave(
  file.path(fig_dir, "05c_umap_d15_highlight_panels.pdf"),
  plot_highlight_panels(obj, "umap_d15", "Subtype highlights | d15"),
  width = 7.2,
  height = 6.6,
  device = cairo_pdf
)

reductions_for_grid <- c("closer4" = "umap_closer4", "d15" = "umap_d15")
if (phate_ok) {
  p_phate <- plot_overview(obj, "phate", "PHATE | Harmony dims=1:30")
  ggsave(file.path(fig_dir, "05c_phate_overview_outlined.pdf"), p_phate, width = 6.8, height = 5.6, device = cairo_pdf)
  ggsave(
    file.path(fig_dir, "05c_phate_highlight_panels.pdf"),
    plot_highlight_panels(obj, "phate", "Subtype highlights | PHATE"),
    width = 7.2,
    height = 6.6,
    device = cairo_pdf
  )
  reductions_for_grid <- c(reductions_for_grid, "PHATE" = "phate")
}

ggsave(
  file.path(fig_dir, "05c_embedding_candidate_overview_grid.pdf"),
  plot_compare_grid(obj, reductions_for_grid),
  width = ifelse(phate_ok, 16, 11),
  height = 5.8,
  device = cairo_pdf
)

param_table <- data.frame(
  reduction = c("umap_closer4", "umap_d15", if (phate_ok) "phate" else character(0)),
  source_reduction = c("harmony", "harmony", if (phate_ok) "harmony" else character(0)),
  dims = c("1:30", "1:15", if (phate_ok) "1:30" else character(0)),
  method = c("UMAP", "UMAP", if (phate_ok) "PHATE" else character(0)),
  n_neighbors_or_knn = c(120, 120, if (phate_ok) 30 else numeric(0)),
  min_dist = c(0.30, 0.30, if (phate_ok) NA_real_ else numeric(0)),
  spread = c(0.80, 0.80, if (phate_ok) NA_real_ else numeric(0)),
  phate_installed = phate_ok
)
write.csv(param_table, file.path(tab_dir, "05c_embedding_candidate_parameters.csv"), row.names = FALSE)

msg("Saving object with candidate reductions:", out_obj)
qs2::qs_save(obj, out_obj)

cat("\nDone. Compare these PDFs:\n")
cat("- figures/05c_embedding_candidate_overview_grid.pdf\n")
cat("- figures/05c_umap_closer4_overview_outlined.pdf\n")
cat("- figures/05c_umap_d15_overview_outlined.pdf\n")
cat("- figures/05c_umap_closer4_highlight_panels.pdf\n")
cat("- figures/05c_umap_d15_highlight_panels.pdf\n")
if (phate_ok) {
  cat("- figures/05c_phate_overview_outlined.pdf\n")
  cat("- figures/05c_phate_highlight_panels.pdf\n")
} else {
  cat("- PHATE skipped because phateR is not installed.\n")
}
