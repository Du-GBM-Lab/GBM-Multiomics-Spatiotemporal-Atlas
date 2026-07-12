# 05_恶性细胞分亚群与Neftel对照/05d_umap_final_round.R
# Final display-only UMAP round: compare closer4 with three honest parameter candidates.
# This script does not change Harmony, clustering, subtype assignment, or metadata.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(Seurat)
  library(qs2)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

proj <- "05_恶性细胞分亚群与Neftel对照"
in_obj <- file.path(proj, "outputs", "GBM.malignant.subtyped.umap_candidates.qs2")
fig_dir <- file.path(proj, "figures")
tab_dir <- file.path(proj, "tables")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

LANCET_K4 <- c(
  "Subtype1" = "#00468B",
  "Subtype2" = "#ED0000",
  "Subtype3" = "#42B540",
  "Subtype4" = "#0099B4"
)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

theme_embed_compare <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 11, face = "bold", hjust = 0),
      plot.margin = margin(5, 5, 5, 5)
    )
}

embedding_df <- function(obj, reduction) {
  emb <- as.data.frame(Embeddings(obj, reduction))
  colnames(emb)[1:2] <- c("dim1", "dim2")
  emb$subtype_k4 <- factor(obj$subtype_k4, levels = names(LANCET_K4))
  emb <- emb[sample(nrow(emb)), ]
  emb
}

plot_one <- function(obj, red, title) {
  emb <- embedding_df(obj, red)
  centers <- aggregate(
    emb[, c("dim1", "dim2")],
    list(subtype_k4 = emb$subtype_k4),
    median
  )
  centers$label <- sub("Subtype", "S", as.character(centers$subtype_k4))

  p <- ggplot(emb, aes(dim1, dim2, color = subtype_k4)) +
    geom_point(size = 0.30, alpha = 0.62, stroke = 0, shape = 16) +
    scale_color_manual(values = LANCET_K4, drop = FALSE) +
    labs(title = title) +
    theme_embed_compare()

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_label_repel(
      data = centers,
      aes(dim1, dim2, label = label),
      inherit.aes = FALSE,
      size = 4.5,
      fontface = "bold",
      color = "black",
      fill = "white",
      alpha = 0.92,
      label.size = 0.20,
      box.padding = 0.30,
      point.padding = 0.20,
      segment.color = NA
    )
  } else {
    p <- p + geom_label(
      data = centers,
      aes(dim1, dim2, label = label),
      inherit.aes = FALSE,
      size = 4.5,
      fontface = "bold",
      color = "black",
      fill = "white",
      alpha = 0.92,
      label.size = 0.20
    )
  }

  p
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

if (!"umap_closer4" %in% Reductions(obj)) {
  msg("umap_closer4 not found; recomputing closer4 display UMAP.")
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
}

msg("Computing e1: dims=1:10, nn=120, md=0.30, spread=0.80.")
obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = 1:10,
  n.neighbors = 120,
  min.dist = 0.30,
  spread = 0.80,
  reduction.name = "umap_e1",
  reduction.key = "UMAPe1_",
  seed.use = 42,
  verbose = FALSE
)

msg("Computing e2: dims=1:15, nn=120, md=0.10, spread=0.80.")
obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = 1:15,
  n.neighbors = 120,
  min.dist = 0.10,
  spread = 0.80,
  reduction.name = "umap_e2",
  reduction.key = "UMAPe2_",
  seed.use = 42,
  verbose = FALSE
)

msg("Computing e3: dims=1:15, nn=50, md=0.30, spread=1.20.")
obj <- RunUMAP(
  obj,
  reduction = "harmony",
  dims = 1:15,
  n.neighbors = 50,
  min.dist = 0.30,
  spread = 1.20,
  reduction.name = "umap_e3",
  reduction.key = "UMAPe3_",
  seed.use = 42,
  verbose = FALSE
)

p_closer4 <- plot_one(obj, "umap_closer4", "closer4: dims=1:30 | nn120 md0.30 sp0.8")
p_e1 <- plot_one(obj, "umap_e1", "e1: dims=1:10 | nn120 md0.30 sp0.8")
p_e2 <- plot_one(obj, "umap_e2", "e2: dims=1:15 | nn120 md0.10 sp0.8")
p_e3 <- plot_one(obj, "umap_e3", "e3: dims=1:15 | nn50 md0.30 sp1.2")

combined <- (p_closer4 | p_e1) / (p_e2 | p_e3) +
  plot_annotation(
    title = "Final round UMAP candidates",
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(
  file.path(fig_dir, "05d_umap_final_round_compare.pdf"),
  combined,
  width = 12,
  height = 10,
  device = cairo_pdf
)

param_table <- data.frame(
  reduction = c("umap_closer4", "umap_e1", "umap_e2", "umap_e3"),
  source_reduction = "harmony",
  dims = c("1:30", "1:10", "1:15", "1:15"),
  n_neighbors = c(120, 120, 120, 50),
  min_dist = c(0.30, 0.30, 0.10, 0.30),
  spread = c(0.80, 0.80, 0.80, 1.20),
  note = c(
    "current closer4 baseline",
    "lower-dimensional input, tighter global summary",
    "d15 with tighter min.dist",
    "d15 with lower neighbors and larger spread"
  )
)
write.csv(param_table, file.path(tab_dir, "05d_umap_final_round_parameters.csv"), row.names = FALSE)

msg("Saving object with e1/e2/e3 reductions:", in_obj)
qs2::qs_save(obj, in_obj)

cat("\nDone. Review: figures/05d_umap_final_round_compare.pdf\n")
cat("After choosing one panel, move to 05b subtype naming.\n")
