# 05_恶性细胞分亚群与Neftel对照/10d_S3_vs_S4_marker_DE.R
# Rule v2 Evidence 2: RNA-assay DE between Subtype3 and Subtype4 plus
# marker heatmap draft and marker x NMF-metaprogram overlap.

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(qs2)
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

set.seed(42)

params <- list(
  input_object = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.v2.qs2"
  ),
  tables_dir = file.path("05_恶性细胞分亚群与Neftel对照", "tables"),
  outputs_dir = file.path("05_恶性细胞分亚群与Neftel对照", "outputs"),
  metaprogram_summary_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10b_metaprogram_summary.csv"
  ),
  metaprogram_signature_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10c_metaprogram_signatures_used.csv"
  ),
  metaprogram_scores_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10c_per_cell_metaprogram_scores.tsv"
  ),
  subtype_levels = paste0("Subtype", 1:4),
  ident_1 = "Subtype3",
  ident_2 = "Subtype4",
  logfc_threshold_for_test = 0,
  min_pct_for_test = 0.1,
  marker_log2fc_cutoff = 0.5,
  marker_q_cutoff = 0.05,
  heatmap_cells_per_subtype = 200L,
  heatmap_markers_per_direction = 15L
)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

write_session_info <- function(path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  sink(con)
  print(sessionInfo())
  sink()
}

split_genes <- function(x) {
  genes <- trimws(unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE))
  genes[!is.na(genes) & genes != ""]
}

rank_markers <- function(df) {
  df |>
    mutate(
      rank_score = abs(avg_log2FC) * -log10(pmax(BH_q, .Machine$double.xmin))
    ) |>
    arrange(desc(rank_score), BH_q, desc(abs(avg_log2FC)))
}

msg("Loading object:", params$input_object)
obj <- qs2::qs_read(params$input_object)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

md <- obj@meta.data
stopifnot("subtype_k4" %in% colnames(md))
md$subtype_k4 <- factor(as.character(md$subtype_k4), levels = params$subtype_levels)
obj@meta.data$subtype_k4 <- md$subtype_k4
Idents(obj) <- "subtype_k4"

msg("Running FindMarkers on RNA assay:", params$ident_1, "vs", params$ident_2)
mk <- FindMarkers(
  obj,
  ident.1 = params$ident_1,
  ident.2 = params$ident_2,
  assay = "RNA",
  logfc.threshold = params$logfc_threshold_for_test,
  min.pct = params$min_pct_for_test,
  test.use = "wilcox",
  only.pos = FALSE,
  verbose = FALSE
)

mk$gene <- rownames(mk)
if (!"avg_log2FC" %in% colnames(mk) && "avg_logFC" %in% colnames(mk)) {
  mk$avg_log2FC <- mk$avg_logFC
}
stopifnot(all(c("p_val", "avg_log2FC", "pct.1", "pct.2", "gene") %in% colnames(mk)))
mk$BH_q <- p.adjust(mk$p_val, method = "BH")
mk <- mk |>
  relocate(gene) |>
  arrange(BH_q, desc(abs(avg_log2FC)))

all_de_file <- file.path(params$tables_dir, "10d_S3_vs_S4_FindMarkers_all.csv")
write_csv(mk, all_de_file)

passing <- mk |>
  filter(BH_q < params$marker_q_cutoff, abs(avg_log2FC) > params$marker_log2fc_cutoff)
passing_file <- file.path(params$tables_dir, "10d_S3_vs_S4_passing_markers.csv")
write_csv(passing, passing_file)

s3_top20 <- passing |>
  filter(avg_log2FC > 0) |>
  rank_markers() |>
  slice_head(n = 20)
s4_top20 <- passing |>
  filter(avg_log2FC < 0) |>
  rank_markers() |>
  slice_head(n = 20)

write_csv(s3_top20, file.path(params$tables_dir, "10d_S3_markers_top20.csv"))
write_csv(s4_top20, file.path(params$tables_dir, "10d_S4_markers_top20.csv"))

mp_signatures <- read_csv(params$metaprogram_signature_file, show_col_types = FALSE)
get_mp_genes <- function(mp_id) {
  genes <- mp_signatures |>
    filter(metaprogram == mp_id, used_in_object) |>
    pull(gene)
  stopifnot(length(genes) >= 20)
  genes
}
mp02_genes <- get_mp_genes("MP02")
mp04_genes <- get_mp_genes("MP04")
mp06_genes <- get_mp_genes("MP06")

s3_top20_genes <- s3_top20$gene
s4_top20_genes <- s4_top20$gene
overlap_df <- tibble(
  direction = c("S3_up", "S3_up", "S4_up"),
  vs_metaprogram = c("MP02_perivascular", "MP06_endothelial", "MP04_inflammatory"),
  n_overlap = c(
    length(intersect(s3_top20_genes, mp02_genes)),
    length(intersect(s3_top20_genes, mp06_genes)),
    length(intersect(s4_top20_genes, mp04_genes))
  ),
  overlap_genes = c(
    paste(intersect(s3_top20_genes, mp02_genes), collapse = ","),
    paste(intersect(s3_top20_genes, mp06_genes), collapse = ","),
    paste(intersect(s4_top20_genes, mp04_genes), collapse = ",")
  )
)
write_csv(overlap_df, file.path(params$tables_dir, "10d_marker_x_metaprogram_overlap.csv"))

evidence2_n_passing <- nrow(passing)
evidence2_decision <- tibble(
  evidence = "evidence2_marker_DE",
  n_tested_genes = nrow(mk),
  n_passing = evidence2_n_passing,
  n_S3_up = sum(passing$avg_log2FC > 0),
  n_S4_up = sum(passing$avg_log2FC < 0),
  statistical_threshold_passed = evidence2_n_passing >= 10,
  visually_separable = NA,
  evidence2_satisfied = NA
)
write_csv(evidence2_decision, file.path(params$tables_dir, "10d_evidence2_decision.csv"))

msg("Building draft marker heatmap")
mp_scores <- read_tsv(params$metaprogram_scores_file, show_col_types = FALSE)
stopifnot(all(c("cell_id", "MP02", "MP04") %in% colnames(mp_scores)))
mp_scores <- as.data.frame(mp_scores)
rownames(mp_scores) <- mp_scores$cell_id

s3_cells_all <- rownames(md)[md$subtype_k4 == params$ident_1]
s4_cells_all <- rownames(md)[md$subtype_k4 == params$ident_2]
stopifnot(length(s3_cells_all) >= params$heatmap_cells_per_subtype)
stopifnot(length(s4_cells_all) >= params$heatmap_cells_per_subtype)
s3_cells <- sample(s3_cells_all, params$heatmap_cells_per_subtype)
s4_cells <- sample(s4_cells_all, params$heatmap_cells_per_subtype)
cells_use <- c(s3_cells, s4_cells)

marker_genes <- unique(c(
  head(s3_top20$gene, params$heatmap_markers_per_direction),
  head(s4_top20$gene, params$heatmap_markers_per_direction)
))
marker_genes <- intersect(marker_genes, rownames(obj))
stopifnot(length(marker_genes) >= 20)

expr <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[marker_genes, cells_use, drop = FALSE])
expr_scaled <- t(scale(t(expr)))
expr_scaled[is.na(expr_scaled)] <- 0
expr_scaled[expr_scaled > 2.5] <- 2.5
expr_scaled[expr_scaled < -2.5] <- -2.5

subtype_vec <- factor(as.character(md[cells_use, "subtype_k4"]), levels = c(params$ident_1, params$ident_2))
row_split <- ifelse(marker_genes %in% s3_top20$gene, "S3_up", "S4_up")
row_split <- factor(row_split, levels = c("S3_up", "S4_up"))

col_anno <- HeatmapAnnotation(
  Subtype = subtype_vec,
  MP02_perivascular = mp_scores[cells_use, "MP02"],
  MP04_inflammatory = mp_scores[cells_use, "MP04"],
  col = list(
    Subtype = c(Subtype3 = "#1F77B4", Subtype4 = "#FF7F0E"),
    MP02_perivascular = colorRamp2(c(0, 1), c("white", "#1F77B4")),
    MP04_inflammatory = colorRamp2(c(0, 1), c("white", "#FF7F0E"))
  ),
  annotation_name_gp = gpar(fontsize = 8),
  show_annotation_name = TRUE
)

heatmap_file <- file.path(params$outputs_dir, "10d_S3_vs_S4_marker_heatmap_draft.pdf")
pdf(heatmap_file, width = 8, height = 7, useDingbats = FALSE)
draw(
  Heatmap(
    expr_scaled,
    name = "Z-score",
    col = colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "white", "#B2182B")),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    row_split = row_split,
    column_split = subtype_vec,
    top_annotation = col_anno,
    column_title = NULL,
    row_title_gp = gpar(fontsize = 9),
    row_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
  ),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

sanity <- tibble(
  n_tested_genes = nrow(mk),
  n_passing = nrow(passing),
  n_S3_up = sum(passing$avg_log2FC > 0),
  n_S4_up = sum(passing$avg_log2FC < 0),
  n_heatmap_cells = length(cells_use),
  n_heatmap_genes = length(marker_genes),
  all_passing_q_lt_0_05 = all(passing$BH_q < params$marker_q_cutoff),
  all_passing_abs_log2fc_gt_0_5 = all(abs(passing$avg_log2FC) > params$marker_log2fc_cutoff)
)
write_csv(sanity, file.path(params$tables_dir, "10d_sanity_checks.csv"))

stopifnot(nrow(mk) > 5000)
stopifnot(all(passing$BH_q < params$marker_q_cutoff))
stopifnot(all(abs(passing$avg_log2FC) > params$marker_log2fc_cutoff))
stopifnot(sum(passing$avg_log2FC > 0) > 0)
stopifnot(sum(passing$avg_log2FC < 0) > 0)
stopifnot(file.exists(heatmap_file))

write_session_info(file.path(params$tables_dir, "10d_session_info.txt"))

msg("10d completed")
print(evidence2_decision)
print(overlap_df)
