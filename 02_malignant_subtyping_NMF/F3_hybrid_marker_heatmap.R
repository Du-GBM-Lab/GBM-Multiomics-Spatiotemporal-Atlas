suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(readr)
  library(Matrix)
  library(grid)
})

# F3 hybrid marker heatmap.
# Reuses existing Wilcoxon one-vs-rest DE table; does not rerun DE/NMF.

step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir) {
  setwd(file.path(getwd(), step_dir))
}

obj_path <- file.path("outputs", "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2")
de_path <- file.path("tables", "10e_one_vs_rest_FindAllMarkers_wilcoxon_downsampled.csv")

fig_dir <- "figures"
source_dir <- file.path(fig_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(fig_dir, "F3_marker_heatmap.pdf")
out_source <- file.path(source_dir, "F3_marker_heatmap_source.csv")
out_z <- file.path(source_dir, "F3_marker_heatmap_zscore_matrix.csv")
out_sanity <- file.path(source_dir, "F3_marker_heatmap_sanity_checks.csv")
out_removed <- file.path(source_dir, "F3_marker_heatmap_removed_genes.csv")

subtype_levels <- paste0("Subtype", 1:4)
subtype_colors <- c(
  "Subtype1" = "#0072B5",
  "Subtype2" = "#E18727",
  "Subtype3" = "#20854E",
  "Subtype4" = "#BC3C29"
)

forced_markers <- tibble::tribble(
  ~gene, ~assigned_subtype, ~forced_reason, ~forced_order,
  "GAP43", "Subtype1", "original_manuscript_representative", 1L,
  "TOP2A", "Subtype1", "original_manuscript_representative", 2L,
  "MAG", "Subtype2", "original_manuscript_representative", 1L,
  "COL1A1", "Subtype3", "original_manuscript_representative", 1L
)

blacklist_regex <- "^(MT-|RPL|RPS)"
removed_gene_rules <- tibble::tribble(
  ~gene, ~removal_reason,
  "AC003090.1", "uncharacterized_AC_transcript_not_used_for_marker_display",
  "AC069277.2", "uncharacterized_AC_transcript_not_used_for_marker_display",
  "AC093590.1", "uncharacterized_AC_transcript_not_used_for_marker_display",
  "AC061992.1", "uncharacterized_AC_transcript_not_used_for_marker_display",
  "RP3-406A7.7", "uncharacterized_RP_transcript_not_used_for_marker_display",
  "RP11-862L9.3", "uncharacterized_RP_transcript_not_used_for_marker_display",
  "CTD-2049O4.1", "uncharacterized_CTD_transcript_not_used_for_marker_display",
  "LINC01445", "lncRNA_not_used_for_marker_display",
  "NTM-AS1", "antisense_lncRNA_not_used_for_marker_display",
  "IGKV3D-15", "immunoglobulin_variable_gene_removed_from_malignant_marker_display",
  "IGKV3-15", "immunoglobulin_variable_gene_removed_from_malignant_marker_display"
)

blacklist_exact <- c(
  "MALAT1", "HSPA1A", "HSPA1B", "HSP90AA1", "HSP90AB1",
  removed_gene_rules$gene
)

de <- readr::read_csv(de_path, show_col_types = FALSE)
stopifnot(all(c("subtype_k4", "gene", "avg_log2FC", "BH_q", "pct.1", "pct.2") %in% colnames(de)))

de_clean <- de %>%
  mutate(
    subtype_k4 = factor(subtype_k4, levels = subtype_levels),
    blacklisted = grepl(blacklist_regex, gene) | gene %in% blacklist_exact
  )

removed_tbl <- de_clean %>%
  filter(gene %in% removed_gene_rules$gene, BH_q < 0.05, avg_log2FC > 0.5) %>%
  left_join(removed_gene_rules, by = "gene") %>%
  arrange(subtype_k4, desc(avg_log2FC), BH_q) %>%
  select(gene, subtype_k4, avg_log2FC, BH_q, pct.1, pct.2, removal_reason)
readr::write_csv(removed_tbl, out_removed)

unbiased <- de_clean %>%
  filter(!blacklisted, BH_q < 0.05, avg_log2FC > 0.5) %>%
  group_by(subtype_k4) %>%
  arrange(desc(avg_log2FC), BH_q, .by_group = TRUE) %>%
  slice_head(n = 12) %>%
  ungroup() %>%
  transmute(
    gene,
    assigned_subtype = as.character(subtype_k4),
    avg_log2FC,
    BH_q,
    pct.1,
    pct.2,
    marker_class = "unbiased",
    forced_reason = NA_character_,
    forced_order = NA_integer_
  )

forced <- forced_markers %>%
  left_join(
    de_clean %>%
      select(gene, subtype_k4, avg_log2FC, BH_q, pct.1, pct.2) %>%
      group_by(gene) %>%
      arrange(desc(avg_log2FC), BH_q, .by_group = TRUE) %>%
      slice_head(n = 1) %>%
      ungroup(),
    by = "gene"
  ) %>%
  transmute(
    gene,
    assigned_subtype,
    avg_log2FC,
    BH_q,
    pct.1,
    pct.2,
    marker_class = "manuscript_representative",
    forced_reason,
    forced_order
  )

marker_tbl <- bind_rows(unbiased, forced) %>%
  mutate(
    is_manuscript_representative = gene %in% forced_markers$gene,
    marker_class = ifelse(is_manuscript_representative, "manuscript_representative", marker_class),
    assigned_subtype = factor(assigned_subtype, levels = subtype_levels),
    marker_class_priority = ifelse(marker_class == "unbiased", 1L, 2L),
    forced_order = ifelse(is.na(forced_order), 999L, forced_order)
  ) %>%
  arrange(assigned_subtype, marker_class_priority, forced_order, desc(avg_log2FC), BH_q) %>%
  group_by(gene) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(assigned_subtype, marker_class_priority, forced_order, desc(avg_log2FC), BH_q) %>%
  mutate(row_order = row_number())

obj <- qs2::qs_read(obj_path)
DefaultAssay(obj) <- "RNA"
stopifnot("subtype_k4" %in% colnames(obj@meta.data))
obj@meta.data$subtype_k4 <- factor(as.character(obj@meta.data$subtype_k4), levels = subtype_levels)
stopifnot(!anyNA(obj@meta.data$subtype_k4))

marker_tbl <- marker_tbl %>% filter(gene %in% rownames(obj))
marker_genes <- marker_tbl$gene
stopifnot(!anyDuplicated(marker_genes))
stopifnot(all(forced_markers$gene %in% marker_genes))

expr <- GetAssayData(obj, assay = "RNA", layer = "data")[marker_genes, , drop = FALSE]
subtype_vec <- as.character(obj@meta.data$subtype_k4)

pseudobulk <- sapply(subtype_levels, function(s) {
  cells <- which(subtype_vec == s)
  Matrix::rowMeans(expr[, cells, drop = FALSE])
})
rownames(pseudobulk) <- marker_genes
colnames(pseudobulk) <- subtype_levels

z <- t(scale(t(pseudobulk)))
z[is.na(z)] <- 0
z[z > 2.5] <- 2.5
z[z < -2.5] <- -2.5

row_split <- factor(marker_tbl$assigned_subtype, levels = subtype_levels)
col_anno <- HeatmapAnnotation(
  Subtype = factor(colnames(z), levels = subtype_levels),
  col = list(Subtype = subtype_colors),
  show_annotation_name = FALSE,
  height = unit(3, "mm"),
  annotation_legend_param = list(
    Subtype = list(
      title = "Subtype",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      grid_height = unit(3, "mm"),
      grid_width = unit(3, "mm")
    )
  )
)

row_anno <- rowAnnotation(
  Subtype = row_split,
  col = list(Subtype = subtype_colors),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  width = unit(2.5, "mm")
)

heat_col <- colorRamp2(c(-2.5, 0, 2.5), c("#4393C3", "white", "#D6604D"))

ht <- Heatmap(
  z,
  name = "Z-score",
  col = heat_col,
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  cluster_columns = FALSE,
  row_split = row_split,
  column_split = factor(colnames(z), levels = subtype_levels),
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 6.5),
  row_names_side = "right",
  row_title = NULL,
  column_title = NULL,
  row_gap = unit(0.9, "mm"),
  column_gap = unit(0.8, "mm"),
  top_annotation = col_anno,
  left_annotation = row_anno,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 0.2),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7),
    legend_height = unit(2.7, "cm"),
    grid_width = unit(3, "mm")
  )
)

pdf(out_pdf, width = 4.8, height = 7.0, useDingbats = FALSE)
draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  padding = unit(c(3, 3, 3, 5), "mm")
)
dev.off()

source_tbl <- marker_tbl %>%
  mutate(
    pseudobulk_Subtype1 = pseudobulk[gene, "Subtype1"],
    pseudobulk_Subtype2 = pseudobulk[gene, "Subtype2"],
    pseudobulk_Subtype3 = pseudobulk[gene, "Subtype3"],
    pseudobulk_Subtype4 = pseudobulk[gene, "Subtype4"],
    z_Subtype1 = z[gene, "Subtype1"],
    z_Subtype2 = z[gene, "Subtype2"],
    z_Subtype3 = z[gene, "Subtype3"],
    z_Subtype4 = z[gene, "Subtype4"]
  ) %>%
  select(-marker_class_priority)

readr::write_csv(source_tbl, out_source)
write.csv(z, out_z)

sanity <- tibble::tibble(
  metric = c(
    "n_markers",
    "n_unbiased_rows",
    "n_manuscript_representative_rows",
    "forced_markers_present",
    "n_blacklisted_removed_from_candidate_pool",
    "n_removed_display_genes_passing_marker_cutoff"
  ),
  value = c(
    nrow(marker_tbl),
    sum(marker_tbl$marker_class == "unbiased"),
    sum(marker_tbl$marker_class == "manuscript_representative"),
    paste(forced_markers$gene[forced_markers$gene %in% marker_tbl$gene], collapse = ";"),
    sum(de_clean$blacklisted),
    nrow(removed_tbl)
  )
)
readr::write_csv(sanity, out_sanity)

cat("F3 markers by subtype:\n")
print(table(marker_tbl$assigned_subtype))
cat("Forced markers present:", paste(forced_markers$gene[forced_markers$gene %in% marker_tbl$gene], collapse = ", "), "\n")
cat("Removed display genes:", paste(unique(removed_tbl$gene), collapse = ", "), "\n")
cat("Output PDF:", out_pdf, "\n")
cat("Source:", out_source, "\n")
