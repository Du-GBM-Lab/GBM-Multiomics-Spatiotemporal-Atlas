# 05_恶性细胞分亚群与Neftel对照/08_patient_subtype_composition.R
# Patient x subtype composition sanity check.
# Purpose: show cohort diversity and defend against single-patient artifact.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

set.seed(42)

proj <- "05_恶性细胞分亚群与Neftel对照"
in_obj <- file.path(proj, "outputs", "GBM.malignant.subtyped.named.qs2")
fig_dir <- file.path(proj, "figures")
tab_dir <- file.path(proj, "tables")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

subtype_levels <- c("NPC-Cycling", "OPC-like", "MES-Perivascular", "MES-Inflammatory")
subtype_palette <- c(
  "NPC-Cycling" = "#00468B",
  "OPC-like" = "#ED0000",
  "MES-Perivascular" = "#42B540",
  "MES-Inflammatory" = "#0099B4"
)

pick_col <- function(md, candidates) {
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit) == 0) {
    stop("None of these patient/sample columns exist: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  hit[1]
}

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

msg("Loading object:", in_obj)
obj <- qs2::qs_read(in_obj)
md <- obj@meta.data

stopifnot("subtype_named" %in% colnames(md))
patient_col <- pick_col(md, c("patient", "Pt_number", "sample", "orig.ident"))
msg("Using patient/sample column:", patient_col)

md$subtype_named <- factor(as.character(md$subtype_named), levels = subtype_levels)
md$patient_for_08 <- as.character(md[[patient_col]])
if (any(is.na(md$subtype_named))) {
  stop("subtype_named contains NA after factor level locking.", call. = FALSE)
}

# 1. Counts matrix.
cnt_mat <- as.matrix(table(md$patient_for_08, md$subtype_named))
cnt <- as.data.frame.matrix(cnt_mat)
cnt$patient <- rownames(cnt)
cnt$total <- rowSums(cnt[, subtype_levels, drop = FALSE])
cnt <- cnt[order(-cnt$total), c("patient", subtype_levels, "total")]

write.csv(cnt, file.path(tab_dir, "08_patient_subtype_counts.csv"), row.names = FALSE)

# 2. Row-normalized: subtype proportion within each patient.
count_only <- as.matrix(cnt[, subtype_levels, drop = FALSE])
rownames(count_only) <- cnt$patient
prop_by_pt <- prop.table(count_only, margin = 1) * 100
write.csv(
  data.frame(patient = rownames(prop_by_pt), prop_by_pt, check.names = FALSE),
  file.path(tab_dir, "08_patient_subtype_prop_rownorm.csv"),
  row.names = FALSE
)

# 3. Column-normalized: patient contribution within each subtype.
prop_by_st <- prop.table(count_only, margin = 2) * 100
write.csv(
  data.frame(patient = rownames(prop_by_st), prop_by_st, check.names = FALSE),
  file.path(tab_dir, "08_subtype_patient_contribution.csv"),
  row.names = FALSE
)

# 4. Sanity report.
sanity <- data.frame(
  subtype = colnames(prop_by_st),
  top1_patient = apply(prop_by_st, 2, function(x) rownames(prop_by_st)[which.max(x)]),
  top1_patient_pct = as.numeric(apply(prop_by_st, 2, max)),
  top3_patients_pct = as.numeric(apply(prop_by_st, 2, function(x) sum(sort(x, decreasing = TRUE)[1:3]))),
  n_patients_gt1pct = as.integer(apply(prop_by_st, 2, function(x) sum(x > 1))),
  n_patients_gt5pct = as.integer(apply(prop_by_st, 2, function(x) sum(x > 5))),
  stringsAsFactors = FALSE
)

write.csv(sanity, file.path(tab_dir, "08_subtype_concentration_sanity.csv"), row.names = FALSE)

cat("\n== Subtype patient concentration (sanity) ==\n")
print(sanity)

# 5. Heatmap, patient x subtype row-normalized percent.
dominant_st <- apply(prop_by_pt, 1, function(x) colnames(prop_by_pt)[which.max(x)])
pt_order <- order(
  factor(dominant_st, levels = subtype_levels),
  -apply(prop_by_pt, 1, max),
  rownames(prop_by_pt)
)
prop_by_pt_ord <- prop_by_pt[pt_order, , drop = FALSE]

dominant_ordered <- dominant_st[pt_order]
row_ha <- rowAnnotation(
  dominant = dominant_ordered,
  col = list(dominant = subtype_palette),
  show_annotation_name = FALSE,
  width = unit(3, "mm")
)

col_fun <- colorRamp2(c(0, 50, 100), c("white", "#FDB863", "#B2182B"))

ht <- Heatmap(
  prop_by_pt_ord,
  name = "% in patient",
  col = col_fun,
  left_annotation = row_ha,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  column_names_rot = 30,
  column_names_gp = gpar(fontsize = 10),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.0f", prop_by_pt_ord[i, j]), x, y, gp = gpar(fontsize = 7, col = "black"))
  },
  width = unit(10, "cm"),
  height = unit(0.4, "cm") * nrow(prop_by_pt_ord)
)

pdf(
  file.path(fig_dir, "08_patient_subtype_heatmap.pdf"),
  width = 7,
  height = max(6, 0.3 * nrow(prop_by_pt_ord) + 2),
  useDingbats = FALSE
)
draw(ht, column_title = "Subtype composition per patient (%)")
dev.off()

png(
  file.path(fig_dir, "08_patient_subtype_heatmap.png"),
  width = 7,
  height = max(6, 0.3 * nrow(prop_by_pt_ord) + 2),
  units = "in",
  res = 300
)
draw(ht, column_title = "Subtype composition per patient (%)")
dev.off()

# 6. Stacked barplot.
prop_long <- as.data.frame(prop_by_pt_ord) |>
  tibble::rownames_to_column("patient") |>
  pivot_longer(-patient, names_to = "subtype", values_to = "pct") |>
  mutate(
    patient = factor(patient, levels = rownames(prop_by_pt_ord)),
    subtype = factor(subtype, levels = subtype_levels)
  )

p_bar <- ggplot(prop_long, aes(patient, pct, fill = subtype)) +
  geom_col(width = 0.85) +
  scale_fill_manual(values = subtype_palette, drop = FALSE) +
  labs(title = "Subtype composition per patient", x = NULL, y = "% of malignant cells", fill = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right",
    legend.key.size = unit(0.38, "cm"),
    plot.title = element_text(face = "bold", size = 12)
  )

ggsave(
  file.path(fig_dir, "08_patient_subtype_stacked_bar.pdf"),
  p_bar,
  width = max(7, 0.3 * nrow(prop_by_pt_ord)),
  height = 4.5,
  device = cairo_pdf
)

ggsave(
  file.path(fig_dir, "08_patient_subtype_stacked_bar.png"),
  p_bar,
  width = max(7, 0.3 * nrow(prop_by_pt_ord)),
  height = 4.5,
  dpi = 300
)

cat("\nDone. Outputs:\n")
cat("- tables/08_patient_subtype_counts.csv\n")
cat("- tables/08_patient_subtype_prop_rownorm.csv\n")
cat("- tables/08_subtype_patient_contribution.csv\n")
cat("- tables/08_subtype_concentration_sanity.csv\n")
cat("- figures/08_patient_subtype_heatmap.pdf/png\n")
cat("- figures/08_patient_subtype_stacked_bar.pdf/png\n")
