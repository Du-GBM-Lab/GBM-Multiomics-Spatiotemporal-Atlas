# 05_恶性细胞分亚群与Neftel对照/04c_subtype3_vs_4_diagnostic.R
# Diagnose whether Subtype3/4 separation reflects MES1/MES2 biology, other programs, or patient skew.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(Seurat)
  library(qs2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

set.seed(42)

proj <- "05_恶性细胞分亚群与Neftel对照"
in_obj <- file.path(proj, "outputs", "GBM.malignant.subtyped.neftel_scored.submodule_labeled.qs2")
fig_dir <- file.path(proj, "figures")
tab_dir <- file.path(proj, "tables")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")
pick_col <- function(md, candidates) {
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit) == 0) {
    stop("None of these metadata columns exist: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  hit[1]
}

msg("Loading object:", in_obj)
obj <- qs2::qs_read(in_obj)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

md <- obj@meta.data
required <- c("subtype_k4", "MES1_UCell", "MES2_UCell")
missing <- setdiff(required, colnames(md))
if (length(missing) > 0) {
  stop("Missing required metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
}
patient_col <- pick_col(md, c("Pt_number", "patient", "sample", "orig.ident"))
msg("Patient/sample column:", patient_col)

if (!all(c("Subtype3", "Subtype4") %in% md$subtype_k4)) {
  stop("Subtype3/Subtype4 not both present in subtype_k4.", call. = FALSE)
}

msg("Check 1: FindMarkers Subtype3 vs Subtype4.")
Idents(obj) <- "subtype_k4"
de34 <- FindMarkers(
  obj,
  ident.1 = "Subtype3",
  ident.2 = "Subtype4",
  min.pct = 0.25,
  logfc.threshold = 0.5
)
de34$gene <- rownames(de34)
fc_col <- pick_col(de34, c("avg_log2FC", "avg_logFC"))
de34 <- de34 |> arrange(desc(.data[[fc_col]]))
write.csv(de34, file.path(tab_dir, "04c_DE_S3_vs_S4.csv"), row.names = FALSE)

top_s3 <- de34 |>
  arrange(desc(.data[[fc_col]])) |>
  select(gene, all_of(fc_col), pct.1, pct.2, p_val_adj) |>
  head(15)
top_s4 <- de34 |>
  arrange(.data[[fc_col]]) |>
  select(gene, all_of(fc_col), pct.1, pct.2, p_val_adj) |>
  head(15)

write.csv(top_s3, file.path(tab_dir, "04c_DE_top15_Subtype3_vs_S4.csv"), row.names = FALSE)
write.csv(top_s4, file.path(tab_dir, "04c_DE_top15_Subtype4_vs_S3.csv"), row.names = FALSE)

cat("\n== Top 15 UP in Subtype3 (vs S4) ==\n")
print(top_s3)
cat("\n== Top 15 UP in Subtype4 (vs S3) ==\n")
print(top_s4)

msg("Check 2: MES1/MES2 Spearman correlations.")
rho_global <- cor(md$MES1_UCell, md$MES2_UCell, method = "spearman", use = "complete.obs")
rho_in_s34 <- with(
  md[md$subtype_k4 %in% c("Subtype3", "Subtype4"), ],
  cor(MES1_UCell, MES2_UCell, method = "spearman", use = "complete.obs")
)
rho_df <- data.frame(
  scope = c("global", "Subtype3+Subtype4"),
  spearman = c(rho_global, rho_in_s34)
)
write.csv(rho_df, file.path(tab_dir, "04c_MES1_MES2_spearman.csv"), row.names = FALSE)
cat(sprintf("\n== Global Spearman MES1 ~ MES2 (UCell): %.3f ==\n", rho_global))
cat(sprintf("Spearman within S3+S4 only: %.3f\n", rho_in_s34))

msg("Check 3: MES1/MES2 violin and medians.")
mes_long <- md |>
  select(subtype_k4, MES1 = MES1_UCell, MES2 = MES2_UCell) |>
  pivot_longer(-subtype_k4, names_to = "module", values_to = "score")

p_mes <- ggplot(mes_long, aes(subtype_k4, score, fill = module)) +
  geom_violin(scale = "width", alpha = 0.70, position = position_dodge(0.8), linewidth = 0.2) +
  geom_boxplot(width = 0.12, position = position_dodge(0.8), outlier.shape = NA, alpha = 0.55, linewidth = 0.25) +
  scale_fill_manual(values = c(MES1 = "#B24745", MES2 = "#2F6FA3")) +
  theme_classic(base_size = 12) +
  labs(title = "MES1 vs MES2 UCell scores per subtype", y = "UCell score", x = NULL, fill = NULL)

ggsave(
  file.path(fig_dir, "04c_MES1_MES2_violin.pdf"),
  p_mes,
  width = 7,
  height = 4.5,
  device = cairo_pdf
)

mes_medians <- md |>
  group_by(subtype_k4) |>
  summarise(
    n_cells = n(),
    MES1_med = median(MES1_UCell),
    MES2_med = median(MES2_UCell),
    MES2_minus_MES1 = MES2_med - MES1_med,
    .groups = "drop"
  )
write.csv(mes_medians, file.path(tab_dir, "04c_MES1_MES2_medians_by_subtype.csv"), row.names = FALSE)
cat("\n== MES1 / MES2 medians per subtype ==\n")
print(mes_medians)

msg("Check 4: patient composition in Subtype3/4.")
pt_in_s34 <- md[md$subtype_k4 %in% c("Subtype3", "Subtype4"), ]
pt_tab <- prop.table(table(pt_in_s34[[patient_col]], pt_in_s34$subtype_k4), margin = 2) * 100
pt_counts <- as.data.frame.matrix(table(pt_in_s34[[patient_col]], pt_in_s34$subtype_k4))
pt_pct <- as.data.frame.matrix(round(pt_tab, 2))
pt_counts <- cbind(patient = rownames(pt_counts), pt_counts)
pt_pct <- cbind(patient = rownames(pt_pct), pt_pct)
rownames(pt_counts) <- NULL
rownames(pt_pct) <- NULL

write.csv(pt_counts, file.path(tab_dir, "04c_patient_composition_S3_S4_counts.csv"), row.names = FALSE)
write.csv(pt_pct, file.path(tab_dir, "04c_patient_composition_S3_S4_pct.csv"), row.names = FALSE)

cat("\n== Patient composition within S3 / S4 (% of subtype) ==\n")
print(round(pt_tab[order(-pt_tab[, "Subtype3"]), ], 2))

top1_s3 <- max(pt_tab[, "Subtype3"])
top1_s4 <- max(pt_tab[, "Subtype4"])
top3_s3 <- sum(sort(pt_tab[, "Subtype3"], decreasing = TRUE)[1:3])
top3_s4 <- sum(sort(pt_tab[, "Subtype4"], decreasing = TRUE)[1:3])

patient_concentration <- data.frame(
  subtype = c("Subtype3", "Subtype4"),
  top1_patient_pct = c(top1_s3, top1_s4),
  top3_patients_pct = c(top3_s3, top3_s4)
)
write.csv(patient_concentration, file.path(tab_dir, "04c_patient_concentration_S3_S4.csv"), row.names = FALSE)

cat(sprintf("\nS3: top1 patient %.1f%%, top3 patients %.1f%%\n", top1_s3, top3_s3))
cat(sprintf("S4: top1 patient %.1f%%, top3 patients %.1f%%\n", top1_s4, top3_s4))

cat("\nWritten files:\n")
cat("- tables/04c_DE_S3_vs_S4.csv\n")
cat("- tables/04c_MES1_MES2_spearman.csv\n")
cat("- tables/04c_MES1_MES2_medians_by_subtype.csv\n")
cat("- tables/04c_patient_composition_S3_S4_pct.csv\n")
cat("- tables/04c_patient_concentration_S3_S4.csv\n")
cat("- figures/04c_MES1_MES2_violin.pdf\n")
msg("Done.")
