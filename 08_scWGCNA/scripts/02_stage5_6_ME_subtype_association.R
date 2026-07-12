# =============================================================================
# R2 · scWGCNA (hdWGCNA) —— STAGE 5-6: ME-subtype association, then stop
# -----------------------------------------------------------------------------
# Input: STAGE 4 object after ConstructNetwork.
# Output: kME-enabled module table, hub genes, harmonized ME by subtype, and
# module-subtype association.
#
# Boundary: do not run STAGE 7 focus-gene lookup here. The MES-V module and kME
# column naming must be reviewed at STOP 4 before querying PLAUR/AP-1-FOSL genes.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(hdWGCNA)
  library(WGCNA)
  library(qs2)
  library(pheatmap)
})

set.seed(12345)
options(stringsAsFactors = FALSE)

## ===== CONFIG ================================================================
OUT_DIR <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/R2_scWGCNA"
WGCNA_NAME <- "malignant_R2"
STAGE4_OBJ <- file.path(OUT_DIR, paste0("obj_", WGCNA_NAME, "_stage4_network_softpower9.qs2"))

PATIENT_COL <- "Pt_number"
ASSAY <- "RNA"
## ===========================================================================

fig_dir <- file.path(OUT_DIR, "figures")
tab_dir <- file.path(OUT_DIR, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n########## STAGE 5: ModuleEigengenes + ModuleConnectivity ##########\n")
cat("Input:", STAGE4_OBJ, "\n")
stopifnot(file.exists(STAGE4_OBJ))
obj <- qs2::qs_read(STAGE4_OBJ)
stopifnot(PATIENT_COL %in% colnames(obj@meta.data))
stopifnot("subtype_bio" %in% colnames(obj@meta.data))

obj <- ModuleEigengenes(
  seurat_obj = obj,
  group.by.vars = PATIENT_COL,
  assay = ASSAY,
  wgcna_name = WGCNA_NAME
)

obj <- ModuleConnectivity(
  seurat_obj = obj,
  group.by = "malignant_all",
  group_name = "malignant",
  assay = ASSAY,
  slot = "data",
  layer = "data",
  harmonized = TRUE,
  wgcna_name = WGCNA_NAME
)

modules <- GetModules(obj, wgcna_name = WGCNA_NAME)
hub_df <- GetHubGenes(obj, n_hubs = 25, wgcna_name = WGCNA_NAME)

write.csv(
  modules,
  file.path(tab_dir, "stage5_modules_gene_assignment_kME.csv"),
  row.names = FALSE
)
write.csv(
  hub_df,
  file.path(tab_dir, "stage5_hub_genes_top25.csv"),
  row.names = FALSE
)

cat("\ncolnames(modules) after ModuleConnectivity:\n")
print(colnames(modules))
cat("\nModule counts after kME:\n")
print(table(modules$module))

cat("\n########## STAGE 6: hME-subtype association ##########\n")
hMEs <- GetMEs(obj, harmonized = TRUE, wgcna_name = WGCNA_NAME)
hMEs <- as.data.frame(hMEs)
cat("dim(hMEs):", paste(dim(hMEs), collapse = " x "), "\n")
cat("colnames(hMEs):\n")
print(colnames(hMEs))

sub <- obj$subtype_bio
if (length(sub) != nrow(hMEs)) {
  common_cells <- intersect(rownames(hMEs), colnames(obj))
  if (!length(common_cells)) {
    stop("Cannot align hMEs rows to object cells.")
  }
  hMEs <- hMEs[common_cells, , drop = FALSE]
  sub <- obj$subtype_bio[match(common_cells, colnames(obj))]
}

module_from_me <- function(x) sub("^ME", "", x)
me_cols <- colnames(hMEs)
module_labels <- module_from_me(me_cols)
keep <- module_labels != "grey"
me_cols <- me_cols[keep]
module_labels <- module_labels[keep]

me_by_sub <- t(vapply(seq_along(me_cols), function(i) {
  tapply(hMEs[[me_cols[i]]], sub, mean, na.rm = TRUE)
}, numeric(length(unique(sub)))))
rownames(me_by_sub) <- module_labels
colnames(me_by_sub) <- names(tapply(hMEs[[me_cols[1]]], sub, mean, na.rm = TRUE))
me_by_sub <- me_by_sub[, c("NPC-P", "OPC-M", "MES-V", "MES-I"), drop = FALSE]

write.csv(me_by_sub, file.path(tab_dir, "ME_mean_by_subtype.csv"))

kw_p <- vapply(seq_along(me_cols), function(i) {
  kruskal.test(hMEs[[me_cols[i]]], factor(sub))$p.value
}, numeric(1))
assoc <- data.frame(
  module = module_labels,
  top_subtype = colnames(me_by_sub)[max.col(me_by_sub)],
  KW_p = kw_p,
  KW_padj = p.adjust(kw_p, method = "BH"),
  stringsAsFactors = FALSE
)
assoc <- assoc[order(assoc$top_subtype, assoc$KW_padj), ]
write.csv(
  assoc,
  file.path(tab_dir, "module_subtype_association.csv"),
  row.names = FALSE
)

cat("\nmodule_subtype_association:\n")
print(assoc)

if ("turquoise" %in% rownames(me_by_sub)) {
  cat("\nturquoise hME means by subtype:\n")
  print(me_by_sub["turquoise", , drop = FALSE])
}

pdf(file.path(fig_dir, "B_ME_subtype_heatmap.pdf"), width = 6.5, height = 7.5)
pheatmap::pheatmap(
  me_by_sub,
  scale = "row",
  cluster_cols = FALSE,
  main = "Module eigengene activity across malignant subtypes",
  color = colorRampPalette(c("#0072B5", "white", "#BC3C29"))(50)
)
dev.off()

qs2::qs_save(
  obj,
  file.path(OUT_DIR, paste0("obj_", WGCNA_NAME, "_stage6_ME_subtype.qs2"))
)

cat("\n########## STOP 4 reached ##########\n")
cat("Report these files before STAGE 7:\n")
cat("  ", file.path(tab_dir, "stage5_modules_gene_assignment_kME.csv"), "\n", sep = "")
cat("  ", file.path(tab_dir, "stage5_hub_genes_top25.csv"), "\n", sep = "")
cat("  ", file.path(tab_dir, "module_subtype_association.csv"), "\n", sep = "")
cat("  ", file.path(tab_dir, "ME_mean_by_subtype.csv"), "\n", sep = "")
cat("  ", file.path(fig_dir, "B_ME_subtype_heatmap.pdf"), "\n", sep = "")
