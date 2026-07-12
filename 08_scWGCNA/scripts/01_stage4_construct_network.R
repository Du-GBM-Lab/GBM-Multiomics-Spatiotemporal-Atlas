# =============================================================================
# R2 · scWGCNA (hdWGCNA) —— STAGE 4: construct network, then stop
# -----------------------------------------------------------------------------
# Input: STAGE 3 object after TestSoftPowers.
# Output: module assignment table, dendrogram PDF, stage4 object.
#
# Boundary: this script only constructs the co-expression network. It does not
# run ModuleEigengenes, subtype association, hub interpretation, or TF/regulon
# functions.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(hdWGCNA)
  library(WGCNA)
  library(qs2)
})

set.seed(12345)
options(stringsAsFactors = FALSE)

## ===== CONFIG ================================================================
OUT_DIR <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/R2_scWGCNA"
WGCNA_NAME <- "malignant_R2"
STAGE3_OBJ <- file.path(OUT_DIR, paste0("obj_", WGCNA_NAME, "_stage3_softpower.qs2"))

SOFT_POWER <- 9
NET_TYPE <- "signed"
MIN_MODSIZE <- 30
MERGE_CUT <- 0.25
## ===========================================================================

fig_dir <- file.path(OUT_DIR, "figures")
tab_dir <- file.path(OUT_DIR, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n########## STAGE 4: ConstructNetwork ##########\n")
cat("Input:", STAGE3_OBJ, "\n")
cat("SOFT_POWER:", SOFT_POWER, "\n")
cat("networkType/TOMType:", NET_TYPE, "\n")
cat("minModuleSize:", MIN_MODSIZE, "\n")
cat("mergeCutHeight:", MERGE_CUT, "\n")

stopifnot(file.exists(STAGE3_OBJ))
obj <- qs2::qs_read(STAGE3_OBJ)

tom_dir <- file.path(OUT_DIR, "TOM")
dir.create(tom_dir, recursive = TRUE, showWarnings = FALSE)

obj <- ConstructNetwork(
  seurat_obj = obj,
  soft_power = SOFT_POWER,
  networkType = NET_TYPE,
  TOMType = NET_TYPE,
  minModuleSize = MIN_MODSIZE,
  mergeCutHeight = MERGE_CUT,
  tom_outdir = tom_dir,
  tom_name = WGCNA_NAME,
  overwrite_tom = TRUE,
  wgcna_name = WGCNA_NAME
)

modules <- GetModules(obj, wgcna_name = WGCNA_NAME)
write.csv(
  modules,
  file.path(tab_dir, "stage4_modules_gene_assignment.csv"),
  row.names = FALSE
)

module_counts <- as.data.frame(table(modules$module), stringsAsFactors = FALSE)
colnames(module_counts) <- c("module", "n_genes")
module_counts <- module_counts[order(module_counts$module), ]
write.csv(
  module_counts,
  file.path(tab_dir, "stage4_module_gene_counts.csv"),
  row.names = FALSE
)

cat("\nModule counts:\n")
print(module_counts)
cat("\ncolnames(modules):\n")
print(colnames(modules))

pdf(file.path(fig_dir, "S_stage4_dendrogram_softpower9.pdf"), width = 10, height = 5.5)
PlotDendrogram(obj, main = "R2 malignant hdWGCNA dendrogram, soft power 9")
dev.off()

qs2::qs_save(
  obj,
  file.path(OUT_DIR, paste0("obj_", WGCNA_NAME, "_stage4_network_softpower9.qs2"))
)

cat("\n########## STOP 3 reached ##########\n")
cat("Report these files before STAGE 5-7:\n")
cat("  ", file.path(tab_dir, "stage4_module_gene_counts.csv"), "\n", sep = "")
cat("  ", file.path(tab_dir, "stage4_modules_gene_assignment.csv"), "\n", sep = "")
cat("  ", file.path(fig_dir, "S_stage4_dendrogram_softpower9.pdf"), "\n", sep = "")
