# =============================================================================
# R2 · scWGCNA (hdWGCNA) —— STAGE 0-3: soft power audit
# -----------------------------------------------------------------------------
# Design: one malignant-only co-expression network. Metacells are grouped by
# patient only; subtype labels are used later for association, not network
# construction.
#
# Boundary: module != regulon; hub gene = high kME co-expression hub, not a
# causal regulator. Do not call hdWGCNA TF-network functions in R2.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(hdWGCNA)
  library(WGCNA)
  library(Matrix)
  library(qs2)
  library(patchwork)
})

set.seed(12345)
options(stringsAsFactors = FALSE)

## ===== CONFIG ================================================================
OBJ_PATH <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/05_恶性细胞分亚群与Neftel对照/outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"
BL_PATH <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/06_恶性细胞拟时序/data/gene_blacklist.txt"
OUT_DIR <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/R2_scWGCNA"
WGCNA_NAME <- "malignant_R2"

SUBTYPE_COL <- "subtype_label_final"
PATIENT_COL <- "Pt_number"
ASSAY <- "RNA"
REDUCTION <- "harmony"

FRAC_EXPR <- 0.05
K_METACELL <- 25
MAX_SHARED <- 10
MIN_CELLS <- 100
NET_TYPE <- "signed"
SAVE_STAGE3_OBJECT <- TRUE

SUBTYPE_MAP <- c(
  "Proliferative-NPC" = "NPC-P",
  "OPC-Myelination" = "OPC-M",
  "Vascular-niche MES" = "MES-V",
  "MES-Antigen-presenting" = "MES-I"
)

FOCUS_GENES <- c(
  "PLAUR", "FOSL1", "FOSL2", "FOS", "FOSB",
  "JUN", "JUNB", "JUND", "BATF", "ATF3"
)
## ===========================================================================

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(OUT_DIR, "figures")
tab_dir <- file.path(OUT_DIR, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n########## STAGE 0: read object + clean genes ##########\n")
cat("R version:", R.version.string, "\n")
for (p in c("qs2", "Seurat", "SeuratObject", "WGCNA", "hdWGCNA", "harmony")) {
  cat(sprintf("  %-13s: %s\n", p, as.character(utils::packageVersion(p))))
}

stopifnot(file.exists(OBJ_PATH), file.exists(BL_PATH))
obj <- qs2::qs_read(OBJ_PATH)
DefaultAssay(obj) <- ASSAY

stopifnot(
  SUBTYPE_COL %in% colnames(obj@meta.data),
  PATIENT_COL %in% colnames(obj@meta.data),
  REDUCTION %in% Reductions(obj)
)

subtype_raw <- as.character(obj@meta.data[[SUBTYPE_COL]])
obj$subtype_bio <- unname(SUBTYPE_MAP[subtype_raw])
if (anyNA(obj$subtype_bio)) {
  stop("Unmapped subtype labels: ", paste(unique(subtype_raw[is.na(obj$subtype_bio)]), collapse = ", "))
}
obj$malignant_all <- "malignant"

cat("Cells:", ncol(obj), "\n")
cat("Genes:", nrow(obj), "\n")
cat("Subtype counts:\n")
print(table(obj$subtype_bio, useNA = "ifany"))
cat("Patient counts summary:\n")
print(summary(as.integer(table(obj@meta.data[[PATIENT_COL]]))))
cat("Patients with < MIN_CELLS:\n")
print(sort(table(obj@meta.data[[PATIENT_COL]])[table(obj@meta.data[[PATIENT_COL]]) < MIN_CELLS]))

# The blacklist file stores regex patterns, not literal gene symbols.
bl_pat <- readLines(BL_PATH, warn = FALSE)
bl_pat <- trimws(bl_pat)
bl_pat <- bl_pat[bl_pat != "" & !startsWith(bl_pat, "#")]
extra_pat <- "^(IGHD|IGHJ|IGKJ|IGLJ|TRBJ|TRAJ|TRGJ|TRDJ|TRDD)"
combined_pat <- paste(c(bl_pat, extra_pat), collapse = "|")

get_assay_data_compat <- function(seurat_obj, assay, layer) {
  tryCatch(
    SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = layer),
    error = function(e) {
      SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = layer)
    }
  )
}

cnts <- get_assay_data_compat(obj, assay = ASSAY, layer = "counts")
frac <- Matrix::rowSums(cnts > 0) / ncol(cnts)
genes_keep <- names(frac)[frac >= FRAC_EXPR]
n_before <- length(genes_keep)
bl_hit <- grep(combined_pat, genes_keep, value = TRUE)
genes_keep <- setdiff(genes_keep, bl_hit)

cat(sprintf(
  "Expression >= %.0f%% genes: %d -> blacklist regex hits: %d -> network genes: %d\n",
  100 * FRAC_EXPR, n_before, length(bl_hit), length(genes_keep)
))
cat("Blacklist hit examples:", paste(head(bl_hit, 20), collapse = ", "), "\n")
cat("Focus genes retained:", paste(intersect(FOCUS_GENES, genes_keep), collapse = ", "), "\n")

write.csv(
  data.frame(gene = genes_keep),
  file.path(tab_dir, "stage0_network_genes_after_blacklist.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(gene = bl_hit),
  file.path(tab_dir, "stage0_blacklist_regex_hits.csv"),
  row.names = FALSE
)

cat("\n########## STAGE 1: SetupForWGCNA + metacells by patient ##########\n")
obj <- SetupForWGCNA(
  seurat_obj = obj,
  wgcna_name = WGCNA_NAME,
  features = genes_keep
)

obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = c("malignant_all", PATIENT_COL),
  ident.group = "malignant_all",
  k = K_METACELL,
  reduction = REDUCTION,
  assay = ASSAY,
  slot = "counts",
  layer = "counts",
  min_cells = MIN_CELLS,
  max_shared = MAX_SHARED,
  wgcna_name = WGCNA_NAME,
  verbose = TRUE
)
obj <- NormalizeMetacells(obj, wgcna_name = WGCNA_NAME)

mc <- GetMetacellObject(obj, wgcna_name = WGCNA_NAME)
mc_pt <- table(mc@meta.data[[PATIENT_COL]])
cat("Metacell count:", ncol(mc), "\n")
cat("Effective patients:", length(mc_pt), "\n")
cat("Metacells per patient:\n")
print(sort(mc_pt))
write.csv(
  data.frame(Pt_number = names(mc_pt), metacells = as.integer(mc_pt)),
  file.path(tab_dir, "stage1_metacells_per_patient.csv"),
  row.names = FALSE
)

cat("\n########## STAGE 2: SetDatExpr, all malignant metacells ##########\n")
obj <- SetDatExpr(
  seurat_obj = obj,
  group_name = "malignant",
  group.by = "malignant_all",
  assay = ASSAY,
  slot = "data",
  layer = "data",
  wgcna_name = WGCNA_NAME
)

cat("\n########## STAGE 3: TestSoftPowers ##########\n")
obj <- TestSoftPowers(
  seurat_obj = obj,
  networkType = NET_TYPE,
  wgcna_name = WGCNA_NAME
)
pt_tbl <- GetPowerTable(obj, wgcna_name = WGCNA_NAME)
print(pt_tbl)
write.csv(pt_tbl, file.path(tab_dir, "softpower_table.csv"), row.names = FALSE)

pdf(file.path(fig_dir, "S_softpower_sweep.pdf"), width = 9, height = 7)
plot_list <- PlotSoftPowers(obj, wgcna_name = WGCNA_NAME)
print(patchwork::wrap_plots(plot_list, ncol = 2))
dev.off()

if (SAVE_STAGE3_OBJECT) {
  qs2::qs_save(
    obj,
    file.path(OUT_DIR, paste0("obj_", WGCNA_NAME, "_stage3_softpower.qs2"))
  )
}

cat("\n########## STOP 2 reached ##########\n")
cat("Report these files before choosing SOFT_POWER:\n")
cat("  ", file.path(tab_dir, "stage1_metacells_per_patient.csv"), "\n", sep = "")
cat("  ", file.path(tab_dir, "softpower_table.csv"), "\n", sep = "")
cat("  ", file.path(fig_dir, "S_softpower_sweep.pdf"), "\n", sep = "")
