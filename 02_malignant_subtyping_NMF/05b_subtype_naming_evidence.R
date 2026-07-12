# 05_恶性细胞分亚群与Neftel对照/05b_subtype_naming_evidence.R
# Multidimensional evidence for final subtype naming.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(42)

proj <- "05_恶性细胞分亚群与Neftel对照"
in_obj <- file.path(proj, "outputs", "GBM.malignant.subtyped.umap_candidates.qs2")
tab_dir <- file.path(proj, "tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) {
    stop("None of these columns exist: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  hit[1]
}

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

msg("Loading object:", in_obj)
obj <- qs2::qs_read(in_obj)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

md <- obj@meta.data
required <- c(
  "subtype_k4",
  "UCell_MES", "UCell_AC", "UCell_OPC", "UCell_NPC", "UCell_cycling",
  "MES1_UCell", "MES2_UCell",
  "Phase"
)
missing <- setdiff(required, colnames(md))
if (length(missing) > 0) {
  stop("Missing metadata columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

obj$subtype_k4 <- factor(as.character(obj$subtype_k4), levels = c("Subtype1", "Subtype2", "Subtype3", "Subtype4"))
md <- obj@meta.data

# 1. Multidimensional subtype evidence.
ev <- md |>
  group_by(subtype_k4) |>
  summarise(
    n_cells = n(),
    MES_mean = mean(UCell_MES),
    AC_mean = mean(UCell_AC),
    OPC_mean = mean(UCell_OPC),
    NPC_mean = mean(UCell_NPC),
    cycling_mean = mean(UCell_cycling),
    pct_G1 = mean(Phase == "G1") * 100,
    pct_G2M = mean(Phase == "G2M") * 100,
    pct_S = mean(Phase == "S") * 100,
    MES1_mean = mean(MES1_UCell),
    MES2_mean = mean(MES2_UCell),
    .groups = "drop"
  )

write.csv(ev, file.path(tab_dir, "05b_subtype_multimodal_evidence.csv"), row.names = FALSE)
cat("\n== Subtype multimodal evidence ==\n")
print(ev)

# 2. Full positive markers for each subtype.
msg("Running FindAllMarkers for subtype_k4.")
Idents(obj) <- "subtype_k4"
mk <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5
)

fc_col <- pick_col(mk, c("avg_log2FC", "avg_logFC"))
write.csv(mk, file.path(tab_dir, "05b_subtype_markers_all.csv"), row.names = FALSE)

top10 <- mk |>
  group_by(cluster) |>
  slice_max(.data[[fc_col]], n = 10, with_ties = FALSE) |>
  ungroup()
write.csv(top10, file.path(tab_dir, "05b_subtype_markers_top10.csv"), row.names = FALSE)

cat("\n== Top 10 markers per subtype ==\n")
print(top10)

# 3. Naming confidence table. Internal reference, not manuscript text.
confidence <- data.frame(
  subtype = c("Subtype1", "Subtype2", "Subtype3", "Subtype4"),
  proposed_name = c("NPC-Cycling", "OPC-like", "MES-Perivascular", "MES-Inflammatory"),
  neftel_anchor = c("NPC2+cycling", "OPC", "MES (refined)", "MES (refined)"),
  evidence_strength = c("medium", "high", "high", "high"),
  notes = c(
    "Mixed NPC2/MES1 dominant plus high S/G2M fraction; cycling is a defining axis.",
    "OPC UCell dominant and almost entirely G1; cleanest Neftel match.",
    "ACTA2/TAGLN/FRZB/MFAP4-like signal suggests perivascular/vascular-mesenchymal program.",
    "C1Q/CD68/APOC1/IL1B-like signal suggests inflammatory/myeloid-like MES program in malignant cells."
  ),
  stringsAsFactors = FALSE
)

write.csv(confidence, file.path(tab_dir, "05b_naming_confidence.csv"), row.names = FALSE)
cat("\n== Naming confidence ==\n")
print(confidence)

cat("\nWritten files:\n")
cat("- tables/05b_subtype_multimodal_evidence.csv\n")
cat("- tables/05b_subtype_markers_all.csv\n")
cat("- tables/05b_subtype_markers_top10.csv\n")
cat("- tables/05b_naming_confidence.csv\n")
msg("Done.")
