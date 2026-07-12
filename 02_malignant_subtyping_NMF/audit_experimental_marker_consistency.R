suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
})

# Input:
#   Primary manuscript evidence:
#     <DATA_ROOT>/项目/分型/gene &dease/manuscript.docx
#   Current malignant Seurat object:
#     outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2
# Output:
#   tables/audit_experimental_markers.csv
#   tables/experimental_marker_consistency.csv
#   tables/audit_experimental_marker_session_info.txt
#
# Notes:
#   The manuscript was treated as the primary source. Wet-lab markers are mainly
#   from Figure 8 and the Methods paragraphs describing Western blotting and
#   PLAUR perturbation. Computational/spatial marker validations are retained
#   as a separate evidence tier, not mixed with wet-lab validation.

workdir <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/05_恶性细胞分亚群与Neftel对照"
setwd(workdir)

obj_path <- file.path("outputs", "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2")
out_marker <- file.path("tables", "audit_experimental_markers.csv")
out_consistency <- file.path("tables", "experimental_marker_consistency.csv")
out_session <- file.path("tables", "audit_experimental_marker_session_info.txt")

source_file <- "<DATA_ROOT>/项目/分型/gene &dease/manuscript.docx"

marker_tbl <- tibble::tribble(
  ~gene_symbol, ~display_name, ~experiment_type, ~evidence_tier, ~hypothesis_in_original, ~original_hypothesis_subtype, ~source_file, ~source_page_or_table, ~manuscript_evidence,
  "PLAUR", "PLAUR", "shRNA knockdown + Western blot + wound healing + Transwell + orthotopic xenograft", "wet_lab_primary", "PLAUR is the targetable effector of the invasive mesenchymal program; knockdown reverses EMT-like markers and suppresses migration/invasion/tumor growth.", "Subtype3", source_file, "Methods paragraphs 25-27; Results paragraphs 98-100; Figure 8I-N", "Stable PLAUR knockdown was performed in U-87 MG. Western blot validated PLAUR knockdown; wound healing, Transwell and xenograft assays tested functional consequences.",
  "VIM", "Vimentin", "Western blot", "wet_lab_primary", "Mesenchymal marker enriched in the invasive mesenchymal model and downregulated after PLAUR knockdown.", "Subtype3", source_file, "Methods paragraph 26; Results paragraphs 98-100; Figure 8E-F, 8J", "Western blot showed high VIM in U-87 MG and downregulation after shPLAUR.",
  "CDH2", "N-cadherin", "Western blot", "wet_lab_primary", "Mesenchymal/EMT marker enriched in the invasive mesenchymal model and downregulated after PLAUR knockdown.", "Subtype3", source_file, "Results paragraphs 99-100; Figure 8E-F, 8J", "Western blot showed high N-cadherin in U-87 MG and downregulation after shPLAUR.",
  "COL1A1", "COL1A1", "Western blot", "wet_lab_primary", "ECM component enriched in the invasive mesenchymal model and downregulated after PLAUR knockdown.", "Subtype3", source_file, "Results paragraphs 99-100; Figure 8E-F, 8J", "Western blot showed high COL1A1 in U-87 MG and downregulation after shPLAUR.",
  "CDH1", "E-cadherin", "Western blot", "wet_lab_primary", "Epithelial/MET marker upregulated after PLAUR knockdown; not expected to define a malignant subtype in scRNA-seq.", NA_character_, source_file, "Methods paragraph 26; Results paragraph 100; Figure 8J", "E-cadherin was used as a MET-response marker after PLAUR knockdown.",
  "TOP2A", "TOP2A", "Western blot", "wet_lab_primary", "Proliferation marker higher in T98G than U-87 MG in cell-line comparison.", "Subtype1", source_file, "Results paragraph 99; Figure 8E-F", "Western blot showed T98G overexpressed TOP2A relative to U-87 MG.",
  "PLAU", "PLAU", "cell-line scRNA signature + CellChat/spatial interaction", "computational_functional", "Ligand component of the PLAU-PLAUR axis linked to mesenchymal-macrophage communication.", "Subtype3", source_file, "Paragraphs 90, 96, 98, 109; Figure 7F-G; Figure 8C/O", "PLAU was part of the invasive-mesenchymal signature and the PLAU-PLAUR interaction model.",
  "MMP9", "MMP9", "cell-line scRNA signature + schematic mechanism", "computational_functional", "Invasive-mesenchymal signature gene downstream of FOSL1/PLAUR axis.", "Subtype3", source_file, "Results paragraphs 98, 101, 110; Figure 8C/O", "MMP9 was included in the invasive-mesenchymal signature and mechanistic schematic.",
  "MMP14", "MMP14", "cell-line scRNA signature", "computational_functional", "Invasive-mesenchymal signature gene used for cell-line model selection.", "Subtype3", source_file, "Results paragraph 98; Figure 8C", "MMP14 was included in the invasive-mesenchymal signature for model selection.",
  "FOSL1", "FOSL1", "scATAC/SCENIC/spatial regulon evidence", "computational_regulatory", "Regulatory gatekeeper of the mesenchymal program and PLAU-PLAUR axis.", "Subtype3", source_file, "Paragraphs 54-61, 80-84, 101, 108; Figures 3, 6, 8O", "FOSL1 regulon/motif/spatial activity was linked to mesenchymal adaptation; not a wet-lab marker in the manuscript.",
  "GAP43", "GAP43", "scRNA/spatial feature validation", "computational_marker", "Astro-like/root marker in the original four-subtype framework.", NA_character_, source_file, "Paragraphs 44, 46-47; Figure S2C-E", "GAP43 was validated by feature/violin/spatial plots in the original computational framework.",
  "MAG", "MAG", "scRNA/spatial feature validation", "computational_marker", "Oligodendrocyte-like marker in the original four-subtype framework.", "Subtype2", source_file, "Paragraphs 42, 44, 47; Figure S2C-E", "MAG was validated by feature/violin/spatial plots in the original computational framework."
)

readr::write_csv(marker_tbl, out_marker)

obj <- qs2::qs_read(obj_path)
DefaultAssay(obj) <- "RNA"
md <- obj@meta.data
stopifnot("subtype_k4" %in% colnames(md))

subtype_levels <- c("Subtype1", "Subtype2", "Subtype3", "Subtype4")
md$subtype_k4 <- factor(as.character(md$subtype_k4), levels = subtype_levels)
stopifnot(!any(is.na(md$subtype_k4)))

genes <- unique(marker_tbl$gene_symbol)
genes_present <- intersect(genes, rownames(obj))
genes_missing <- setdiff(genes, rownames(obj))

expr <- GetAssayData(obj, assay = "RNA", layer = "data")
mean_by_subtype <- matrix(NA_real_, nrow = length(genes), ncol = length(subtype_levels),
                          dimnames = list(genes, subtype_levels))

for (g in genes_present) {
  vals <- expr[g, , drop = TRUE]
  for (s in subtype_levels) {
    cells <- which(md$subtype_k4 == s)
    mean_by_subtype[g, s] <- Matrix::mean(vals[cells])
  }
}

consistency <- marker_tbl %>%
  mutate(
    gene_present = gene_symbol %in% genes_present,
    mean_Subtype1 = mean_by_subtype[gene_symbol, "Subtype1"],
    mean_Subtype2 = mean_by_subtype[gene_symbol, "Subtype2"],
    mean_Subtype3 = mean_by_subtype[gene_symbol, "Subtype3"],
    mean_Subtype4 = mean_by_subtype[gene_symbol, "Subtype4"],
    new_top_subtype = apply(mean_by_subtype[gene_symbol, subtype_levels, drop = FALSE], 1, function(x) {
      if (all(is.na(x))) return(NA_character_)
      subtype_levels[which.max(x)]
    }),
    consistent = dplyr::case_when(
      is.na(original_hypothesis_subtype) ~ NA,
      !gene_present ~ NA,
      TRUE ~ new_top_subtype == original_hypothesis_subtype
    ),
    consistency_note = dplyr::case_when(
      is.na(original_hypothesis_subtype) ~ "No subtype-level consistency expected; marker was used for perturbation/MET response or old framework only.",
      !gene_present ~ "Gene not detected in current object.",
      consistent ~ "Top mean expression matches the original hypothesis subtype.",
      TRUE ~ "Top mean expression does not match original hypothesis; inspect before integrating into figures."
    )
  )

readr::write_csv(consistency, out_consistency)
writeLines(capture.output(sessionInfo()), out_session)

cat("Audit marker rows:", nrow(marker_tbl), "\n")
cat("Genes present:", length(genes_present), "of", length(genes), "\n")
if (length(genes_missing) > 0) {
  cat("Genes missing:", paste(genes_missing, collapse = ", "), "\n")
}
cat("Evidence tiers:\n")
print(table(marker_tbl$evidence_tier))
cat("Consistency by tier:\n")
print(with(consistency, table(evidence_tier, consistent, useNA = "ifany")))
cat("Outputs:\n")
cat(out_marker, "\n")
cat(out_consistency, "\n")
