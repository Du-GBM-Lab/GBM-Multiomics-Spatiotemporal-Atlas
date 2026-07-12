# 05_恶性细胞分亚群与Neftel对照/10f_apply_subtype_final_labels.R
# Apply final subtype labels after Rule v2.2 merge decision.
# Does not overwrite previous qs2 objects.

Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS = "1")
Sys.setenv(MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(qs2)
  library(Seurat)
  library(dplyr)
  library(readr)
})

params <- list(
  input_object = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.v2.metaprogram_labeled.qs2"
  ),
  output_object = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"
  ),
  mapping_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10f_subtype_label_final.csv"
  ),
  metadata_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10f_subtype_final_label_metadata.csv"
  ),
  sanity_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10f_apply_subtype_final_labels_sanity_checks.csv"
  ),
  session_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10f_apply_subtype_final_labels_session_info.txt"
  )
)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

write_session_info <- function(path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  sink(con)
  print(sessionInfo())
  sink()
}

subtype_map <- tibble::tribble(
  ~subtype_k4, ~subtype_label_final, ~label_status, ~naming_rationale, ~evidence_sources,
  "Subtype1", "Proliferative-NPC subtype", "final", "Cycling signal is dominant while NPC2/NPC-like signal is retained; final name places proliferation first to avoid overclaiming a pure NPC identity.", "Neftel G1S/G2M and NPC2 scores; MP05 Cell-cycle program; downsampled one-vs-rest GSEA G2M/E2F/MYC/DNA repair.",
  "Subtype2", "OPC-Myelination subtype", "final", "OPC-like state is supported by a dominant myelination-related metaprogram and neuronal/myelination pathway overview.", "Neftel OPC score; MP03 Myelination program; one-vs-rest GSEA neuronal-system / synapse / ion-channel terms.",
  "Subtype3", "MES-ECM subtype", "final", "MES-like subtype distinguished by ECM organization and angiogenesis-signaling programs; concise label uses ECM as the strongest and most reviewer-safe process term.", "MP02 ECM-organization program; MP06 Angiogenesis-signaling program; S3 vs S4 GSEA EMT/OXPHOS/myogenesis/hypoxia; marker heatmap S3-up block.",
  "Subtype4", "MES-Antigen-presenting subtype", "final", "MES-like subtype distinguished by MHC-II antigen-presentation and inflammatory signaling; label uses the most specific ORA-supported process.", "MP04 MHC-II Antigen-presentation program; S3 vs S4 GSEA TNFA/IFN/allograft/complement; marker heatmap S4-up block."
)

write_csv(subtype_map, params$mapping_file)

msg("Loading object:", params$input_object)
obj <- qs2::qs_read(params$input_object)
md <- obj@meta.data
stopifnot("subtype_k4" %in% colnames(md))

md$cell_id <- rownames(md)
new_md <- md |>
  left_join(subtype_map |> select(subtype_k4, subtype_label_final), by = "subtype_k4")
rownames(new_md) <- new_md$cell_id
new_md$cell_id <- NULL
stopifnot(identical(rownames(obj@meta.data), rownames(new_md)))
stopifnot(!any(is.na(new_md$subtype_label_final)))
obj@meta.data <- new_md

write_csv(
  tibble(cell_id = rownames(obj@meta.data)) |>
    bind_cols(obj@meta.data[, c("subtype_k4", "subtype_label_final", "metaprogram_id_final", "metaprogram_label_final")]),
  params$metadata_file
)

qs2::qs_save(obj, params$output_object)

sanity <- tibble(
  n_cells = ncol(obj),
  n_subtype_ids = n_distinct(obj@meta.data$subtype_k4),
  n_subtype_final_labels = n_distinct(obj@meta.data$subtype_label_final),
  n_missing_subtype_label_final = sum(is.na(obj@meta.data$subtype_label_final)),
  output_object = params$output_object
)
write_csv(sanity, params$sanity_file)
write_session_info(params$session_file)

msg("Applied final subtype labels")
print(sanity)
print(table(obj@meta.data$subtype_k4, obj@meta.data$subtype_label_final))
