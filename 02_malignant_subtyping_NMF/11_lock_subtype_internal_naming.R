suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(dplyr)
  library(readr)
})

# Lock internal biological names after non-malignant signature audit.
# Manuscript-facing figures must still use Subtype1-4.

step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir) {
  setwd(file.path(getwd(), step_dir))
}

obj_path <- file.path("outputs", "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2")
label_csv <- file.path("tables", "11_subtype_internal_naming_locked.csv")
sanity_csv <- file.path("tables", "11_subtype_internal_naming_sanity_checks.csv")
session_txt <- file.path("tables", "11_subtype_internal_naming_session_info.txt")

label_map <- tibble::tribble(
  ~subtype_k4, ~subtype_label_final, ~rationale, ~mandatory_disclaimer,
  "Subtype1", "Proliferative-NPC", "NPC + proliferation anchor; NMF cycling program support.", "Internal biological label only; manuscript-facing figures use Subtype1.",
  "Subtype2", "OPC-Myelination", "OPC/myelination anchor; MAG rank 11 in original-vs-revised remapping audit.", "Internal biological label only; manuscript-facing figures use Subtype2.",
  "Subtype3", "Vascular-niche MES", "Niche-coupled mesenchymal program with pericyte/vSMC-mimicking signature; 100% inferCNV high-confidence malignant.", "Vascular-niche refers to malignant transcriptional coupling to perivascular cues, not vascular cell identity.",
  "Subtype4", "MES-Antigen-presenting", "MHC-II antigen-presentation malignant program supported by NMF MP04 and GO immune/cytokine terms.", "MES-Antigen-presenting denotes MHC class II expression program in CNV-positive malignant cells, not professional antigen-presenting cell identity."
)

obj <- qs2::qs_read(obj_path)
stopifnot("subtype_k4" %in% colnames(obj@meta.data))
obj@meta.data$subtype_k4 <- as.character(obj@meta.data$subtype_k4)
stopifnot(setequal(unique(obj@meta.data$subtype_k4), label_map$subtype_k4))

map_vec <- setNames(label_map$subtype_label_final, label_map$subtype_k4)
obj@meta.data$subtype_label_final <- unname(map_vec[obj@meta.data$subtype_k4])
stopifnot(!any(is.na(obj@meta.data$subtype_label_final)))

qs2::qs_save(obj, obj_path)
readr::write_csv(label_map, label_csv)

sanity <- obj@meta.data %>%
  as_tibble() %>%
  count(subtype_k4, subtype_label_final, name = "n_cells") %>%
  arrange(subtype_k4)
readr::write_csv(sanity, sanity_csv)
writeLines(capture.output(sessionInfo()), session_txt)

cat("Updated object:", obj_path, "\n")
cat("Label map:", label_csv, "\n")
cat("Sanity:\n")
print(sanity)
