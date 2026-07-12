# 05_恶性细胞分亚群与Neftel对照/10e_apply_metaprogram_final_labels.R
# Apply final ORA-reviewed metaprogram labels to a new qs2 Seurat object.
# Does not overwrite the original v2 object.

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

set.seed(42)

params <- list(
  input_object = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.v2.qs2"
  ),
  metaprogram_scores_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10c_per_cell_metaprogram_scores.tsv"
  ),
  proposed_label_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10b_metaprogram_label_proposed.csv"
  ),
  final_label_file = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10b_metaprogram_label_final.csv"
  ),
  output_object = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.v2.metaprogram_labeled.qs2"
  ),
  metadata_out = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10e_metaprogram_final_label_metadata.csv"
  ),
  session_out = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "tables",
    "10e_apply_metaprogram_final_labels_session_info.txt"
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

final_labels <- tibble::tribble(
  ~metaprogram_id, ~metaprogram_label_final, ~label_status, ~label_basis, ~notes,
  "MP01", "IGFBP-signaling / Stress-response program", "final", "ORA main plus Evidence-3 convergence", "Final metaprogram label; keep MP ID in figures.",
  "MP02", "ECM-organization program", "final", "ORA main: EMT, ECM proteoglycans, integrin, extracellular matrix organization; S3 GSEA EMT enriched", "Pathway/process label.",
  "MP03", "Myelination program", "final", "ORA main: ensheathment of neurons, axon development, gliogenesis; S2 pathway overview neuronal system", "Pathway/process label.",
  "MP04", "MHC-II Antigen-presentation program", "final", "ORA main: MHC-II antigen processing/presentation; S4 GSEA inflammatory/TNFA/IFN/allograft enriched", "Pathway/process label; MP04 audit retained with biology note.",
  "MP05", "Cell-cycle program", "final", "ORA main: G2M checkpoint, E2F targets, mitotic cell cycle; S1 GSEA G2M/E2F enriched", "Cycling supplement label.",
  "MP06", "Angiogenesis-signaling program", "final", "FLT1, KDR, ANGPT2 and vascular/angiogenesis-related ORA; renamed from Vasculature-development after review", "Pathway/process label; MP06 audit retained with biology note."
)
write_csv(final_labels, params$final_label_file)

proposed_labels <- read_csv(params$proposed_label_file, show_col_types = FALSE) |>
  select(metaprogram_id, metaprogram_label_proposed)

msg("Loading object")
obj <- qs2::qs_read(params$input_object)
md <- obj@meta.data
md$cell_id <- rownames(md)

scores <- read_tsv(params$metaprogram_scores_file, show_col_types = FALSE)
score_cols_raw <- c("MP01", "MP02", "MP03", "MP04", "MP05_cycling_supplement", "MP06")
stopifnot(all(c("cell_id", score_cols_raw) %in% colnames(scores)))
score_df <- scores |>
  select(cell_id, all_of(score_cols_raw)) |>
  rename(MP05 = MP05_cycling_supplement)
score_cols <- c("MP01", "MP02", "MP03", "MP04", "MP05", "MP06")

stopifnot(identical(sort(md$cell_id), sort(score_df$cell_id)))
mp_scores <- as.matrix(score_df[, score_cols])
rownames(mp_scores) <- score_df$cell_id
dominant_id <- score_cols[max.col(mp_scores, ties.method = "first")]
dominant_df <- tibble(
  cell_id = rownames(mp_scores),
  metaprogram_id_final = dominant_id,
  metaprogram_dominant_score = mp_scores[cbind(seq_len(nrow(mp_scores)), match(dominant_id, score_cols))]
) |>
  left_join(proposed_labels, by = c("metaprogram_id_final" = "metaprogram_id")) |>
  left_join(final_labels |> select(metaprogram_id, metaprogram_label_final), by = c("metaprogram_id_final" = "metaprogram_id"))

score_out <- score_df |>
  rename(
    MP01_UCell = MP01,
    MP02_UCell = MP02,
    MP03_UCell = MP03,
    MP04_UCell = MP04,
    MP05_UCell = MP05,
    MP06_UCell = MP06
  )

new_md <- md |>
  left_join(score_out, by = "cell_id") |>
  left_join(dominant_df, by = "cell_id")
rownames(new_md) <- new_md$cell_id
new_md$cell_id <- NULL
stopifnot(identical(rownames(obj@meta.data), rownames(new_md)))
stopifnot(!any(is.na(new_md$metaprogram_label_final)))
obj@meta.data <- new_md

write_csv(
  tibble(cell_id = rownames(new_md)) |>
    bind_cols(new_md[, c(
      "subtype_k4",
      "metaprogram_id_final",
      "metaprogram_label_proposed",
      "metaprogram_label_final",
      "metaprogram_dominant_score"
    )]),
  params$metadata_out
)

qs2::qs_save(obj, params$output_object)
write_session_info(params$session_out)

sanity <- tibble(
  n_cells = ncol(obj),
  n_final_labels = n_distinct(obj@meta.data$metaprogram_label_final),
  output_object = params$output_object,
  n_missing_final_label = sum(is.na(obj@meta.data$metaprogram_label_final))
)
write_csv(
  sanity,
  file.path("05_恶性细胞分亚群与Neftel对照", "tables", "10e_apply_metaprogram_final_labels_sanity_checks.csv")
)

msg("Applied final metaprogram labels")
print(sanity)
print(table(obj@meta.data$metaprogram_id_final))
