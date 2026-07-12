# 05_恶性细胞分亚群与Neftel对照/09_prepare_reviewer_ai_handoff.R
# Prepare a read-only handoff package for external review of Neftel scoring,
# subtype status, and current malignant object metadata.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

set.seed(42)

params <- list(
  project_dir = "05_恶性细胞分亚群与Neftel对照",
  input_object = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.submodule_labeled.qs2"
  ),
  neftel_script = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "04_neftel_signature_scoring.R"
  ),
  submodule_script = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "04b_neftel_submodule_and_cycle_summary.R"
  ),
  tables_dir = file.path("05_恶性细胞分亚群与Neftel对照", "tables"),
  handoff_dir = file.path("05_恶性细胞分亚群与Neftel对照", "review_handoff"),
  subtype_col = "subtype_k4",
  sample_col = "Pt_number",
  provisional_labels = c(
    Subtype1 = "NPC2/Cycling-mixed",
    Subtype2 = "OPC-like",
    Subtype3 = "MES-Perivascular-like",
    Subtype4 = "MES-Inflammatory-like"
  )
)

dir.create(params$handoff_dir, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

write_session_info <- function(path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  sink(con)
  print(sessionInfo())
  sink()
}

summarise_score <- function(md, group_col, score_cols) {
  md |>
    select(all_of(c(group_col, score_cols))) |>
    pivot_longer(
      cols = all_of(score_cols),
      names_to = "score",
      values_to = "value"
    ) |>
    group_by(.data[[group_col]], score) |>
    summarise(
      n_cells = n(),
      n_non_missing = sum(!is.na(value)),
      n_missing = sum(is.na(value)),
      missing_pct = round(100 * n_missing / n_cells, 3),
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      q25 = quantile(value, 0.25, na.rm = TRUE, names = FALSE),
      q75 = quantile(value, 0.75, na.rm = TRUE, names = FALSE),
      iqr = IQR(value, na.rm = TRUE),
      min = min(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    rename(subtype_k4 = all_of(group_col))
}

read_script_excerpt <- function(path, n = 220) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  paste(utils::head(lines, n), collapse = "\n")
}

msg("Input object:", params$input_object)
stopifnot(file.exists(params$input_object))
stopifnot(file.exists(params$neftel_script))
stopifnot(file.exists(params$submodule_script))

obj <- qs2::qs_read(params$input_object)
md <- obj@meta.data

required_cols <- c(params$subtype_col, params$sample_col)
missing_required <- setdiff(required_cols, colnames(md))
if (length(missing_required) > 0) {
  stop("Missing metadata columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
}

score_cols <- c(
  "AMS_MES1", "AMS_MES2", "AMS_AC", "AMS_OPC", "AMS_NPC1", "AMS_NPC2", "AMS_G1S", "AMS_G2M",
  "MES1_UCell", "MES2_UCell", "AC_UCell", "OPC_UCell", "NPC1_UCell", "NPC2_UCell", "G1S_UCell", "G2M_UCell",
  "AMS_MES", "AMS_NPC", "AMS_cycling", "UCell_MES", "UCell_AC", "UCell_OPC", "UCell_NPC", "UCell_cycling"
)
score_cols <- intersect(score_cols, colnames(md))
if (length(score_cols) < 16) {
  stop("Too few Neftel/cell-cycle score columns found. Found: ", paste(score_cols, collapse = ", "), call. = FALSE)
}

subtype_levels <- paste0("Subtype", 1:4)
md[[params$subtype_col]] <- factor(as.character(md[[params$subtype_col]]), levels = subtype_levels)
md$provisional_label <- unname(params$provisional_labels[as.character(md[[params$subtype_col]])])

metadata_columns <- tibble(
  column_index = seq_along(colnames(md)),
  column_name = colnames(md)
)

sample_counts <- md |>
  count(.data[[params$sample_col]], name = "n_cells") |>
  arrange(.data[[params$sample_col]]) |>
  rename(sample = all_of(params$sample_col))

subtype_counts <- md |>
  count(.data[[params$subtype_col]], provisional_label, name = "n_cells") |>
  mutate(pct = round(100 * n_cells / sum(n_cells), 2)) |>
  arrange(.data[[params$subtype_col]])

sample_subtype_counts <- md |>
  count(.data[[params$sample_col]], .data[[params$subtype_col]], provisional_label, name = "n_cells") |>
  group_by(.data[[params$sample_col]]) |>
  mutate(sample_pct = round(100 * n_cells / sum(n_cells), 2)) |>
  ungroup() |>
  rename(sample = all_of(params$sample_col), subtype_k4 = all_of(params$subtype_col)) |>
  arrange(sample, subtype_k4)

score_summary <- summarise_score(md, params$subtype_col, score_cols) |>
  mutate(provisional_label = unname(params$provisional_labels[as.character(subtype_k4)])) |>
  relocate(provisional_label, .after = subtype_k4)

dominant_state_cols <- intersect(
  c("neftel_state_AMS", "neftel_state_UCell", "neftel_submodule_AMS", "neftel_submodule_UCell", "Phase"),
  colnames(md)
)

dominant_summary <- bind_rows(lapply(dominant_state_cols, function(col) {
  md |>
    count(.data[[params$subtype_col]], .data[[col]], name = "n_cells") |>
    group_by(.data[[params$subtype_col]]) |>
    mutate(row_pct = round(100 * n_cells / sum(n_cells), 2)) |>
    ungroup() |>
    transmute(
      subtype_k4 = .data[[params$subtype_col]],
      category = col,
      level = as.character(.data[[col]]),
      n_cells = n_cells,
      row_pct = row_pct
    )
}))

sanity_checks <- tibble(
  check = c(
    "input_object_exists",
    "n_cells",
    "n_features",
    "sample_col_present",
    "subtype_col_present",
    "score_columns_found",
    "subtype_na",
    "sample_na",
    "score_na_total",
    "provisional_label_na"
  ),
  value = c(
    file.exists(params$input_object),
    ncol(obj),
    nrow(obj),
    params$sample_col %in% colnames(md),
    params$subtype_col %in% colnames(md),
    length(score_cols),
    sum(is.na(md[[params$subtype_col]])),
    sum(is.na(md[[params$sample_col]])),
    sum(is.na(md[, score_cols])),
    sum(is.na(md$provisional_label))
  )
)

write_csv(metadata_columns, file.path(params$handoff_dir, "09_metadata_columns.csv"))
write_csv(sample_counts, file.path(params$handoff_dir, "09_cell_counts_by_sample.csv"))
write_csv(subtype_counts, file.path(params$handoff_dir, "09_cell_counts_by_subtype.csv"))
write_csv(sample_subtype_counts, file.path(params$handoff_dir, "09_cell_counts_by_sample_and_subtype.csv"))
write_csv(score_summary, file.path(params$handoff_dir, "09_neftel_cycle_score_summary_by_subtype.csv"))
write_csv(dominant_summary, file.path(params$handoff_dir, "09_dominant_state_summary_by_subtype.csv"))
write_csv(sanity_checks, file.path(params$handoff_dir, "09_handoff_sanity_checks.csv"))
write_session_info(file.path(params$handoff_dir, "09_session_info.txt"))

handoff_md <- c(
  "# Reviewer AI handoff: malignant subtype / Neftel status",
  "",
  "## Role boundary",
  "",
  "- This package is for method and claim review only.",
  "- Current subtype names are provisional. New outputs should use `Subtype1`-`Subtype4` as primary display labels.",
  "- `provisional_label` is provided as an explanatory metadata-style label, not as final nomenclature.",
  "- This script does not modify the Seurat object and does not call `RenameIdents()`.",
  "",
  "## Input",
  "",
  paste0("- Object: `", params$input_object, "`"),
  paste0("- Cells: ", ncol(obj)),
  paste0("- Features: ", nrow(obj)),
  paste0("- Sample column: `", params$sample_col, "`"),
  paste0("- Subtype column: `", params$subtype_col, "`"),
  "",
  "## Output files",
  "",
  "- `09_metadata_columns.csv`: all metadata columns in the inspected object.",
  "- `09_cell_counts_by_sample.csv`: malignant cells per patient/sample.",
  "- `09_cell_counts_by_subtype.csv`: malignant cells per Subtype1-4.",
  "- `09_cell_counts_by_sample_and_subtype.csv`: sample x subtype counts and within-sample percentages.",
  "- `09_neftel_cycle_score_summary_by_subtype.csv`: per-cell score summary, mean/median/IQR/range by subtype.",
  "- `09_dominant_state_summary_by_subtype.csv`: dominant Neftel state/submodule/Phase percentages by subtype.",
  "- `09_handoff_sanity_checks.csv`: input and NA sanity checks.",
  "- `09_session_info.txt`: R session info.",
  "",
  "## Provisional label map",
  "",
  paste0("- Subtype1: ", params$provisional_labels[["Subtype1"]]),
  paste0("- Subtype2: ", params$provisional_labels[["Subtype2"]]),
  paste0("- Subtype3: ", params$provisional_labels[["Subtype3"]]),
  paste0("- Subtype4: ", params$provisional_labels[["Subtype4"]]),
  "",
  "## Key completed status to review",
  "",
  "- scRNA QC, DoubletFinder, sample-wise inferCNV, immune-reference inferCNV rerun, and high-confidence malignant extraction have already been completed.",
  "- Malignant subtype work has completed malignant-only reclustering, theta sweep, k=4 candidate selection, Neftel AMS/UCell scoring, submodule/cell-cycle audit, subtype naming evidence, pathway enrichment, and patient composition sanity check.",
  "- Current reviewer-risk point: cluster-first subtype is not a one-to-one recapitulation of Neftel states. Subtype3/4 are both MES1-dominant and require marker/pathway/NMF support for separation.",
  "",
  "## Neftel scoring script excerpt",
  "",
  "```r",
  read_script_excerpt(params$neftel_script, n = 260),
  "```",
  "",
  "## Submodule/cell-cycle audit script excerpt",
  "",
  "```r",
  read_script_excerpt(params$submodule_script, n = 220),
  "```"
)

writeLines(handoff_md, file.path(params$handoff_dir, "09_reviewer_ai_handoff.md"), useBytes = TRUE)

msg("Sanity checks:")
print(sanity_checks)
msg("Subtype counts:")
print(subtype_counts)
msg("Outputs written to:", params$handoff_dir)
msg("Done.")
