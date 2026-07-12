# Extract high-confidence malignant cells for malignant-cell reclustering.
# This script is deterministic: no subtype decision is made here.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(42)

config <- list(
  input_qs2 = normalizePath(
    file.path("04_inferCNV_免疫参考验证", "outputs", "GBM.RNA.qc_doubletfinder.infercnv_immune_reference_calls.qs2"),
    winslash = "\\",
    mustWork = FALSE
  ),
  out_dir = normalizePath(file.path("05_恶性细胞分亚群与Neftel对照", "outputs"), winslash = "\\", mustWork = FALSE),
  tab_dir = normalizePath(file.path("05_恶性细胞分亚群与Neftel对照", "tables"), winslash = "\\", mustWork = FALSE),
  fig_dir = normalizePath(file.path("05_恶性细胞分亚群与Neftel对照", "figures"), winslash = "\\", mustWork = FALSE),
  sample_col = "Pt_number",
  annotation_col = "anno_ident",
  malignant_bool_col = "is_malignant_for_downstream_immune",
  call_col = "infercnv_immune_call",
  status_col = "malignant_call_status_immune",
  min_malignant_cells_per_patient = 50
)

dir.create(config$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(config$fig_dir, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

if (!file.exists(config$input_qs2)) {
  stop("Input object not found: ", config$input_qs2, call. = FALSE)
}

msg("Loading:", config$input_qs2)
obj <- qs2::qs_read(config$input_qs2)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

msg("Input object:", ncol(obj), "cells x", nrow(obj), "genes")

required_cols <- c(
  config$sample_col,
  config$annotation_col,
  config$malignant_bool_col,
  config$call_col,
  config$status_col
)
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

metadata_snapshot <- data.frame(
  column = colnames(obj@meta.data),
  class = vapply(obj@meta.data, function(x) paste(class(x), collapse = ";"), character(1)),
  n_unique = vapply(obj@meta.data, function(x) length(unique(as.character(x))), integer(1)),
  stringsAsFactors = FALSE
)
write.csv(
  metadata_snapshot,
  file.path(config$tab_dir, "01_input_metadata_columns_snapshot.csv"),
  row.names = FALSE
)

call_summary <- obj@meta.data |>
  count(.data[[config$call_col]], name = "n_cells") |>
  mutate(percent = round(100 * n_cells / sum(n_cells), 3)) |>
  rename(infercnv_immune_call = 1)
write.csv(
  call_summary,
  file.path(config$tab_dir, "01_infercnv_immune_call_overall.csv"),
  row.names = FALSE
)

status_summary <- obj@meta.data |>
  count(.data[[config$status_col]], name = "n_cells") |>
  mutate(percent = round(100 * n_cells / sum(n_cells), 3)) |>
  rename(malignant_call_status_immune = 1)
write.csv(
  status_summary,
  file.path(config$tab_dir, "01_malignant_call_status_overall.csv"),
  row.names = FALSE
)

status_by_sample <- obj@meta.data |>
  count(.data[[config$sample_col]], .data[[config$status_col]], name = "n_cells") |>
  rename(Pt_number = 1, malignant_call_status_immune = 2) |>
  tidyr::pivot_wider(
    names_from = malignant_call_status_immune,
    values_from = n_cells,
    values_fill = 0
  )
write.csv(
  status_by_sample,
  file.path(config$tab_dir, "01_malignant_call_status_by_sample.csv"),
  row.names = FALSE
)

msg("Malignant status distribution:")
print(status_summary)

main_cells <- colnames(obj)[which(obj@meta.data[[config$malignant_bool_col]] %in% TRUE)]
if (length(main_cells) < 1000) {
  warning("Only ", length(main_cells), " high-confidence malignant cells found. Check upstream calls.")
}

malig <- subset(obj, cells = main_cells)
malig@meta.data[[config$sample_col]] <- droplevels(factor(malig@meta.data[[config$sample_col]]))

msg("Extracted high-confidence malignant cells:", ncol(malig))
msg("Patients represented:", nlevels(malig@meta.data[[config$sample_col]]))

per_patient <- malig@meta.data |>
  count(.data[[config$sample_col]], name = "n_malignant") |>
  rename(Pt_number = 1) |>
  arrange(n_malignant)
write.csv(
  per_patient,
  file.path(config$tab_dir, "01_high_confidence_malignant_cells_per_patient.csv"),
  row.names = FALSE
)

low_n_patients <- per_patient$Pt_number[
  per_patient$n_malignant < config$min_malignant_cells_per_patient
]
low_n_table <- data.frame(
  Pt_number = low_n_patients,
  n_malignant = per_patient$n_malignant[
    per_patient$Pt_number %in% low_n_patients
  ],
  threshold = config$min_malignant_cells_per_patient,
  stringsAsFactors = FALSE
)
write.csv(
  low_n_table,
  file.path(config$tab_dir, "01_patients_below_50_malignant_cells.csv"),
  row.names = FALSE
)
if (length(low_n_patients) > 0) {
  warning(
    "Patients with < ", config$min_malignant_cells_per_patient,
    " high-confidence malignant cells: ",
    paste(low_n_patients, collapse = ", ")
  )
}

sens_cells <- colnames(obj)[
  obj@meta.data[[config$call_col]] %in% c(
    "malignant_like_CNV_high_confidence",
    "malignant_like_CNV_burden_only"
  )
]
malig_sens <- subset(obj, cells = sens_cells)
malig_sens@meta.data[[config$sample_col]] <- droplevels(factor(malig_sens@meta.data[[config$sample_col]]))
msg("Sensitivity subset, high-confidence + burden-only:", ncol(malig_sens), "cells")

main_out <- file.path(config$out_dir, "GBM.malignant.high_confidence.qs2")
sens_out <- file.path(config$out_dir, "GBM.malignant.sensitivity_with_burden_only.qs2")
qs2::qs_save(malig, main_out)
qs2::qs_save(malig_sens, sens_out)

plot_data <- per_patient |>
  mutate(
    below_threshold = n_malignant < config$min_malignant_cells_per_patient,
    Pt_number = factor(Pt_number, levels = Pt_number)
  )

p <- ggplot(plot_data, aes(x = Pt_number, y = n_malignant, fill = below_threshold)) +
  geom_col(width = 0.78) +
  geom_hline(
    yintercept = config$min_malignant_cells_per_patient,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.35
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c("FALSE" = "#B24745", "TRUE" = "#9A9A9A"),
    labels = c("FALSE" = "Included", "TRUE" = "Below 50"),
    name = NULL
  ) +
  labs(
    x = "Patient",
    y = "High-confidence malignant cells",
    title = "High-confidence malignant cells per patient",
    subtitle = "Dashed line: n = 50 threshold for downstream inclusion checks"
  ) +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey35")
  )

ggsave(
  file.path(config$fig_dir, "01_high_confidence_malignant_cells_per_patient.pdf"),
  p,
  width = 6.2,
  height = 7,
  device = grDevices::cairo_pdf
)

msg("Done.")
msg("Main object:", main_out)
msg("Sensitivity object:", sens_out)
