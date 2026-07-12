options(stringsAsFactors = FALSE)
set.seed(42)

if (basename(getwd()) == "06_恶性细胞拟时序") {
  project_root <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(file.path(project_root, "06_恶性细胞拟时序"))

dir.create("tables", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/source_data", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

log_file <- file.path("logs", "07_validation.log")
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  txt <- paste(...)
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", txt)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_msg("Phase 3 supplemental terminal marker validation started")
log_msg("set.seed = 42")

suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

source("scripts/_naming.R")

final_obj_path <- file.path(project_root, "05_恶性细胞分亚群与Neftel对照", "outputs", "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2")
full_obj_path <- file.path(project_root, "GBM.RNA.integrated.24.rds")
panel_d_path <- file.path("figures/source_data", "06_panel_D_endpoint_boxplot.csv")
raw_endpoint_path <- file.path("tables", "terminal_endpoint_composition_pseudotime_top20_raw.csv")

log_msg("Final malignant object:", final_obj_path)
log_msg("Full object:", full_obj_path)
log_msg("Endpoint source:", panel_d_path)

get_data_matrix <- function(obj, genes) {
  genes <- intersect(genes, rownames(obj[["RNA"]]))
  if (length(genes) == 0) return(matrix(nrow = 0, ncol = ncol(obj)))
  mat <- tryCatch(
    Seurat::GetAssayData(obj, assay = "RNA", layer = "data")[genes, , drop = FALSE],
    error = function(e) Seurat::GetAssayData(obj, assay = "RNA", slot = "data")[genes, , drop = FALSE]
  )
  mat
}

summarise_expr <- function(df, group_col = "validation_group") {
  df |>
    tidyr::pivot_longer(cols = dplyr::all_of(marker_genes), names_to = "marker", values_to = "expression") |>
    dplyr::group_by(.data[[group_col]], .data$marker) |>
    dplyr::summarise(
      n_cells = dplyr::n(),
      mean = mean(.data$expression, na.rm = TRUE),
      median = stats::median(.data$expression, na.rm = TRUE),
      pct_pos = mean(.data$expression > 0, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::rename(validation_group = dplyr::all_of(group_col))
}

save_plot <- function(plot, stem, width = 5.2, height = 3.6) {
  pdf_path <- file.path("figures", paste0(stem, ".pdf"))
  png_path <- file.path("figures", paste0(stem, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf)
  ggsave(png_path, plot, width = width, height = height, dpi = 300)
  log_msg("Wrote", pdf_path)
  log_msg("Wrote", png_path)
}

mhc_genes <- c("CD74", "HLA-DRA", "HLA-DPB1")
vascular_genes <- c("RGS5", "ACTA2", "TAGLN", "HIGD1B")
malignant_marker_genes <- c("SOX2", "SOX4", "SOX9", "OLIG1", "OLIG2")
marker_genes <- unique(c(mhc_genes, vascular_genes, malignant_marker_genes))

load_start <- Sys.time()
obj <- qs2::qs_read(final_obj_path)
log_msg("Final object load seconds:", round(as.numeric(difftime(Sys.time(), load_start, units = "secs")), 2))

meta <- obj@meta.data |>
  tibble::rownames_to_column("cellID") |>
  recode_subtype()

endpoint_source <- readr::read_csv(panel_d_path, show_col_types = FALSE) |>
  dplyr::distinct(.data$cellID, .data$endpoint_group, .data$subtype_k4, .data$subtype_label_final, .data$subtype_short, .data$lineage_id, .data$pseudotime)

endpoint_cells <- endpoint_source |>
  dplyr::filter(
    (.data$endpoint_group == "MES-I terminal" & .data$subtype_short == "MES-I") |
      (.data$endpoint_group == "MES-V terminal" & .data$subtype_short == "MES-V")
  )

if (nrow(endpoint_cells) == 0) stop("No strict endpoint cells found from Panel D source data.")
log_msg("Strict endpoint cells:", nrow(endpoint_cells))

all_mes_cells <- meta |>
  dplyr::filter(.data$subtype_short %in% c("MES-I", "MES-V")) |>
  dplyr::transmute(
    cellID,
    endpoint_group = paste0(as.character(.data$subtype_short), " all"),
    subtype_k4,
    subtype_label_final,
    subtype_short,
    lineage_id = NA_character_,
    pseudotime = NA_real_
  )

expr_mat <- get_data_matrix(obj, marker_genes)
available_markers <- rownames(expr_mat)
missing_markers <- setdiff(marker_genes, available_markers)
if (length(missing_markers) > 0) log_msg("Missing markers in malignant object:", paste(missing_markers, collapse = ", "))
expr_df <- as.data.frame(t(as.matrix(expr_mat))) |>
  tibble::rownames_to_column("cellID")
for (g in setdiff(marker_genes, colnames(expr_df))) expr_df[[g]] <- NA_real_

validation_cells <- dplyr::bind_rows(
  endpoint_cells |> dplyr::mutate(validation_group = as.character(.data$endpoint_group)),
  all_mes_cells |> dplyr::mutate(validation_group = as.character(.data$endpoint_group))
) |>
  dplyr::left_join(meta, by = c("cellID", "subtype_k4", "subtype_label_final", "subtype_short")) |>
  dplyr::left_join(expr_df, by = "cellID") |>
  dplyr::mutate(
    is_cnv_positive = .data$is_malignant_for_downstream_immune %in% TRUE | .data$malignant_call_status_immune %in% c("malignant", "high_confidence_malignant"),
    is_doublet = .data$doubletfinder_class %in% c("Doublet", "doublet")
  )

mhc_present <- intersect(mhc_genes, colnames(validation_cells))
vascular_present <- intersect(vascular_genes, colnames(validation_cells))
malignant_marker_present <- intersect(malignant_marker_genes, colnames(validation_cells))
validation_cells$any_mhc_pos <- if (length(mhc_present) > 0) rowSums(validation_cells[, mhc_present, drop = FALSE] > 0, na.rm = TRUE) > 0 else NA
validation_cells$all_mhc_pos <- if (length(mhc_present) > 0) rowSums(validation_cells[, mhc_present, drop = FALSE] > 0, na.rm = TRUE) == length(mhc_present) else NA
validation_cells$any_vascular_pos <- if (length(vascular_present) > 0) rowSums(validation_cells[, vascular_present, drop = FALSE] > 0, na.rm = TRUE) > 0 else NA
validation_cells$any_malignant_marker_pos <- if (length(malignant_marker_present) > 0) rowSums(validation_cells[, malignant_marker_present, drop = FALSE] > 0, na.rm = TRUE) > 0 else NA
validation_cells$mhc_cnv_coexpr <- validation_cells$is_cnv_positive & validation_cells$any_mhc_pos
validation_cells$mhc_malignant_marker_coexpr <- validation_cells$any_malignant_marker_pos & validation_cells$any_mhc_pos
validation_cells$vascular_cnv_coexpr <- validation_cells$is_cnv_positive & validation_cells$any_vascular_pos
validation_cells$vascular_malignant_marker_coexpr <- validation_cells$any_malignant_marker_pos & validation_cells$any_vascular_pos

ambient_candidates <- list.files(project_root, pattern = "(SoupX|soupx|CellBender|cellbender|DecontX|decontx|ambient)", recursive = TRUE, full.names = TRUE)
ambient_available <- length(ambient_candidates) > 0
if (ambient_available) {
  log_msg("Ambient-related candidate files detected:", paste(ambient_candidates, collapse = " | "))
} else {
  log_msg("No SoupX / CellBender / DecontX / ambient correction files detected by filename search.")
}

summary_table <- validation_cells |>
  dplyr::group_by(.data$validation_group) |>
  dplyr::summarise(
    n = dplyr::n_distinct(.data$cellID),
    CNV_positive_pct = mean(.data$is_cnv_positive, na.rm = TRUE),
    infercnv_cnv_burden_median = stats::median(.data$infercnv_immune_cnv_burden, na.rm = TRUE),
    infercnv_cnv_correlation_ref_z_median = stats::median(.data$infercnv_immune_cnv_correlation_ref_z, na.rm = TRUE),
    doublet_pct = mean(.data$is_doublet, na.rm = TRUE),
    doubletfinder_pANN_median = stats::median(.data$doubletfinder_pANN, na.rm = TRUE),
    MHC_any_pos_pct = mean(.data$any_mhc_pos, na.rm = TRUE),
    MHC_all_pos_pct = mean(.data$all_mhc_pos, na.rm = TRUE),
    MHC_CNV_coexpr_pct = mean(.data$mhc_cnv_coexpr, na.rm = TRUE),
    MHC_malignant_marker_coexpr_pct = mean(.data$mhc_malignant_marker_coexpr, na.rm = TRUE),
    vascular_any_pos_pct = mean(.data$any_vascular_pos, na.rm = TRUE),
    vascular_CNV_coexpr_pct = mean(.data$vascular_cnv_coexpr, na.rm = TRUE),
    vascular_malignant_marker_coexpr_pct = mean(.data$vascular_malignant_marker_coexpr, na.rm = TRUE),
    ambient_correction_available = ambient_available,
    .groups = "drop"
  )

marker_summary <- validation_cells |>
  dplyr::select(.data$cellID, .data$validation_group, dplyr::all_of(marker_genes)) |>
  summarise_expr("validation_group")

readr::write_csv(validation_cells, file.path("tables", "terminal_marker_validation_per_cell.csv"))
readr::write_csv(summary_table, file.path("tables", "terminal_marker_validation.csv"))
readr::write_csv(marker_summary, file.path("tables", "terminal_marker_validation_marker_summary.csv"))

full_reference_summary <- tibble::tibble()
full_reference_expr <- tibble::tibble()
if (file.exists(full_obj_path)) {
  full_start <- Sys.time()
  full_obj <- readRDS(full_obj_path)
  log_msg("Full object load seconds:", round(as.numeric(difftime(Sys.time(), full_start, units = "secs")), 2))
  full_meta <- full_obj@meta.data |> tibble::rownames_to_column("cellID")
  ref_cells <- full_meta |>
    dplyr::mutate(
      validation_group = dplyr::case_when(
        .data$anno_ident %in% c("Macrophages", "Microglial", "Monocytes", "cDCs", "pDCs") ~ "TAM/myeloid ref",
        .data$anno_ident %in% c("Mural cells") ~ "Mural/pericyte ref",
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(!is.na(.data$validation_group)) |>
    dplyr::select(.data$cellID, .data$validation_group, .data$anno_ident)
  full_expr_mat <- get_data_matrix(full_obj, marker_genes)
  full_expr_df <- as.data.frame(t(as.matrix(full_expr_mat))) |>
    tibble::rownames_to_column("cellID")
  for (g in setdiff(marker_genes, colnames(full_expr_df))) full_expr_df[[g]] <- NA_real_
  full_reference_expr <- ref_cells |>
    dplyr::left_join(full_expr_df, by = "cellID")
  full_reference_summary <- full_reference_expr |>
    dplyr::select(.data$cellID, .data$validation_group, dplyr::all_of(marker_genes)) |>
    summarise_expr("validation_group")
  readr::write_csv(full_reference_summary, file.path("tables", "terminal_marker_validation_reference_marker_summary.csv"))
  readr::write_csv(full_reference_expr, file.path("figures/source_data", "07_reference_marker_expression.csv"))
  rm(full_obj)
  gc()
} else {
  log_msg("Full object not found; reference comparison skipped.")
}

combined_marker_summary <- dplyr::bind_rows(marker_summary, full_reference_summary)
readr::write_csv(combined_marker_summary, file.path("figures/source_data", "07_marker_summary_all_groups.csv"))

cnv_plot_data <- summary_table |>
  dplyr::filter(.data$validation_group %in% c("MES-I terminal", "MES-V terminal", "MES-I all", "MES-V all")) |>
  dplyr::mutate(validation_group = factor(.data$validation_group, levels = c("MES-I terminal", "MES-I all", "MES-V terminal", "MES-V all")))
readr::write_csv(cnv_plot_data, file.path("figures/source_data", "07_CNV比例.csv"))
p_cnv <- ggplot(cnv_plot_data, aes(validation_group, CNV_positive_pct, fill = validation_group)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = paste0(round(CNV_positive_pct * 100, 1), "%\nn=", n)), vjust = -0.15, size = 2.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.08)) +
  scale_fill_manual(values = c("MES-I terminal" = "#BC3C29", "MES-I all" = "#D77A70", "MES-V terminal" = "#20854E", "MES-V all" = "#7DBD87"), guide = "none") +
  labs(title = "CNV-confirmed malignant fraction", x = NULL, y = "CNV+ fraction") +
  theme_bw(base_size = 7) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid.minor = element_blank())
save_plot(p_cnv, "07_CNV比例", width = 4.5, height = 3.2)

doublet_plot_data <- validation_cells |>
  dplyr::filter(.data$validation_group %in% c("MES-I terminal", "MES-V terminal", "MES-I all", "MES-V all")) |>
  dplyr::mutate(validation_group = factor(.data$validation_group, levels = c("MES-I terminal", "MES-I all", "MES-V terminal", "MES-V all")))
readr::write_csv(doublet_plot_data |> dplyr::select(.data$cellID, .data$validation_group, .data$doubletfinder_class, .data$doubletfinder_pANN), file.path("figures/source_data", "07_doublet分布.csv"))
p_doublet <- ggplot(doublet_plot_data, aes(validation_group, doubletfinder_pANN, fill = validation_group)) +
  geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
  geom_boxplot(width = 0.13, outlier.size = 0.2, linewidth = 0.2) +
  scale_fill_manual(values = c("MES-I terminal" = "#BC3C29", "MES-I all" = "#D77A70", "MES-V terminal" = "#20854E", "MES-V all" = "#7DBD87"), guide = "none") +
  labs(title = "DoubletFinder pANN distribution", x = NULL, y = "pANN") +
  theme_bw(base_size = 7) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid.minor = element_blank())
save_plot(p_doublet, "07_doublet分布", width = 4.5, height = 3.2)

coexpr_plot_data <- summary_table |>
  dplyr::filter(.data$validation_group %in% c("MES-I terminal", "MES-V terminal")) |>
  dplyr::select(.data$validation_group, .data$MHC_CNV_coexpr_pct, .data$MHC_malignant_marker_coexpr_pct, .data$vascular_CNV_coexpr_pct, .data$vascular_malignant_marker_coexpr_pct) |>
  tidyr::pivot_longer(-.data$validation_group, names_to = "metric", values_to = "fraction") |>
  dplyr::mutate(
    metric = dplyr::recode(
      .data$metric,
      MHC_CNV_coexpr_pct = "CNV+ & any MHC-II",
      MHC_malignant_marker_coexpr_pct = "Tumor marker+ & any MHC-II",
      vascular_CNV_coexpr_pct = "CNV+ & any vascular marker",
      vascular_malignant_marker_coexpr_pct = "Tumor marker+ & any vascular marker"
    )
  )
readr::write_csv(coexpr_plot_data, file.path("figures/source_data", "07_共表达比例.csv"))
p_coexpr <- ggplot(coexpr_plot_data, aes(metric, fraction, fill = validation_group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.05)) +
  scale_fill_manual(values = c("MES-I terminal" = "#BC3C29", "MES-V terminal" = "#20854E"), name = NULL) +
  labs(title = "Single-cell co-expression sanity", x = NULL, y = "Fraction of cells") +
  theme_bw(base_size = 7) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid.minor = element_blank(), legend.position = "bottom")
save_plot(p_coexpr, "07_共表达比例", width = 5.4, height = 3.3)

if (nrow(full_reference_summary) > 0) {
  terminal_expr_for_plot <- validation_cells |>
    dplyr::filter(.data$validation_group %in% c("MES-I terminal", "MES-V terminal")) |>
    dplyr::select(.data$cellID, .data$validation_group, dplyr::all_of(marker_genes))
  ref_expr_for_plot <- full_reference_expr |>
    dplyr::filter(.data$validation_group %in% c("TAM/myeloid ref", "Mural/pericyte ref")) |>
    dplyr::select(.data$cellID, .data$validation_group, dplyr::all_of(marker_genes))
  expr_plot_all <- dplyr::bind_rows(terminal_expr_for_plot, ref_expr_for_plot)

  mhc_plot_data <- expr_plot_all |>
    dplyr::select(.data$cellID, .data$validation_group, dplyr::all_of(mhc_genes)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(mhc_genes), names_to = "marker", values_to = "expression") |>
    dplyr::filter(.data$validation_group %in% c("MES-I terminal", "TAM/myeloid ref"))
  readr::write_csv(mhc_plot_data, file.path("figures/source_data", "07_MHCII_vs_TAM.csv"))
  p_mhc <- ggplot(mhc_plot_data, aes(validation_group, expression, fill = validation_group)) +
    geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.size = 0.15, linewidth = 0.2) +
    facet_wrap(~marker, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("MES-I terminal" = "#BC3C29", "TAM/myeloid ref" = "#666666"), guide = "none") +
    labs(title = "MHC-II expression: MES-I terminal vs TAM/myeloid reference", x = NULL, y = "Normalized expression") +
    theme_bw(base_size = 7) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid.minor = element_blank())
  save_plot(p_mhc, "07_MHCII_vs_TAM", width = 5.4, height = 3.2)

  vascular_plot_data <- expr_plot_all |>
    dplyr::select(.data$cellID, .data$validation_group, dplyr::all_of(vascular_genes)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(vascular_genes), names_to = "marker", values_to = "expression") |>
    dplyr::filter(.data$validation_group %in% c("MES-V terminal", "Mural/pericyte ref"))
  readr::write_csv(vascular_plot_data, file.path("figures/source_data", "07_血管marker_vs_mural.csv"))
  p_vascular <- ggplot(vascular_plot_data, aes(validation_group, expression, fill = validation_group)) +
    geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.size = 0.15, linewidth = 0.2) +
    facet_wrap(~marker, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("MES-V terminal" = "#20854E", "Mural/pericyte ref" = "#666666"), guide = "none") +
    labs(title = "Vascular marker expression: MES-V terminal vs mural reference", x = NULL, y = "Normalized expression") +
    theme_bw(base_size = 7) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid.minor = element_blank())
  save_plot(p_vascular, "07_血管marker_vs_mural", width = 6.1, height = 3.2)
}

log_msg("Wrote tables/terminal_marker_validation.csv")
log_msg("Wrote tables/terminal_marker_validation_per_cell.csv")
log_msg("Wrote tables/terminal_marker_validation_marker_summary.csv")
if (nrow(full_reference_summary) > 0) log_msg("Wrote tables/terminal_marker_validation_reference_marker_summary.csv")
log_msg("Validation summary:")
capture.output(print(summary_table)) |> paste(collapse = "\n") |> log_msg()
if (nrow(full_reference_summary) > 0) {
  log_msg("Reference expression summary:")
  capture.output(print(full_reference_summary |> dplyr::filter(.data$marker %in% c(mhc_genes, vascular_genes)))) |> paste(collapse = "\n") |> log_msg()
}
writeLines(capture.output(sessionInfo()), file.path("logs", "07_session_info.txt"))
log_msg("Session info: logs/07_session_info.txt")
log_msg("STOP: terminal marker validation complete. Await audit before Phase 4.")
