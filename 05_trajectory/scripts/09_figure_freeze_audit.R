options(stringsAsFactors = FALSE)

if (basename(getwd()) == "06_恶性细胞拟时序") {
  project_root <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(file.path(project_root, "06_恶性细胞拟时序"))

dir.create("tables", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

log_file <- file.path("logs", "09_figure_freeze_audit.log")
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(...))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_msg("Figure freeze audit started.")

fig_files <- list.files("figures", pattern = "\\.(pdf|png)$", recursive = TRUE, full.names = TRUE)
fig_files <- fig_files[!grepl("figures/source_data", fig_files, fixed = TRUE)]
fig_info <- tibble(
  source_path = normalizePath(fig_files, winslash = "/", mustWork = TRUE),
  relative_path = sub(paste0(normalizePath(getwd(), winslash = "/", mustWork = TRUE), "/"), "", normalizePath(fig_files, winslash = "/", mustWork = TRUE), fixed = TRUE),
  file_name = basename(fig_files),
  extension = tools::file_ext(fig_files),
  bytes = file.info(fig_files)$size,
  modified_time = as.character(file.info(fig_files)$mtime)
) |>
  mutate(
    stem = sub("\\.(pdf|png)$", "", file_name, ignore.case = TRUE),
    is_preview = grepl("preview", file_name, ignore.case = TRUE),
    is_combined_or_multipanel = grepl("^(01_|02_|03_|04_|05_|06_|08_driver_genes_preview)", stem),
    is_single_panel = !is_combined_or_multipanel & !is_preview,
    phase = case_when(
      grepl("^01_", stem) ~ "Phase0_root",
      grepl("^02_", stem) ~ "Phase1_slingshot",
      grepl("^03_", stem) ~ "Phase2A_monocle3",
      grepl("^04_", stem) ~ "Phase2B_root_sensitivity",
      grepl("^05_", stem) ~ "Phase2B_cytotrace",
      grepl("^06_", stem) ~ "Phase3_overlay",
      grepl("^07_", stem) ~ "Phase3_validation",
      grepl("^08_", stem) ~ "Phase4_driver",
      TRUE ~ "other"
    )
  )

source_files <- list.files("figures/source_data", pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
source_info <- tibble(
  source_data_path = normalizePath(source_files, winslash = "/", mustWork = TRUE),
  source_relative_path = sub(paste0(normalizePath(getwd(), winslash = "/", mustWork = TRUE), "/"), "", normalizePath(source_files, winslash = "/", mustWork = TRUE), fixed = TRUE),
  source_file = basename(source_files),
  source_bytes = file.info(source_files)$size,
  source_modified_time = as.character(file.info(source_files)$mtime),
  source_stem = sub("\\.csv$", "", basename(source_files), ignore.case = TRUE)
) |>
  mutate(
    source_phase = case_when(
      grepl("^01_", source_stem) ~ "Phase0_root",
      grepl("^02_", source_stem) ~ "Phase1_slingshot",
      grepl("^03_", source_stem) ~ "Phase2A_monocle3",
      grepl("^04_", source_stem) ~ "Phase2B_root_sensitivity",
      grepl("^05_", source_stem) ~ "Phase2B_cytotrace",
      grepl("^06_", source_stem) ~ "Phase3_overlay",
      grepl("^07_", source_stem) ~ "Phase3_validation",
      grepl("^08_", source_stem) ~ "Phase4_driver",
      TRUE ~ "other"
    )
  )

decision <- tribble(
  ~item_id, ~recommended_location, ~panel_role, ~primary_files, ~source_data_hint, ~action_needed, ~risk_note,
  "M_A_root_CytoTRACE2", "main", "A", "Need regenerate single panel from 01_CytoTRACE2_by_subtype.csv", "01_CytoTRACE2_by_subtype.csv; 拟时序_root候选评估.csv", "regenerate_single_panel", "Use root summary only; full Phase0 multipanel goes supplementary.",
  "M_B_slingshot_curves", "main", "B", "Need split from 02_trajectory_primary or regenerate UMAP+curves", "02_panel_A.csv; 02_panel_B.csv; 02_panel_B_branch_points.csv", "regenerate_single_panel", "Must use malignant-only umap_closer4; not old full-cell UMAP.",
  "M_C_lineage_density", "main", "C", "Need split from 02_trajectory_primary Panel E", "02_panel_E.csv", "regenerate_single_panel", "Use lineage-specific assigned cells only.",
  "M_D_topology_robustness", "main", "D", "Need new compact topology panel from Phase2A/2B", "04_curve_coordinates.csv; 04_lineage_meta.csv; monocle3_terminal_detection.csv", "regenerate_or_redesign", "Do not lead with cross-run pseudotime correlation heatmap.",
  "M_E_MES_endpoint", "main", "E", "Need regenerate single endpoint boxplot from 06_panel_D_endpoint_boxplot.csv", "06_panel_D_endpoint_boxplot.csv; terminal_endpoint_comparison.csv; terminal_marker_validation.csv", "regenerate_single_panel", "Main R1 contribution. Include Fig07 validation caveat.",
  "M_F_diffEnd_volcano", "main", "F", "08_MES_diffEnd_volcano.pdf/png", "08_MES_diffEnd_volcano.csv; tradeseq_diffEnd_MESI_vs_MESV.csv", "polish_existing_single_panel", "Do not claim CD74/RGS5 are top tradeSeq genes.",
  "S1_root_full_metrics", "supplement", "S1", "01_root候选评估.pdf", "01_* source csv; 拟时序_root候选统计检验.csv", "polish_or_keep", "Replace p=0 style if still present; independent y-scale caveat.",
  "S2_monocle3_validation", "supplement", "S2", "03_monocle3_vs_slingshot.pdf", "03_* source csv; monocle3_branchlen_diagnostic.csv", "polish_or_keep", "Explain L1 fine-ordering sensitivity.",
  "S3_root_sensitivity", "supplement", "S3", "04_root_sensitivity.pdf", "04_* source csv", "caption_reframe_or_redesign", "Pseudotime correlations are root-relative; topology is the point.",
  "S4_cytotrace_pseudotime", "supplement", "S4", "05_cytotrace2_pseudotime.pdf", "05_* source csv", "keep_if_title_correct", "Title already corrected; restrict monotonic claim to L3.",
  "S5_marker_validation", "supplement", "S5", "07_CNV比例;07_doublet分布;07_共表达比例;07_MHCII_vs_TAM;07_血管marker_vs_mural", "07_* source csv; terminal_marker_validation.csv", "keep_single_panels", "No formal ambient correction; state limitation.",
  "S6_biology_overlay", "supplement", "S6", "Need split from 06_biology_overlay or use source data", "06_panel_A_metaprogram.csv;06_panel_B_neftel.csv;06_panel_C_MP02_MP04.csv", "regenerate_single_panels", "L1/L2 trends supporting only.",
  "S7_driver_heatmap_GO", "supplement", "S7", "08_L3_driver_heatmap;08_driver_enrichment_dotplot", "08_L3_driver_heatmap.csv;08_driver_enrichment_dotplot.csv;driver_enrichment_background.csv", "keep_single_panels", "GO background is 600-gene tradeSeq universe; not genome-wide discovery.",
  "DROP_L1L2_assoc_bar", "delete_or_table_only", "none", "08_L1L2_association_caveat.pdf/png", "08_L1L2_association_caveat.csv", "exclude_from_figures", "Keep table only if needed; do not include as figure."
)

readr::write_csv(fig_info, file.path("tables", "figure_freeze_inventory.csv"))
readr::write_csv(source_info, file.path("tables", "figure_source_data_inventory.csv"))
readr::write_csv(decision, file.path("tables", "figure_freeze_decision_draft.csv"))

log_msg("Figure files:", nrow(fig_info))
log_msg("Source data files:", nrow(source_info))
log_msg("Decision draft rows:", nrow(decision))
log_msg("Wrote tables/figure_freeze_inventory.csv")
log_msg("Wrote tables/figure_source_data_inventory.csv")
log_msg("Wrote tables/figure_freeze_decision_draft.csv")
log_msg("STOP: audit tables generated; no files archived or deleted.")
