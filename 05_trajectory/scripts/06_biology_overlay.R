suppressPackageStartupMessages({
  .libPaths(c("<R_LIBS>", "<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(qs2)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# Phase 3: biology overlay along trajectory.
# Read-only: no trajectory rerun and no qs2 writeback.

started_at <- Sys.time()
set.seed(42)

cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(cwd) == "scripts") {
  step_dir <- normalizePath(file.path(cwd, ".."), winslash = "/", mustWork = TRUE)
} else if (file.exists(file.path(cwd, "scripts", "_naming.R"))) {
  step_dir <- cwd
} else {
  step_dir <- normalizePath(file.path(cwd, "06_恶性细胞拟时序"), winslash = "/", mustWork = TRUE)
}
setwd(step_dir)

dir.create("tables", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/source_data", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

source(file.path("scripts", "_naming.R"))

log_file <- file.path("logs", "06_biology_overlay.log")
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_msg("Phase 3 biology overlay started")
log_msg("set.seed = 42")

obj_path <- normalizePath(file.path("..", "05_恶性细胞分亚群与Neftel对照", "outputs", "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"), winslash = "/", mustWork = TRUE)
slingshot_path <- file.path("tables", "slingshot_pseudotime_per_cell.csv")
lineage_weight_path <- file.path("tables", "slingshot_lineage_weight_per_cell.csv")
cytotrace_path <- file.path("outputs", "CytoTRACE2_per_cell.csv")
if (!file.exists(slingshot_path)) stop("Missing Slingshot pseudotime table.")
if (!file.exists(lineage_weight_path)) stop("Missing lineage weight table.")
if (!file.exists(cytotrace_path)) stop("Missing cached CytoTRACE2 table.")

log_msg("Input object:", obj_path)
log_msg("Slingshot table:", slingshot_path)
log_msg("Lineage weight table:", lineage_weight_path)
log_msg("Cached CytoTRACE2:", cytotrace_path)
log_msg("MP02 program: ECM / vascular-niche MES support")
log_msg("MP04 program: MHC-II antigen-presentation / immune-interacting MES support")

load_start <- Sys.time()
obj <- qs2::qs_read(obj_path)
log_msg("Object load seconds:", round(as.numeric(difftime(Sys.time(), load_start, units = "secs")), 2))

mp_cols <- paste0("MP0", 1:6, "_UCell")
neftel_cols <- c("AMS_NPC1", "AMS_NPC2", "AMS_OPC", "AMS_MES1", "AMS_MES2", "AMS_AC")
required_cols <- c("subtype_k4", mp_cols, neftel_cols)
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))

md <- obj@meta.data |>
  rownames_to_column("cellID") |>
  recode_subtype() |>
  dplyr::select(cellID, subtype_k4, subtype_label_final, subtype_short, dplyr::all_of(mp_cols), dplyr::all_of(neftel_cols))
cytotrace <- readr::read_csv(cytotrace_path, show_col_types = FALSE)
md <- md |>
  dplyr::left_join(cytotrace, by = c("cellID" = "cell_id"))

slingshot <- readr::read_csv(slingshot_path, show_col_types = FALSE)
trajectory_cells <- slingshot |>
  dplyr::select(cellID, subtype_k4, subtype_label_final, subtype_short, lineage_assignment, dplyr::starts_with("lineage")) |>
  dplyr::left_join(md, by = c("cellID", "subtype_k4", "subtype_label_final", "subtype_short"))

lineage_ids <- paste0("lineage", 1:3)
assigned_long <- lapply(lineage_ids, function(lin) {
  pst_col <- paste0(lin, "_pst")
  weight_col <- paste0(lin, "_weight")
  trajectory_cells |>
    dplyr::filter(lineage_assignment == lin, .data[[weight_col]] > 0.5, !is.na(.data[[pst_col]])) |>
    dplyr::mutate(
      lineage_id = lin,
      lineage_label = factor(lin, levels = lineage_ids, labels = paste0("L", 1:3)),
      pseudotime = .data[[pst_col]],
      lineage_weight = .data[[weight_col]]
    )
}) |>
  bind_rows()
if (nrow(assigned_long) != 28213) {
  log_msg("WARNING assigned cell total != 28213:", nrow(assigned_long))
}

bin_scores <- function(df, score_cols, score_type) {
  df |>
    dplyr::group_by(lineage_id) |>
    dplyr::mutate(pseudotime_bin = dplyr::ntile(pseudotime, 20)) |>
    dplyr::ungroup() |>
    dplyr::select(cellID, subtype_k4, subtype_label_final, subtype_short, lineage_id, lineage_label, pseudotime, pseudotime_bin, dplyr::all_of(score_cols)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(score_cols), names_to = "score_name", values_to = "score_value") |>
    dplyr::group_by(lineage_id, lineage_label, pseudotime_bin, score_name) |>
    dplyr::summarise(
      n_cells = dplyr::n(),
      pseudotime_mean = mean(pseudotime, na.rm = TRUE),
      score_mean = mean(score_value, na.rm = TRUE),
      score_median = median(score_value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(score_type = score_type)
}

mp_binned <- bin_scores(assigned_long, mp_cols, "NMF metaprogram")
neftel_binned <- bin_scores(assigned_long, neftel_cols, "Neftel AMS")
readr::write_csv(mp_binned, file.path("tables", "metaprogram_along_pseudotime.csv"))
readr::write_csv(neftel_binned, file.path("tables", "neftel_along_pseudotime.csv"))

marker_genes <- c("RGS5", "ACTA2", "TAGLN", "HIGD1B", "HLA-DRA", "HLA-DPB1", "CD74", "B2M")
available_markers <- intersect(marker_genes, rownames(obj))
missing_markers <- setdiff(marker_genes, available_markers)
if (length(missing_markers) > 0) log_msg("WARNING missing marker genes:", paste(missing_markers, collapse = ", "))
rna_data <- tryCatch(
  Seurat::GetAssayData(obj, assay = "RNA", layer = "data"),
  error = function(e) Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
)
expr <- Matrix::t(rna_data[available_markers, assigned_long$cellID, drop = FALSE]) |>
  as.matrix() |>
  as.data.frame() |>
  rownames_to_column("cellID")

endpoint_cells_raw <- assigned_long |>
  dplyr::filter(lineage_id %in% c("lineage1", "lineage2")) |>
  dplyr::group_by(lineage_id) |>
  dplyr::mutate(endpoint_cutoff = stats::quantile(pseudotime, 0.80, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::filter(pseudotime >= endpoint_cutoff) |>
  dplyr::mutate(
    endpoint_group = dplyr::case_when(
      lineage_id == "lineage1" ~ "MES-I terminal",
      lineage_id == "lineage2" ~ "MES-V terminal",
      TRUE ~ NA_character_
    ),
    endpoint_group = factor(endpoint_group, levels = c("MES-I terminal", "MES-V terminal"))
  ) |>
  dplyr::left_join(expr, by = "cellID")

endpoint_composition_raw <- endpoint_cells_raw |>
  dplyr::count(endpoint_group, lineage_id, subtype_short, name = "n_cells") |>
  dplyr::group_by(endpoint_group) |>
  dplyr::mutate(endpoint_fraction = n_cells / sum(n_cells)) |>
  dplyr::ungroup()
readr::write_csv(endpoint_composition_raw, file.path("tables", "terminal_endpoint_composition_pseudotime_top20_raw.csv"))

endpoint_cells <- assigned_long |>
  dplyr::filter(
    (lineage_id == "lineage1" & subtype_short == "MES-I") |
      (lineage_id == "lineage2" & subtype_short == "MES-V")
  ) |>
  dplyr::group_by(lineage_id, subtype_short) |>
  dplyr::mutate(endpoint_cutoff = stats::quantile(pseudotime, 0.80, na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::filter(pseudotime >= endpoint_cutoff) |>
  dplyr::mutate(
    endpoint_group = dplyr::case_when(
      lineage_id == "lineage1" & subtype_short == "MES-I" ~ "MES-I terminal",
      lineage_id == "lineage2" & subtype_short == "MES-V" ~ "MES-V terminal",
      TRUE ~ NA_character_
    ),
    endpoint_group = factor(endpoint_group, levels = c("MES-I terminal", "MES-V terminal"))
  ) |>
  dplyr::left_join(expr, by = "cellID")

endpoint_composition <- endpoint_cells |>
  dplyr::count(endpoint_group, lineage_id, subtype_short, name = "n_cells") |>
  dplyr::group_by(endpoint_group) |>
  dplyr::mutate(endpoint_fraction = n_cells / sum(n_cells)) |>
  dplyr::ungroup()
readr::write_csv(endpoint_composition, file.path("tables", "terminal_endpoint_composition.csv"))

endpoint_feature_df <- endpoint_cells |>
  dplyr::select(cellID, endpoint_group, subtype_k4, subtype_label_final, subtype_short, lineage_id, pseudotime, MP02_UCell, MP04_UCell, dplyr::all_of(available_markers)) |>
  tidyr::pivot_longer(cols = c("MP02_UCell", "MP04_UCell", dplyr::all_of(available_markers)), names_to = "feature", values_to = "value") |>
  dplyr::mutate(
    feature_class = dplyr::case_when(
      feature %in% c("MP02_UCell", "MP04_UCell") ~ "NMF metaprogram",
      feature %in% c("RGS5", "ACTA2", "TAGLN", "HIGD1B") ~ "Vascular-niche marker",
      TRUE ~ "Immune-interacting marker"
    ),
    feature_label = dplyr::recode(feature, MP02_UCell = "MP02", MP04_UCell = "MP04")
  )

endpoint_stats <- endpoint_feature_df |>
  dplyr::group_by(feature, feature_label, feature_class, endpoint_group) |>
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    pct_pos = mean(value > 0, na.rm = TRUE),
    .groups = "drop"
  )
endpoint_tests <- endpoint_feature_df |>
  dplyr::group_by(feature, feature_label, feature_class) |>
  dplyr::summarise(
    p_value = suppressWarnings(wilcox.test(value ~ endpoint_group, exact = FALSE)$p.value),
    mean_MESI = mean(value[endpoint_group == "MES-I terminal"], na.rm = TRUE),
    mean_MESV = mean(value[endpoint_group == "MES-V terminal"], na.rm = TRUE),
    delta_MESV_minus_MESI = mean_MESV - mean_MESI,
    .groups = "drop"
  ) |>
  dplyr::mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
endpoint_comparison <- endpoint_stats |>
  dplyr::left_join(endpoint_tests, by = c("feature", "feature_label", "feature_class")) |>
  dplyr::arrange(feature_class, feature)
readr::write_csv(endpoint_comparison, file.path("tables", "terminal_endpoint_comparison.csv"))

mp_labels <- c(
  MP01_UCell = "MP01",
  MP02_UCell = "MP02",
  MP03_UCell = "MP03",
  MP04_UCell = "MP04",
  MP05_UCell = "MP05",
  MP06_UCell = "MP06"
)
neftel_labels <- c(
  AMS_NPC1 = "NPC1",
  AMS_NPC2 = "NPC2",
  AMS_OPC = "OPC",
  AMS_MES1 = "MES1",
  AMS_MES2 = "MES2",
  AMS_AC = "AC"
)
mp_colors <- c(MP01_UCell = "#4E79A7", MP02_UCell = "#20854E", MP03_UCell = "#59A14F", MP04_UCell = "#BC3C29", MP05_UCell = "#8C6D31", MP06_UCell = "#9467BD")
neftel_colors <- c(AMS_NPC1 = "#0072B5", AMS_NPC2 = "#56B4E9", AMS_OPC = "#E18727", AMS_MES1 = "#BC3C29", AMS_MES2 = "#20854E", AMS_AC = "#7F7F7F")

mp_plot <- mp_binned |>
  dplyr::mutate(score_label = factor(mp_labels[score_name], levels = unname(mp_labels)))
neftel_plot <- neftel_binned |>
  dplyr::mutate(score_label = factor(neftel_labels[score_name], levels = unname(neftel_labels)))

p_a <- ggplot(mp_plot, aes(pseudotime_mean, score_mean, color = score_name, linetype = lineage_id)) +
  geom_line(linewidth = 0.45) +
  facet_wrap(~ lineage_label, scales = "free_x", nrow = 1) +
  scale_color_manual(values = mp_colors, labels = mp_labels, name = "MP") +
  scale_linetype_manual(values = c(lineage1 = "22", lineage2 = "22", lineage3 = "solid"), guide = "none") +
  labs(title = "A. NMF metaprograms along pseudotime", subtitle = "L1/L2 fine-ordering is uncertain; L3 trend is more reliable", x = "Pseudotime", y = "Binned mean score") +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title = element_text(size = 8), plot.subtitle = element_text(size = 6.4), strip.text = element_text(size = 7))

p_b <- ggplot(neftel_plot, aes(pseudotime_mean, score_mean, color = score_name, linetype = lineage_id)) +
  geom_line(linewidth = 0.45) +
  facet_wrap(~ lineage_label, scales = "free_x", nrow = 1) +
  scale_color_manual(values = neftel_colors, labels = neftel_labels, name = "Neftel") +
  scale_linetype_manual(values = c(lineage1 = "22", lineage2 = "22", lineage3 = "solid"), guide = "none") +
  labs(title = "B. Neftel AMS states along pseudotime", x = "Pseudotime", y = "Binned mean AMS") +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title = element_text(size = 8), strip.text = element_text(size = 7))

panel_c <- mp_plot |>
  dplyr::filter(lineage_id %in% c("lineage1", "lineage2"), score_name %in% c("MP02_UCell", "MP04_UCell")) |>
  dplyr::mutate(curve = paste0(as.character(lineage_label), "-", mp_labels[score_name]))
curve_colors <- c("L1-MP02" = "#7DBD87", "L1-MP04" = "#BC3C29", "L2-MP02" = "#20854E", "L2-MP04" = "#D77A70")
p_c <- ggplot(panel_c, aes(pseudotime_mean, score_mean, color = curve)) +
  geom_line(linewidth = 0.65) +
  scale_color_manual(values = curve_colors, name = NULL) +
  labs(title = "C. MP02/MP04 divergence in MES-like lineages", subtitle = "Endpoint interpretation is prioritized over fine-scale ordering", x = "Pseudotime", y = "Binned mean score") +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", plot.title = element_text(size = 8), plot.subtitle = element_text(size = 6.4))

plot_features <- c("MP02_UCell", "MP04_UCell", "RGS5", "ACTA2", "TAGLN", "HLA-DRA", "HLA-DPB1", "CD74")
panel_d <- endpoint_feature_df |>
  dplyr::filter(feature %in% plot_features) |>
  dplyr::mutate(feature_label = factor(feature_label, levels = c("MP02", "MP04", setdiff(plot_features, c("MP02_UCell", "MP04_UCell")))))
test_labels <- endpoint_tests |>
  dplyr::filter(feature %in% plot_features) |>
  dplyr::mutate(
    feature_label = dplyr::recode(feature, MP02_UCell = "MP02", MP04_UCell = "MP04"),
    feature_label = factor(feature_label, levels = levels(panel_d$feature_label)),
    label = paste0("BH p=", format.pval(p_adj_BH, digits = 2, eps = 1e-300))
  )
p_d <- ggplot(panel_d, aes(endpoint_group, value, fill = endpoint_group)) +
  geom_boxplot(outlier.size = 0.12, linewidth = 0.25) +
  facet_wrap(~ feature_label, scales = "free_y", nrow = 2) +
  geom_text(data = test_labels, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.2, size = 1.8) +
  scale_fill_manual(values = c("MES-I terminal" = "#BC3C29", "MES-V terminal" = "#20854E"), guide = "none") +
  labs(title = "D. MES-I vs MES-V endpoint comparison", subtitle = "Strict endpoint subtype cells: MES-I in L1 and MES-V in L2, top 20% within each subtype-lineage", x = NULL, y = "Score / normalized expression") +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1, size = 5.8), plot.title = element_text(size = 8), plot.subtitle = element_text(size = 6.4), strip.text = element_text(size = 6.5))

readr::write_csv(mp_plot, file.path("figures/source_data", "06_panel_A_metaprogram.csv"))
readr::write_csv(neftel_plot, file.path("figures/source_data", "06_panel_B_neftel.csv"))
readr::write_csv(panel_c, file.path("figures/source_data", "06_panel_C_MP02_MP04.csv"))
readr::write_csv(panel_d, file.path("figures/source_data", "06_panel_D_endpoint_boxplot.csv"))

p_all <- (p_a / p_b / (p_c | p_d)) +
  patchwork::plot_layout(heights = c(1, 1, 1.25)) +
  patchwork::plot_annotation(
    caption = "Continuous trends for L1/L2 are supporting only because MES-like fine ordering is method-sensitive. Endpoint comparisons provide the primary evidence for MES-V/MES-I program divergence.",
    theme = theme(plot.caption = element_text(size = 6, hjust = 0))
  )

pdf(file.path("figures", "06_biology_overlay.pdf"), width = 11.0, height = 11.5, useDingbats = FALSE)
print(p_all)
dev.off()
ggsave(file.path("figures", "06_biology_overlay_preview.png"), p_all, width = 11.0, height = 11.5, dpi = 180, bg = "white")

log_msg("Wrote tables/metaprogram_along_pseudotime.csv")
log_msg("Wrote tables/neftel_along_pseudotime.csv")
log_msg("Wrote tables/terminal_endpoint_comparison.csv")
log_msg("Wrote figures/06_biology_overlay.pdf")
log_msg("Endpoint composition:")
capture.output(print(endpoint_composition)) |> paste(collapse = "\n") |> log_msg()
log_msg("Endpoint tests:")
capture.output(print(endpoint_tests)) |> paste(collapse = "\n") |> log_msg()
session_info <- capture.output(sessionInfo())
cat(session_info, sep = "\n", file = file.path("logs", "06_session_info.txt"))
log_msg("Session info:", file.path("logs", "06_session_info.txt"))
log_msg("Total elapsed seconds:", round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 2))
log_msg("STOP: Phase 3 complete. Await audit before Phase 4.")
