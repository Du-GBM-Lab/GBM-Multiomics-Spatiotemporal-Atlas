options(stringsAsFactors = FALSE)
set.seed(42)

if (basename(getwd()) == "06_恶性细胞拟时序") {
  project_root <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(file.path(project_root, "06_恶性细胞拟时序"))

dir.create("figures/final_panels/main", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/final_panels/source_data/main", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggridges)
  library(ggrepel)
})

log_file <- "logs/10_make_main_single_panels.log"
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(...))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_msg("Main single-panel figure generation started.")

pal <- c("NPC-P" = "#0072B5", "OPC-M" = "#E18727", "MES-V" = "#20854E", "MES-I" = "#BC3C29")
lineage_labs <- c("lineage1" = "L1: NPC-P to MES-I", "lineage2" = "L2: NPC-P to MES-V", "lineage3" = "L3: NPC-P to OPC-M")
lineage_pal <- c("lineage1" = "#8B1A1A", "lineage2" = "#006D2C", "lineage3" = "#7A4A00")

theme_main <- function(base_size = 7) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.18, color = "grey90"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = base_size + 1, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.5),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.4),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.4),
      strip.text = element_text(size = base_size, face = "bold")
    )
}

theme_umap <- function(base_size = 7) {
  theme_void(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = base_size + 1, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.5),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.4)
    )
}

save_panel <- function(plot, stem, width, height) {
  pdf_path <- file.path("figures/final_panels/main", paste0(stem, ".pdf"))
  png_path <- file.path("figures/final_panels/main", paste0(stem, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(png_path, plot, width = width, height = height, dpi = 400, bg = "white")
  log_msg("Wrote", pdf_path)
  log_msg("Wrote", png_path)
}

copy_source <- function(src, dest_name) {
  dest <- file.path("figures/final_panels/source_data/main", dest_name)
  file.copy(src, dest, overwrite = TRUE)
  log_msg("Copied source data", dest)
}

# Panel A: root summary / CytoTRACE2
cyto <- readr::read_csv("figures/source_data/01_CytoTRACE2_by_subtype.csv", show_col_types = FALSE) |>
  mutate(subtype_short = factor(subtype_short, levels = names(pal)))
root_metrics <- readr::read_csv("tables/拟时序_root候选评估.csv", show_col_types = FALSE)
root_summary <- root_metrics |>
  filter(metric_name %in% c("CytoTRACE2_score", "SOX2", "SOX11", "AMS_NPC1", "AMS_NPC2")) |>
  select(subtype_short, metric_name, mean, median)
readr::write_csv(cyto, "figures/final_panels/source_data/main/Main_A_CytoTRACE2_by_subtype.csv")
readr::write_csv(root_summary, "figures/final_panels/source_data/main/Main_A_root_metric_summary.csv")
p_a <- ggplot(cyto, aes(subtype_short, CytoTRACE2_score, fill = subtype_short)) +
  geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.08, linewidth = 0.2) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(
    title = "A. NPC-P shows the highest stemness score",
    subtitle = "Root choice is supported by CytoTRACE2, cell cycle, SOX family and Neftel NPC scores",
    x = NULL,
    y = "CytoTRACE2 score"
  ) +
  theme_main(7) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_panel(p_a, "Main_A_root_CytoTRACE2", 3.6, 3.1)

# Panel B: UMAP + Slingshot curves
panel_a_cells <- readr::read_csv("figures/source_data/02_panel_A.csv", show_col_types = FALSE) |>
  mutate(subtype_short = factor(subtype_short, levels = names(pal)))
panel_b <- readr::read_csv("figures/source_data/02_panel_B.csv", show_col_types = FALSE)
curve <- panel_b |>
  filter(layer %in% c("curve", "slingshot_curve")) |>
  mutate(lineage_id = factor(lineage_id, levels = names(lineage_labs)))
lineage_terminal <- c("lineage1" = "MES-I", "lineage2" = "MES-V", "lineage3" = "OPC-M")
curve_start <- curve |>
  mutate(point_id_num = suppressWarnings(as.numeric(point_id))) |>
  group_by(lineage_id) |>
  arrange(point_id_num, .by_group = TRUE) |>
  slice_head(n = 1) |>
  ungroup()
curve_start_one <- curve_start |>
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    label = "Start: NPC-P"
  )
curve_endpoints <- curve |>
  mutate(point_id_num = suppressWarnings(as.numeric(point_id))) |>
  group_by(lineage_id) |>
  arrange(point_id_num, .by_group = TRUE) |>
  slice_tail(n = 1) |>
  ungroup() |>
  mutate(
    lineage_short = recode(as.character(lineage_id), lineage1 = "L1", lineage2 = "L2", lineage3 = "L3"),
    end_cluster = lineage_terminal[as.character(lineage_id)],
    endpoint_label = paste0(lineage_short, " End: ", end_cluster)
  )
root_centroid <- panel_a_cells |>
  filter(subtype_short == "NPC-P") |>
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
readr::write_csv(panel_a_cells, "figures/final_panels/source_data/main/Main_B_umap_cells.csv")
readr::write_csv(curve, "figures/final_panels/source_data/main/Main_B_slingshot_curves.csv")
readr::write_csv(curve_start_one, "figures/final_panels/source_data/main/Main_B_curve_start.csv")
readr::write_csv(curve_endpoints, "figures/final_panels/source_data/main/Main_B_curve_endpoints.csv")
p_b <- ggplot() +
  geom_point(data = panel_a_cells, aes(UMAP_1, UMAP_2, color = subtype_short), size = 0.08, alpha = 0.45) +
  geom_path(data = curve, aes(UMAP_1, UMAP_2, group = lineage_id, linewidth = lineage_id), color = "white", alpha = 0.85, lineend = "round") +
  geom_path(
    data = curve,
    aes(UMAP_1, UMAP_2, group = lineage_id, color = lineage_id),
    linewidth = 0.72,
    lineend = "round",
    arrow = arrow(type = "closed", length = unit(0.065, "inches"))
  ) +
  geom_point(data = curve_start_one, aes(UMAP_1, UMAP_2), shape = 21, fill = "white", color = "black", size = 2.1, stroke = 0.5) +
  geom_point(data = curve_endpoints, aes(UMAP_1, UMAP_2, fill = end_cluster), shape = 24, color = "black", size = 1.9, stroke = 0.35) +
  geom_label(
    data = curve_start_one,
    aes(UMAP_1, UMAP_2 - 1.25, label = label),
    size = 1.9,
    linewidth = 0.15,
    label.padding = unit(0.10, "lines"),
    fill = "white",
    color = "black"
  ) +
  ggrepel::geom_label_repel(
    data = curve_endpoints,
    aes(UMAP_1, UMAP_2, label = endpoint_label),
    size = 2.1,
    linewidth = 0.15,
    label.padding = unit(0.10, "lines"),
    min.segment.length = 0,
    segment.size = 0.15,
    fill = "white",
    color = "black",
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c(pal, lineage_pal),
    breaks = names(pal),
    name = NULL
  ) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_linewidth_manual(values = c("lineage1" = 1.6, "lineage2" = 1.6, "lineage3" = 1.6), guide = "none") +
  labs(title = "B. NPC-P-rooted three-lineage trajectory", subtitle = "Malignant-only UMAP (umap_closer4)") +
  coord_equal() +
  theme_umap(7) +
  theme(legend.position = "right")
save_panel(p_b, "Main_B_slingshot_curves", 4.2, 3.6)

# Panel C: lineage-specific pseudotime density
density_lineage_labs <- c(
  "lineage1" = "L1: NPC-P to MES-I (n = 11,660)",
  "lineage2" = "L2: NPC-P to MES-V (n = 9,232)",
  "lineage3" = "L3: NPC-P to OPC-M (n = 7,321)"
)
density_df <- readr::read_csv("figures/source_data/02_panel_E.csv", show_col_types = FALSE) |>
  mutate(
    subtype_short = factor(subtype_short, levels = names(pal)),
    lineage_label = factor(density_lineage_labs[lineage_id], levels = density_lineage_labs)
  )
density_counts <- density_df |>
  count(lineage_id, lineage_label, subtype_short, name = "n_cells")
density_plot <- density_df |>
  left_join(density_counts, by = c("lineage_id", "lineage_label", "subtype_short")) |>
  filter(n_cells >= 50)
readr::write_csv(density_df, "figures/final_panels/source_data/main/Main_C_lineage_pseudotime_density.csv")
readr::write_csv(density_counts, "figures/final_panels/source_data/main/Main_C_lineage_density_counts.csv")
readr::write_csv(density_plot, "figures/final_panels/source_data/main/Main_C_lineage_pseudotime_density_plotted.csv")
p_c <- ggplot(density_plot, aes(pseudotime, subtype_short, fill = subtype_short)) +
  geom_density_ridges(alpha = 0.78, scale = 1.1, rel_min_height = 0.01, linewidth = 0.25, color = "white") +
  facet_wrap(~lineage_label, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = pal, guide = "none") +
  labs(
    title = "C. Pseudotime density by assigned lineage",
    subtitle = "NPC-P-to-terminal ordering; subtype-lineage combinations with <50 cells are omitted",
    x = "Pseudotime",
    y = NULL
  ) +
  theme_main(7) +
  theme(
    panel.grid.major.y = element_blank(),
    strip.text = element_text(size = 6.7, face = "bold")
  )
save_panel(p_c, "Main_C_lineage_density", 4.2, 5.0)

# Panel D: topology robustness summary
curve_multi <- readr::read_csv("figures/source_data/04_curve_coordinates.csv", show_col_types = FALSE)
lineage_meta <- readr::read_csv("figures/source_data/04_lineage_meta.csv", show_col_types = FALSE)
run_labels <- c(
  "run1_NPCP" = "Root: NPC-P",
  "run2_OPCM" = "Root: OPC-M",
  "run3_MESV" = "Root: MES-V",
  "run4_auto" = "Algorithm-selected root"
)
curve_multi <- curve_multi |>
  mutate(run_label = factor(run_labels[run_id], levels = run_labels)) |>
  left_join(
    lineage_meta |> select(run_id, lineage_id, start_cluster, end_cluster, lineage_path),
    by = c("run_id", "lineage_id")
  ) |>
  group_by(run_id, lineage_id) |>
  arrange(point_id, .by_group = TRUE) |>
  mutate(point_order = row_number(), n_points = n()) |>
  ungroup()
curve_start <- curve_multi |>
  filter(point_order == 1) |>
  mutate(label = paste0("Root: ", start_cluster))
curve_start_one <- curve_start |>
  group_by(run_id, run_label, start_cluster) |>
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    label = paste0("Root: ", unique(start_cluster)[1]),
    .groups = "drop"
  )
curve_end <- curve_multi |>
  filter(point_order == n_points) |>
  mutate(label = paste0("End: ", end_cluster))
lineage_path_pal <- c(
  "MES-I" = "#BC3C29",
  "MES-V" = "#20854E",
  "OPC-M" = "#E18727"
)
hub_df <- tibble(
  run_label = factor(run_labels, levels = run_labels),
  UMAP_1 = root_centroid$UMAP_1,
  UMAP_2 = root_centroid$UMAP_2
)
umap_down <- panel_a_cells |>
  group_by(subtype_short) |>
  group_modify(~ dplyr::slice_sample(.x, n = min(1800, nrow(.x)))) |>
  ungroup()
readr::write_csv(curve_multi, "figures/final_panels/source_data/main/Main_D_multiroot_curves.csv")
readr::write_csv(lineage_meta, "figures/final_panels/source_data/main/Main_D_multiroot_lineage_meta.csv")
p_d <- ggplot() +
  geom_point(data = umap_down, aes(UMAP_1, UMAP_2, color = subtype_short), size = 0.07, alpha = 0.20) +
  geom_path(
    data = curve_multi,
    aes(UMAP_1, UMAP_2, group = interaction(run_id, lineage_id)),
    color = "white",
    linewidth = 0.95,
    alpha = 0.85
  ) +
  geom_path(
    data = curve_multi,
    aes(UMAP_1, UMAP_2, group = interaction(run_id, lineage_id), color = end_cluster),
    linewidth = 0.55,
    alpha = 0.92,
    arrow = arrow(type = "closed", length = unit(0.06, "inches"))
  ) +
  geom_point(data = curve_start_one, aes(UMAP_1, UMAP_2), shape = 21, fill = "white", color = "black", size = 1.45, stroke = 0.35) +
  geom_point(data = curve_end, aes(UMAP_1, UMAP_2, fill = end_cluster), shape = 24, color = "black", size = 1.45, stroke = 0.25) +
  ggrepel::geom_label_repel(
    data = curve_start_one,
    aes(UMAP_1, UMAP_2, label = label),
    size = 1.55,
    label.size = 0.10,
    label.padding = unit(0.06, "lines"),
    fill = "white",
    color = "black",
    min.segment.length = 0,
    segment.size = 0.10,
    max.overlaps = Inf
  ) +
  ggrepel::geom_label_repel(
    data = curve_end,
    aes(UMAP_1, UMAP_2, label = label),
    size = 1.45,
    label.size = 0.10,
    label.padding = unit(0.06, "lines"),
    fill = "white",
    color = "black",
    min.segment.length = 0,
    segment.size = 0.10,
    max.overlaps = Inf
  ) +
  geom_point(data = hub_df, aes(UMAP_1, UMAP_2), shape = 21, fill = "white", color = "black", size = 1.4, stroke = 0.35) +
  geom_label(
    data = hub_df,
    aes(UMAP_1, UMAP_2 - 1.0, label = "NPC-P hub"),
    size = 1.65,
    label.size = 0.12,
    label.padding = unit(0.08, "lines"),
    fill = "white",
    color = "black"
  ) +
  facet_wrap(~run_label, ncol = 2) +
  scale_color_manual(values = c(pal, lineage_path_pal), breaks = names(pal), name = NULL) +
  scale_fill_manual(values = lineage_path_pal, guide = "none") +
  labs(title = "D. Root sensitivity supports NPC-P as central topology hub", subtitle = "Pseudotime is root-relative; topology is compared here") +
  coord_equal() +
  theme_umap(7) +
  theme(legend.position = "bottom", strip.text = element_text(size = 6.5, face = "bold"))
save_panel(p_d, "Main_D_topology_robustness", 5.2, 4.4)

# Panel E: MES endpoint comparison
endpoint_df <- readr::read_csv("figures/source_data/06_panel_D_endpoint_boxplot.csv", show_col_types = FALSE)
plot_features <- c("MP02_UCell", "MP04_UCell", "RGS5", "ACTA2", "TAGLN", "HLA-DRA", "HLA-DPB1", "CD74")
endpoint_plot <- endpoint_df |>
  filter(feature %in% plot_features) |>
  mutate(
    feature_label = recode(feature, MP02_UCell = "MP02", MP04_UCell = "MP04"),
    feature_label = factor(feature_label, levels = c("MP02", "MP04", "RGS5", "ACTA2", "TAGLN", "HLA-DRA", "HLA-DPB1", "CD74")),
    endpoint_group = factor(endpoint_group, levels = c("MES-I terminal", "MES-V terminal")),
    endpoint_short = factor(
      ifelse(endpoint_group == "MES-I terminal", "MES-I", "MES-V"),
      levels = c("MES-I", "MES-V")
    )
  )
readr::write_csv(endpoint_plot, "figures/final_panels/source_data/main/Main_E_MES_endpoint_boxplot.csv")
p_e <- ggplot(endpoint_plot, aes(endpoint_short, value, fill = endpoint_short)) +
  geom_boxplot(outlier.size = 0.12, linewidth = 0.25, width = 0.55) +
  facet_wrap(~feature_label, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = c("MES-I" = "#BC3C29", "MES-V" = "#20854E"), guide = "none") +
  labs(
    title = "E. MES-like termini diverge in vascular and immune-reactive programs",
    subtitle = "Strict endpoints; marker genuineness validated by CNV/doublet/reference checks",
    x = NULL,
    y = "Score / normalized expression"
  ) +
  theme_main(7) +
  theme(axis.text.x = element_text(size = 6.2), strip.text = element_text(size = 6.4))
save_panel(p_e, "Main_E_MES_endpoint_comparison", 6.2, 4.2)

# Panel F: diffEnd volcano
volcano <- readr::read_csv("figures/source_data/08_MES_diffEnd_volcano.csv", show_col_types = FALSE) |>
  mutate(
    direction = recode(direction, no_direction = "Not directional"),
    direction = factor(direction, levels = c("MES-I high", "MES-V high", "Not directional"))
  )
label_genes <- c("HLA-DRA", "HLA-DPB1", "HLA-DPA1", "ACTA2", "TAGLN", "HIGD1B", "PLAU", "MFAP4", "SOX2", "PTPRZ1")
volcano <- volcano |> mutate(label_gene_main = gene %in% label_genes)
readr::write_csv(volcano, "figures/final_panels/source_data/main/Main_F_diffEnd_volcano.csv")
p_f <- ggplot(volcano, aes(log2FC_MESV_vs_MESI, neglog10_BH_capped)) +
  geom_point(aes(color = direction), alpha = 0.65, size = 0.85) +
  ggrepel::geom_text_repel(
    data = filter(volcano, label_gene_main),
    aes(label = gene),
    size = 2.1,
    max.overlaps = 30,
    box.padding = 0.25,
    min.segment.length = 0,
    segment.size = 0.15
  ) +
  scale_color_manual(values = c("MES-I high" = "#BC3C29", "MES-V high" = "#20854E", "Not directional" = "grey70"), name = NULL) +
  labs(
    title = "F. Candidate switch genes distinguish MES-like termini",
    subtitle = "Focused tradeSeq screen; GO uses 600-gene testing universe",
    x = "Endpoint log2FC (MES-V / MES-I)",
    y = "-log10(BH), capped at 60"
  ) +
  theme_main(7) +
  theme(legend.position = "bottom")
save_panel(p_f, "Main_F_MES_diffEnd_volcano", 4.8, 3.8)

manifest <- tibble(
  panel = LETTERS[1:6],
  file_stem = c(
    "Main_A_root_CytoTRACE2",
    "Main_B_slingshot_curves",
    "Main_C_lineage_density",
    "Main_D_topology_robustness",
    "Main_E_MES_endpoint_comparison",
    "Main_F_MES_diffEnd_volcano"
  ),
  pdf = file.path("figures/final_panels/main", paste0(file_stem, ".pdf")),
  png = file.path("figures/final_panels/main", paste0(file_stem, ".png"))
)
readr::write_csv(manifest, "tables/main_single_panel_manifest.csv")
log_msg("Wrote tables/main_single_panel_manifest.csv")
log_msg("STOP: main single panels generated.")
