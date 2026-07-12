options(stringsAsFactors = FALSE)
set.seed(42)

if (basename(getwd()) == "06_恶性细胞拟时序") {
  project_root <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(file.path(project_root, "06_恶性细胞拟时序"))

dir.create("figures/final_panels/supplement", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/final_panels/source_data/supplement", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggridges)
})

log_file <- "logs/11_make_supp_single_panels.log"
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(...))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_msg("Supplement single-panel figure generation started.")

pal <- c("NPC-P" = "#0072B5", "OPC-M" = "#E18727", "MES-V" = "#20854E", "MES-I" = "#BC3C29")
lineage_labs <- c("lineage1" = "L1: NPC-P to MES-I", "lineage2" = "L2: NPC-P to MES-V", "lineage3" = "L3: NPC-P to OPC-M")
run_labels <- c("run1_NPCP" = "NPC-P", "run2_OPCM" = "OPC-M", "run3_MESV" = "MES-V", "run4_auto" = "Auto")
mp_labels <- c(
  "MP01" = "MP01 IGFBP/stress",
  "MP02" = "MP02 ECM/mural",
  "MP03" = "MP03 myelination",
  "MP04" = "MP04 MHC-II",
  "MP05" = "MP05 cycling",
  "MP06" = "MP06 endothelial"
)
mp_colors <- c(
  "MP01" = "#E64B35",
  "MP02" = "#7E6148",
  "MP03" = "#00A087",
  "MP04" = "#4DBBD5",
  "MP05" = "#3C5488",
  "MP06" = "#F39B7F"
)
ams_colors <- c(
  "AC" = "#E64B35",
  "MES1" = "#7E6148",
  "MES2" = "#00A087",
  "NPC1" = "#4DBBD5",
  "NPC2" = "#3C5488",
  "OPC" = "#F564E3"
)

theme_supp <- function(base_size = 7) {
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
      legend.text = element_text(size = base_size - 0.4),
      strip.text = element_text(size = base_size, face = "bold")
    )
}

save_panel <- function(plot, stem, width, height) {
  pdf_path <- file.path("figures/final_panels/supplement", paste0(stem, ".pdf"))
  png_path <- file.path("figures/final_panels/supplement", paste0(stem, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(png_path, plot, width = width, height = height, dpi = 400, bg = "white")
  log_msg("Wrote", pdf_path)
  log_msg("Wrote", png_path)
}

copy_source <- function(src, dest_name = basename(src)) {
  dest <- file.path("figures/final_panels/source_data/supplement", dest_name)
  file.copy(src, dest, overwrite = TRUE)
  log_msg("Copied source data", dest)
}

copy_existing_panel <- function(stem, new_stem) {
  for (ext in c("pdf", "png")) {
    src <- file.path("figures", paste0(stem, ".", ext))
    dest <- file.path("figures/final_panels/supplement", paste0(new_stem, ".", ext))
    if (file.exists(src)) {
      file.copy(src, dest, overwrite = TRUE)
      log_msg("Copied panel", dest)
    }
  }
}

# S1: root justification full metrics
cell_cycle <- read_csv("figures/source_data/01_cell_cycle_phase_composition.csv", show_col_types = FALSE) |>
  mutate(subtype_short = factor(subtype_short, levels = names(pal)), Phase = factor(Phase, levels = c("G1", "S", "G2M")))
copy_source("figures/source_data/01_cell_cycle_phase_composition.csv")
p_s1a <- ggplot(cell_cycle, aes(subtype_short, prop, fill = Phase)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = function(x) paste0(round(x * 100), "%")) +
  scale_fill_manual(values = c("G1" = "#3C5488", "S" = "#4DBBD5", "G2M" = "#E64B35")) +
  labs(title = "S1A. Cell-cycle phase composition", x = NULL, y = "Cell fraction") +
  theme_supp(7)
save_panel(p_s1a, "S1A_cell_cycle_phase", 3.5, 2.8)

sox <- read_csv("figures/source_data/01_SOX_GAP43_expression.csv", show_col_types = FALSE) |>
  filter(gene %in% c("SOX2", "SOX4", "SOX9", "SOX11")) |>
  mutate(subtype_short = factor(subtype_short, levels = names(pal)), gene = factor(gene, levels = c("SOX2", "SOX4", "SOX9", "SOX11")))
copy_source("figures/source_data/01_SOX_GAP43_expression.csv")
sox_dot <- sox |>
  group_by(subtype_short, gene) |>
  summarise(
    mean_expression = mean(expression, na.rm = TRUE),
    pct_expressing = mean(expression > 0, na.rm = TRUE) * 100,
    n_cells = n(),
    .groups = "drop"
  )
readr::write_csv(sox_dot, "figures/final_panels/source_data/supplement/S1B_SOX_dotplot_summary.csv")
p_s1b <- ggplot(sox_dot, aes(gene, subtype_short)) +
  geom_point(aes(size = pct_expressing, fill = mean_expression), shape = 21, color = "black", stroke = 0.18) +
  scale_size_area(max_size = 7, breaks = c(20, 40, 60), name = "Positive cells (%)") +
  scale_fill_gradient(low = "#F4F6F7", high = "#0072B5", name = "Mean expression") +
  labs(title = "S1B. SOX family expression", x = NULL, y = NULL) +
  theme_supp(7) +
  theme(
    axis.text.x = element_text(size = 6.5, face = "bold"),
    axis.text.y = element_text(size = 6.5),
    legend.position = "right",
    panel.grid.major = element_line(linewidth = 0.18, color = "grey88")
  )
save_panel(p_s1b, "S1B_SOX_expression", 4.2, 2.8)

neftel <- read_csv("figures/source_data/01_Neftel_6state_AMS.csv", show_col_types = FALSE) |>
  mutate(subtype_short = factor(subtype_short, levels = names(pal)), module = factor(module, levels = c("AMS_NPC1", "AMS_NPC2", "AMS_OPC", "AMS_MES1", "AMS_MES2", "AMS_AC")))
copy_source("figures/source_data/01_Neftel_6state_AMS.csv")
neftel_heat <- neftel |>
  group_by(subtype_short, module) |>
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    pct_positive = mean(score > 0, na.rm = TRUE) * 100,
    n_cells = n(),
    .groups = "drop"
  )
readr::write_csv(neftel_heat, "figures/final_panels/source_data/supplement/S1C_Neftel_AMS_heatmap_summary.csv")
p_s1c <- ggplot(neftel_heat, aes(module, subtype_short, fill = mean_score)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", mean_score)), size = 2.2) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    name = "Mean AMS"
  ) +
  labs(title = "S1C. Neftel six-state AMS alignment", x = NULL, y = NULL) +
  theme_supp(7) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6.0, face = "bold"), axis.text.y = element_text(size = 6.5))
save_panel(p_s1c, "S1C_Neftel_AMS", 4.8, 2.8)

# S2: Monocle3 validation
m3_cells <- read_csv("figures/source_data/03_panel_A_B_cells.csv", show_col_types = FALSE) |>
  mutate(subtype_short = factor(subtype_short, levels = names(pal)))
m3_edges <- read_csv("figures/source_data/03_panel_A_graph_edges.csv", show_col_types = FALSE)
m3_nodes <- read_csv("figures/source_data/03_panel_A_graph_nodes.csv", show_col_types = FALSE)
copy_source("figures/source_data/03_panel_A_B_cells.csv")
copy_source("figures/source_data/03_panel_A_graph_edges.csv")
copy_source("figures/source_data/03_panel_A_graph_nodes.csv")
node_pst <- lapply(seq_len(nrow(m3_nodes)), function(i) {
  dx <- m3_cells$UMAP_1 - m3_nodes$UMAP_1[i]
  dy <- m3_cells$UMAP_2 - m3_nodes$UMAP_2[i]
  nearest <- order(dx * dx + dy * dy)[seq_len(min(50, nrow(m3_cells)))]
  tibble(
    node_id = m3_nodes$node_id[i],
    node_pst = median(m3_cells$monocle3_pst[nearest], na.rm = TRUE)
  )
}) |>
  bind_rows()
m3_edges_directed <- m3_edges |>
  left_join(node_pst, by = c("from_node" = "node_id")) |>
  rename(from_pst = node_pst) |>
  left_join(node_pst, by = c("to_node" = "node_id")) |>
  rename(to_pst = node_pst) |>
  mutate(
    x = ifelse(from_pst <= to_pst, from_UMAP_1, to_UMAP_1),
    y = ifelse(from_pst <= to_pst, from_UMAP_2, to_UMAP_2),
    xend = ifelse(from_pst <= to_pst, to_UMAP_1, from_UMAP_1),
    yend = ifelse(from_pst <= to_pst, to_UMAP_2, from_UMAP_2)
  )
m3_edges_arrow <- m3_edges_directed |>
  mutate(
    edge_len = sqrt((xend - x)^2 + (yend - y)^2),
    pst_delta = abs(to_pst - from_pst)
  ) |>
  filter(edge_len >= 0.45 | pst_delta >= quantile(pst_delta, 0.70, na.rm = TRUE)) |>
  arrange(desc(edge_len + pst_delta)) |>
  slice(seq(1, n(), by = 3))
readr::write_csv(node_pst, "figures/final_panels/source_data/supplement/S2A_monocle3_node_pseudotime.csv")
readr::write_csv(m3_edges_directed, "figures/final_panels/source_data/supplement/S2A_monocle3_directed_graph_edges.csv")
readr::write_csv(m3_edges_arrow, "figures/final_panels/source_data/supplement/S2A_monocle3_arrow_edges.csv")
root_node <- m3_nodes |> filter(node_role == "root")
terminal_nodes <- m3_nodes |> filter(node_role == "terminal")
subtype_labels <- m3_cells |>
  group_by(subtype_short) |>
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    UMAP_2 = case_when(
      subtype_short == "NPC-P" ~ UMAP_2 - 1.25,
      subtype_short == "OPC-M" ~ UMAP_2 + 0.95,
      subtype_short == "MES-V" ~ UMAP_2 - 1.00,
      subtype_short == "MES-I" ~ UMAP_2 + 1.00,
      TRUE ~ UMAP_2
    )
  )
readr::write_csv(subtype_labels, "figures/final_panels/source_data/supplement/S2A_subtype_label_positions.csv")
p_s2a <- ggplot() +
  geom_point(data = m3_cells, aes(UMAP_1, UMAP_2, color = monocle3_pst), size = 0.075, alpha = 0.60) +
  geom_segment(
    data = m3_edges_directed,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.70,
    color = "white",
    alpha = 0.90,
    lineend = "round"
  ) +
  geom_segment(
    data = m3_edges_directed,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.34,
    color = "black",
    alpha = 0.92,
    lineend = "round"
  ) +
  geom_segment(
    data = m3_edges_arrow,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.45,
    color = "black",
    alpha = 0.95,
    lineend = "round",
    arrow = arrow(type = "closed", length = unit(0.052, "inches"))
  ) +
  geom_point(data = terminal_nodes, aes(UMAP_1, UMAP_2), shape = 24, fill = "black", color = "black", size = 1.15) +
  geom_point(data = root_node, aes(UMAP_1, UMAP_2), shape = 21, fill = "white", color = "black", size = 1.35, stroke = 0.35) +
  geom_label(
    data = root_node,
    aes(UMAP_1, UMAP_2 - 1.0, label = "Root"),
    size = 1.7,
    linewidth = 0.12,
    label.padding = unit(0.08, "lines"),
    fill = "white",
    color = "black"
  ) +
  ggrepel::geom_label_repel(
    data = subtype_labels,
    aes(UMAP_1, UMAP_2, label = subtype_short),
    size = 2.0,
    linewidth = 0.14,
    label.padding = unit(0.10, "lines"),
    fill = "white",
    color = "black",
    min.segment.length = 0,
    segment.size = 0.12,
    max.overlaps = Inf
  ) +
  scale_color_viridis_c(name = "Monocle3\npseudotime", option = "viridis", na.value = "grey85") +
  labs(
    title = "S2A. Monocle3 pseudotime with principal graph",
    subtitle = "Arrows follow increasing local Monocle3 pseudotime; terminal nodes shown as triangles"
  ) +
  coord_equal() +
  theme_umap(7) +
  theme(legend.position = "right")
save_panel(p_s2a, "S2A_monocle3_graph", 4.2, 3.6)

p_s2b <- ggplot(m3_cells, aes(UMAP_1, UMAP_2, color = monocle3_pst)) +
  geom_point(size = 0.07, alpha = 0.55) +
  scale_color_viridis_c(name = "Monocle3\npseudotime", option = "viridis", na.value = "grey85") +
  labs(title = "S2B. Monocle3 pseudotime") +
  coord_equal() +
  theme_umap(7) +
  theme(legend.position = "right")
save_panel(p_s2b, "S2B_monocle3_pseudotime", 4.2, 3.6)

m3_corr <- read_csv("figures/source_data/03_panel_C.csv", show_col_types = FALSE) |>
  mutate(lineage_id = factor(lineage_id, levels = c("L1", "L2", "L3")), subtype_short = factor(subtype_short, levels = names(pal)))
copy_source("figures/source_data/03_panel_C.csv")
m3_corr_stats <- read_csv("tables/pseudotime_method_correlation.csv", show_col_types = FALSE) |>
  mutate(
    lineage_id = recode(lineage_id, lineage1 = "L1", lineage2 = "L2", lineage3 = "L3"),
    lineage_id = factor(lineage_id, levels = c("L1", "L2", "L3")),
    label = paste0("rho = ", sprintf("%.2f", spearman_rho), "\n", "n = ", format(n_cells, big.mark = ","))
  )
stats_pos <- m3_corr |>
  group_by(lineage_id) |>
  summarise(
    x = quantile(slingshot_pst, 0.05, na.rm = TRUE),
    y = quantile(monocle3_pst, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(m3_corr_stats, by = "lineage_id")
readr::write_csv(m3_corr_stats, "figures/final_panels/source_data/supplement/S2B_method_correlation_stats.csv")
p_s2c <- ggplot(m3_corr, aes(slingshot_pst, monocle3_pst, color = subtype_short)) +
  geom_point(size = 0.08, alpha = 0.18) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "black", linewidth = 0.42, span = 0.7) +
  geom_label(
    data = stats_pos,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.0,
    linewidth = 0.12,
    label.padding = unit(0.10, "lines"),
    fill = "white",
    color = "black",
    hjust = 0
  ) +
  facet_wrap(~lineage_id, scales = "free", nrow = 1) +
  scale_color_manual(values = pal, name = NULL) +
  labs(title = "S2B. Cross-method pseudotime concordance", x = "Slingshot pseudotime", y = "Monocle3 pseudotime") +
  theme_supp(7) +
  theme(legend.position = "bottom")
save_panel(p_s2c, "S2C_method_pseudotime_correlation", 6.2, 2.7)

m3_flow <- read_csv("figures/source_data/03_panel_D.csv", show_col_types = FALSE) |>
  mutate(slingshot_assigned_lineage = recode(slingshot_assigned_lineage, lineage1 = "L1", lineage2 = "L2", lineage3 = "L3"))
terminal_labels <- read_csv("tables/monocle3_terminal_detection.csv", show_col_types = FALSE) |>
  group_by(monocle3_terminal) |>
  slice_max(terminal_fraction, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    terminal_group = if_else(terminal_fraction < 0.6, "mixed", subtype_short),
    terminal_label = paste0(monocle3_terminal, "\n", terminal_group, " ", round(terminal_fraction * 100), "%")
  ) |>
  select(monocle3_terminal, terminal_label)
terminal_order <- c("Y_26", "Y_70", "Y_72", "Y_41", "Y_43", "Y_6", "Y_98", "Y_31")
terminal_labels <- terminal_labels |>
  mutate(monocle3_terminal = factor(monocle3_terminal, levels = terminal_order)) |>
  arrange(monocle3_terminal)
m3_flow <- m3_flow |>
  left_join(terminal_labels, by = "monocle3_terminal") |>
  mutate(
    terminal_label = factor(terminal_label, levels = terminal_labels$terminal_label),
    slingshot_assigned_lineage = factor(slingshot_assigned_lineage, levels = c("L1", "L2", "L3"))
  )
copy_source("figures/source_data/03_panel_D.csv")
p_s2d <- ggplot(m3_flow, aes(terminal_label, slingshot_assigned_lineage, fill = as.numeric(n_cells))) +
  geom_tile(color = "white", linewidth = 0.25) +
  annotate("rect", xmin = 6.5, xmax = 7.5, ymin = 1.5, ymax = 2.5, fill = NA, color = "#20854E", linewidth = 0.45, linetype = "22") +
  annotate("rect", xmin = 7.5, xmax = 8.5, ymin = 2.5, ymax = 3.5, fill = NA, color = "#E18727", linewidth = 0.45, linetype = "22") +
  geom_vline(xintercept = c(6.5, 7.5), linetype = "22", color = "grey45", linewidth = 0.25) +
  geom_text(aes(label = n_cells), size = 2.0) +
  scale_fill_gradient(name = "Cells", low = "#F7FBFF", high = "#F16913") +
  labs(title = "S2C. Slingshot lineage versus Monocle3 terminal", x = "Monocle3 terminal", y = "Slingshot lineage") +
  theme_supp(7) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p_s2d, "S2C_lineage_terminal_heatmap", 5.6, 2.9)

# S3: root sensitivity correlation heatmap
root_cor <- read_csv("figures/source_data/04_panel_E_correlation_heatmap.csv", show_col_types = FALSE) |>
  mutate(run_x = factor(run_labels[run_x], levels = run_labels), run_y = factor(run_labels[run_y], levels = run_labels))
copy_source("figures/source_data/04_panel_E_correlation_heatmap.csv")
p_s3a <- ggplot(root_cor, aes(run_x, run_y, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 2.1) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = c(-1, 1), name = "Spearman\nrho") +
  labs(title = "S3A. Root-relative pseudotime correlation", subtitle = "Low/negative values are expected when roots are reversed", x = NULL, y = NULL) +
  theme_supp(7) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p_s3a, "S3A_root_relative_correlation", 3.4, 3.0)

# S4: CytoTRACE2 versus pseudotime
cyto_pst <- read_csv("figures/source_data/05_cytotrace2_pseudotime_cells.csv", show_col_types = FALSE) |>
  mutate(lineage_label = factor(lineage_label, levels = c("L1", "L2", "L3")), subtype_short = factor(subtype_short, levels = names(pal)))
copy_source("figures/source_data/05_cytotrace2_pseudotime_cells.csv")
p_s4a <- ggplot(cyto_pst, aes(assigned_pseudotime, CytoTRACE2_score, color = subtype_short)) +
  geom_point(size = 0.08, alpha = 0.20) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.45, span = 0.6) +
  facet_wrap(~lineage_label, scales = "free_x", nrow = 1) +
  scale_color_manual(values = pal, name = NULL) +
  labs(title = "S4A. CytoTRACE2 versus primary pseudotime", x = "Assigned pseudotime", y = "CytoTRACE2 score") +
  theme_supp(7) +
  theme(legend.position = "bottom")
save_panel(p_s4a, "S4A_cytotrace2_pseudotime", 6.2, 2.8)

# S5: validation panels already single
validation_map <- c(
  "07_CNV比例" = "S5A_CNV_positive_fraction",
  "07_doublet分布" = "S5B_doublet_score",
  "07_共表达比例" = "S5C_marker_coexpression",
  "07_MHCII_vs_TAM" = "S5D_MHCII_vs_TAM",
  "07_血管marker_vs_mural" = "S5E_vascular_marker_vs_mural"
)
for (i in seq_along(validation_map)) {
  copy_existing_panel(names(validation_map)[i], validation_map[[i]])
}
for (src in c("07_CNV比例.csv", "07_doublet分布.csv", "07_共表达比例.csv", "07_MHCII_vs_TAM.csv", "07_血管marker_vs_mural.csv", "07_marker_summary_all_groups.csv")) {
  copy_source(file.path("figures/source_data", src))
}

# S6: biology overlay
mp <- read_csv("figures/source_data/06_panel_A_metaprogram.csv", show_col_types = FALSE) |>
  mutate(lineage_label = factor(lineage_label, levels = c("L1", "L2", "L3")), score_label = factor(score_label, levels = paste0("MP", sprintf("%02d", 1:6))))
copy_source("figures/source_data/06_panel_A_metaprogram.csv")
p_s6a <- ggplot(mp, aes(pseudotime_mean, score_mean, color = score_label, group = score_label)) +
  geom_line(linewidth = 0.62) +
  facet_wrap(~lineage_label, scales = "free_x", nrow = 1) +
  scale_color_manual(values = mp_colors, labels = mp_labels, drop = FALSE) +
  labs(title = "S6A. NMF metaprograms along pseudotime", x = "Pseudotime bin mean", y = "Mean score", color = "Program") +
  theme_supp(7) +
  theme(legend.position = "bottom", legend.key.width = unit(0.9, "lines"))
save_panel(p_s6a, "S6A_metaprogram_pseudotime", 6.8, 3.05)

ams <- read_csv("figures/source_data/06_panel_B_neftel.csv", show_col_types = FALSE) |>
  mutate(lineage_label = factor(lineage_label, levels = c("L1", "L2", "L3")))
copy_source("figures/source_data/06_panel_B_neftel.csv")
p_s6b <- ggplot(ams, aes(pseudotime_mean, score_mean, color = score_label, group = score_label)) +
  geom_line(linewidth = 0.62) +
  facet_wrap(~lineage_label, scales = "free_x", nrow = 1) +
  scale_color_manual(values = ams_colors, drop = FALSE) +
  labs(title = "S6B. Neftel AMS scores along pseudotime", x = "Pseudotime bin mean", y = "Mean AMS score", color = "State") +
  theme_supp(7) +
  theme(legend.position = "bottom", legend.key.width = unit(0.9, "lines"))
save_panel(p_s6b, "S6B_neftel_pseudotime", 6.8, 3.05)

mp_mes <- read_csv("figures/source_data/06_panel_C_MP02_MP04.csv", show_col_types = FALSE) |>
  mutate(
    lineage_label = factor(lineage_label, levels = c("L1", "L2")),
    score_label = factor(score_label, levels = c("MP02", "MP04"))
  )
copy_source("figures/source_data/06_panel_C_MP02_MP04.csv")
p_s6c <- ggplot(mp_mes, aes(pseudotime_mean, score_mean, color = score_label, linetype = lineage_label, group = interaction(lineage_label, score_label))) +
  geom_line(linewidth = 0.55) +
  scale_color_manual(values = c("MP02" = "#20854E", "MP04" = "#BC3C29")) +
  labs(title = "S6C. MP02 and MP04 across MES-like lineages", x = "Pseudotime bin mean", y = "Mean score", color = NULL, linetype = NULL) +
  theme_supp(7) +
  theme(legend.position = "bottom")
save_panel(p_s6c, "S6C_MP02_MP04_MES_lineages", 4.6, 3.2)

# S7: driver panels already single
copy_existing_panel("08_L3_driver_heatmap", "S7A_L3_driver_heatmap")
copy_existing_panel("08_driver_enrichment_dotplot", "S7B_driver_GO_dotplot")
copy_source("figures/source_data/08_L3_driver_heatmap.csv")
copy_source("figures/source_data/08_driver_enrichment_dotplot.csv")
copy_source("tables/driver_enrichment_background.csv")

manifest <- tibble(
  file_stem = tools::file_path_sans_ext(basename(list.files("figures/final_panels/supplement", pattern = "\\.pdf$", full.names = FALSE))),
  pdf = file.path("figures/final_panels/supplement", paste0(file_stem, ".pdf")),
  png = file.path("figures/final_panels/supplement", paste0(file_stem, ".png"))
) |>
  arrange(file_stem)
write_csv(manifest, "tables/supp_single_panel_manifest.csv")
log_msg("Wrote tables/supp_single_panel_manifest.csv")
log_msg("STOP: supplement single panels generated.")
