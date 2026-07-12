suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

config <- list(
  annotated_qs2 = file.path("outputs", "GBM.RNA.qc_doubletfinder_annotated.qs2"),
  filtered_qs2 = file.path("outputs", "GBM.RNA.qc_doubletfinder.filtered.qs2"),
  table_dir = file.path(".", "tables"),
  figure_dir = file.path(".", "figures"),
  sample_col = "Pt_number",
  annotation_col = "anno_ident"
)

theme_qc_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size * 0.9),
      legend.key.size = grid::unit(0.38, "cm"),
      legend.position = "right",
      plot.margin = margin(6, 8, 6, 6)
    )
}

sample_order <- function(x) {
  x[order(as.integer(sub("^Pt", "", x)))]
}

dir.create(config$table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$figure_dir, showWarnings = FALSE, recursive = TRUE)

obj <- qs2::qs_read(config$annotated_qs2)
md <- obj@meta.data
md$sample <- as.character(md[[config$sample_col]])
md$annotation <- as.character(md[[config$annotation_col]])

sample_levels <- sample_order(unique(md$sample))

sample_retention <- md |>
  group_by(sample) |>
  summarise(
    n_input = n(),
    n_pass_pre_doublet = sum(qc_pass_pre_doublet, na.rm = TRUE),
    n_doublet = sum(doubletfinder_class == "Doublet", na.rm = TRUE),
    n_skipped_low_cells = sum(doubletfinder_class == "Skipped_low_cells", na.rm = TRUE),
    n_pass_final = sum(qc_pass_final, na.rm = TRUE),
    n_filtered_total = n_input - n_pass_final,
    retain_rate = n_pass_final / n_input * 100,
    pre_doublet_qc_rate = n_pass_pre_doublet / n_input * 100,
    doublet_rate_among_pre_qc = n_doublet / pmax(n_pass_pre_doublet, 1) * 100,
    .groups = "drop"
  ) |>
  mutate(sample = factor(sample, levels = sample_levels)) |>
  arrange(sample)

write.csv(
  sample_retention,
  file.path(config$table_dir, "sample_retention_summary.csv"),
  row.names = FALSE
)

p_retention <- ggplot(sample_retention, aes(x = sample, y = retain_rate)) +
  geom_col(width = 0.74, fill = "#2F6F73", color = "white", linewidth = 0.15) +
  geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.35, color = "#B24745") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Retained cells after QC (%)") +
  theme_qc_pub()

ggsave(
  file.path(config$figure_dir, "sample_retention_rate_after_QC.pdf"),
  p_retention,
  width = 8.4,
  height = 4.2,
  units = "in",
  device = cairo_pdf
)

celltype_stage_counts <- md |>
  mutate(qc_stage = if_else(qc_pass_final, "After_QC", "Before_QC")) |>
  count(qc_stage, annotation, name = "n_cells") |>
  group_by(qc_stage) |>
  mutate(percent = n_cells / sum(n_cells) * 100) |>
  ungroup()

before_counts <- md |>
  count(annotation, name = "n_before") |>
  mutate(percent_before = n_before / sum(n_before) * 100)

after_counts <- md |>
  filter(qc_pass_final) |>
  count(annotation, name = "n_after") |>
  mutate(percent_after = n_after / sum(n_after) * 100)

celltype_retention <- full_join(before_counts, after_counts, by = "annotation") |>
  mutate(
    n_after = if_else(is.na(n_after), 0L, n_after),
    percent_after = if_else(is.na(percent_after), 0, percent_after),
    n_filtered = n_before - n_after,
    retain_rate = n_after / n_before * 100,
    delta_percent_point = percent_after - percent_before
  ) |>
  arrange(retain_rate)

write.csv(
  celltype_retention,
  file.path(config$table_dir, "celltype_retention_before_after_QC.csv"),
  row.names = FALSE
)

celltype_levels <- celltype_retention$annotation
celltype_stage_counts$annotation <- factor(celltype_stage_counts$annotation, levels = celltype_levels)

p_celltype_retain <- ggplot(celltype_retention, aes(x = reorder(annotation, retain_rate), y = retain_rate)) +
  geom_col(width = 0.72, fill = "#3E6C8E", color = "white", linewidth = 0.12) +
  coord_flip() +
  geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.35, color = "#B24745") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Retained cells after QC (%)") +
  theme_qc_pub(base_size = 9) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(
  file.path(config$figure_dir, "celltype_retention_rate_after_QC.pdf"),
  p_celltype_retain,
  width = 6.2,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)

p_delta <- ggplot(celltype_retention, aes(x = reorder(annotation, delta_percent_point), y = delta_percent_point)) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey45") +
  geom_col(aes(fill = delta_percent_point > 0), width = 0.72, color = "white", linewidth = 0.12) +
  coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#2F6F73", `FALSE` = "#B24745")) +
  labs(x = NULL, y = "After - before QC composition (percentage points)") +
  theme_qc_pub(base_size = 9) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "none")

ggsave(
  file.path(config$figure_dir, "celltype_composition_delta_after_QC.pdf"),
  p_delta,
  width = 6.4,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)

message("Retention sanity checks completed.")
