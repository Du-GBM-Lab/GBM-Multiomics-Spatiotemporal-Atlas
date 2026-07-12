# 05_恶性细胞分亚群与Neftel对照/10_neftel_submodule_cycle_draft_figures.R
# Draft manuscript-style figures for Neftel state, submodule, and cycling audits.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

set.seed(42)

params <- list(
  project_dir = "05_恶性细胞分亚群与Neftel对照",
  input_object = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.submodule_labeled.qs2"
  ),
  tables_dir = file.path("05_恶性细胞分亚群与Neftel对照", "tables"),
  figures_dir = file.path("05_恶性细胞分亚群与Neftel对照", "figures"),
  plotdata_dir = file.path("05_恶性细胞分亚群与Neftel对照", "plotting_data"),
  run_tag = "20260517_draft10",
  subtype_col = "subtype_k4",
  sample_col = "Pt_number",
  primary_label_levels = paste0("Subtype", 1:4)
)

dir.create(params$figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(params$plotdata_dir, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.3, color = "black"),
      legend.title = element_text(size = base_size * 0.95),
      legend.text = element_text(size = base_size * 0.9),
      legend.key.size = grid::unit(0.36, "cm"),
      plot.title = element_text(face = "bold", size = base_size * 1.05),
      plot.margin = margin(6, 8, 6, 6)
    )
}

write_session_info <- function(path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  sink(con)
  print(sessionInfo())
  sink()
}

save_pdf <- function(plot, filename, width, height) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::pdf,
    useDingbats = FALSE
  )
}

state_levels <- c("AC", "OPC", "NPC", "MES")
submodule_levels <- c("AC", "OPC", "NPC1", "NPC2", "MES1", "MES2")
cycle_scores <- c("G1S_UCell", "G2M_UCell", "UCell_cycling")

msg("Input object:", params$input_object)
stopifnot(file.exists(params$input_object))

obj <- qs2::qs_read(params$input_object)
md <- obj@meta.data

required_cols <- c(
  params$subtype_col,
  params$sample_col,
  "neftel_state_UCell",
  "neftel_submodule_UCell",
  cycle_scores
)
missing_cols <- setdiff(required_cols, colnames(md))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

md <- md |>
  mutate(
    subtype_k4 = factor(as.character(.data[[params$subtype_col]]), levels = params$primary_label_levels),
    neftel_state_UCell = factor(as.character(neftel_state_UCell), levels = state_levels),
    neftel_submodule_UCell = factor(as.character(neftel_submodule_UCell), levels = submodule_levels)
  )

sanity_checks <- tibble(
  check = c(
    "n_cells",
    "n_features",
    "n_samples",
    "subtype_na",
    "state_na",
    "submodule_na",
    "cycle_score_na",
    "cycle_score_min",
    "cycle_score_max"
  ),
  value = c(
    ncol(obj),
    nrow(obj),
    n_distinct(md[[params$sample_col]]),
    sum(is.na(md$subtype_k4)),
    sum(is.na(md$neftel_state_UCell)),
    sum(is.na(md$neftel_submodule_UCell)),
    sum(is.na(md[, cycle_scores])),
    min(as.matrix(md[, cycle_scores]), na.rm = TRUE),
    max(as.matrix(md[, cycle_scores]), na.rm = TRUE)
  )
)

state_plot_data <- md |>
  count(subtype_k4, neftel_state_UCell, name = "n_cells") |>
  group_by(subtype_k4) |>
  mutate(row_pct = 100 * n_cells / sum(n_cells)) |>
  ungroup()

submodule_plot_data <- md |>
  count(subtype_k4, neftel_submodule_UCell, name = "n_cells") |>
  group_by(subtype_k4) |>
  mutate(row_pct = 100 * n_cells / sum(n_cells)) |>
  ungroup()

cycle_plot_data <- md |>
  select(subtype_k4, all_of(cycle_scores)) |>
  pivot_longer(
    cols = all_of(cycle_scores),
    names_to = "score",
    values_to = "value"
  ) |>
  mutate(
    score = recode(
      score,
      G1S_UCell = "G1/S",
      G2M_UCell = "G2/M",
      UCell_cycling = "Cycling"
    ),
    score = factor(score, levels = c("G1/S", "G2/M", "Cycling"))
  )

cycle_summary <- cycle_plot_data |>
  group_by(subtype_k4, score) |>
  summarise(
    n_cells = n(),
    n_missing = sum(is.na(value)),
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE, names = FALSE),
    q75 = quantile(value, 0.75, na.rm = TRUE, names = FALSE),
    iqr = IQR(value, na.rm = TRUE),
    .groups = "drop"
  )

plot_data <- list(
  params = params,
  sanity_checks = sanity_checks,
  state_plot_data = state_plot_data,
  submodule_plot_data = submodule_plot_data,
  cycle_plot_data = cycle_plot_data,
  cycle_summary = cycle_summary
)

plotdata_path <- file.path(params$plotdata_dir, paste0(params$run_tag, "_neftel_submodule_cycle_plotting_data.rds"))
saveRDS(plot_data, plotdata_path)

write_csv(sanity_checks, file.path(params$tables_dir, paste0(params$run_tag, "_neftel_cycle_figure_sanity_checks.csv")))
write_csv(state_plot_data, file.path(params$tables_dir, paste0(params$run_tag, "_neftel_state_heatmap_plot_data.csv")))
write_csv(submodule_plot_data, file.path(params$tables_dir, paste0(params$run_tag, "_neftel_submodule_heatmap_plot_data.csv")))
write_csv(cycle_summary, file.path(params$tables_dir, paste0(params$run_tag, "_cycling_score_summary.csv")))
write_session_info(file.path(params$tables_dir, paste0(params$run_tag, "_session_info.txt")))

state_heatmap <- ggplot(state_plot_data, aes(x = neftel_state_UCell, y = subtype_k4, fill = row_pct)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = sprintf("%.1f", row_pct)), size = 3.1, color = "black") +
  scale_fill_gradientn(
    colours = c("#F7FBFF", "#9ECAE1", "#3182BD", "#08519C"),
    limits = c(0, 100),
    name = "Cells (%)"
  ) +
  scale_x_discrete(drop = FALSE, position = "top") +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL, title = "Neftel 4-state dominance") +
  coord_fixed(ratio = 0.8) +
  theme_pub(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

submodule_heatmap <- ggplot(submodule_plot_data, aes(x = neftel_submodule_UCell, y = subtype_k4, fill = row_pct)) +
  geom_tile(color = "white", linewidth = 0.65) +
  geom_text(aes(label = sprintf("%.1f", row_pct)), size = 2.8, color = "black") +
  scale_fill_gradientn(
    colours = c("#F7FCF5", "#A1D99B", "#31A354", "#006D2C"),
    limits = c(0, 100),
    name = "Cells (%)"
  ) +
  scale_x_discrete(drop = FALSE, position = "top") +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL, title = "Neftel submodule dominance") +
  coord_fixed(ratio = 0.75) +
  theme_pub(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

subtype_palette <- c(
  Subtype1 = "#3C5488",
  Subtype2 = "#00A087",
  Subtype3 = "#F39B7F",
  Subtype4 = "#E64B35"
)

cycle_boxplot <- ggplot(cycle_plot_data, aes(x = subtype_k4, y = value, fill = subtype_k4)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.35) +
  facet_wrap(~score, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = subtype_palette, guide = "none") +
  labs(x = NULL, y = "UCell score", title = "Cell-cycle score distribution") +
  theme_pub(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

save_pdf(
  state_heatmap,
  file.path(params$figures_dir, paste0(params$run_tag, "_neftel_4state_heatmap.pdf")),
  width = 4.4,
  height = 3.2
)
save_pdf(
  submodule_heatmap,
  file.path(params$figures_dir, paste0(params$run_tag, "_neftel_submodule_heatmap.pdf")),
  width = 5.4,
  height = 3.2
)
save_pdf(
  cycle_boxplot,
  file.path(params$figures_dir, paste0(params$run_tag, "_cycling_score_boxplot.pdf")),
  width = 7.2,
  height = 3.6
)

msg("Sanity checks:")
print(sanity_checks)
msg("Plotting data:", plotdata_path)
msg("Draft PDFs written to:", params$figures_dir)
msg("Done.")
