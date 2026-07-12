suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# Redraw Figure 05 from existing source data only.
# This fixes the over-strong "decreases" title without rerunning trajectory.

cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(cwd) == "scripts") {
  step_dir <- normalizePath(file.path(cwd, ".."), winslash = "/", mustWork = TRUE)
} else if (file.exists(file.path(cwd, "scripts", "_naming.R"))) {
  step_dir <- cwd
} else {
  step_dir <- normalizePath(file.path(cwd, "06_恶性细胞拟时序"), winslash = "/", mustWork = TRUE)
}
setwd(step_dir)
source(file.path("scripts", "_naming.R"))

run1_plot <- readr::read_csv(file.path("figures/source_data", "05_cytotrace2_pseudotime_cells.csv"), show_col_types = FALSE)
cyto_corr <- readr::read_csv(file.path("tables", "cytotrace2_pseudotime_correlation.csv"), show_col_types = FALSE)
short_colors <- setNames(subtype_naming_mapping$color, subtype_naming_mapping$abbreviation)

cyto_labels <- cyto_corr |>
  dplyr::mutate(
    lineage_label = factor(lineage_id, levels = paste0("lineage", 1:3), labels = paste0("L", 1:3)),
    label = paste0("rho = ", sprintf("%.2f", spearman_rho), "\np = ", format.pval(p_value, digits = 2, eps = 1e-300), "\nn = ", n_cells)
  )

p_cyto <- ggplot(run1_plot, aes(assigned_pseudotime, CytoTRACE2_score, color = subtype_short)) +
  geom_point(size = 0.12, alpha = 0.28) +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, linewidth = 0.55, color = "black", span = 0.65) +
  facet_wrap(~ lineage_label, scales = "free_x", nrow = 1) +
  geom_text(data = cyto_labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE, hjust = -0.05, vjust = 1.05, size = 2.1) +
  scale_color_manual(values = short_colors, name = NULL) +
  labs(
    title = "CytoTRACE2 vs primary pseudotime by lineage",
    subtitle = "Monotonic stemness loss is strongest for the OPC-M lineage (L3); MES-like lineages show non-monotonic ordering.",
    x = "Primary-run lineage pseudotime",
    y = "CytoTRACE2 score",
    caption = "L1/L2 trends should be interpreted cautiously because MES-like fine ordering is method-sensitive; monotonic stemness-loss claims are restricted to L3."
  ) +
  theme_bw(base_size = 7) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(size = 8),
    plot.subtitle = element_text(size = 6.4),
    plot.caption = element_text(size = 6, hjust = 0),
    strip.text = element_text(size = 7)
  )

pdf(file.path("figures", "05_cytotrace2_pseudotime.pdf"), width = 8.8, height = 3.65, useDingbats = FALSE)
print(p_cyto)
dev.off()
ggsave(file.path("figures", "05_cytotrace2_pseudotime_preview.png"), p_cyto, width = 8.8, height = 3.65, dpi = 180, bg = "white")

message("Redrew figures/05_cytotrace2_pseudotime.pdf and preview with corrected title.")
