# S1_k_selection_diagnostics_theta5.R
# Purpose: Draw reviewer-facing diagnostics supporting the coarse-grained
# k=4 malignant subtype framework.
#
# Input:
#   tables/02c_sweep_theta5/gap_curve.csv
#   tables/02c_sweep_theta5/silhouette.csv
#   tables/02c_sweep_theta5/subtype_risk_per_k.csv
#
# Output:
#   figures/S1_k_selection_diagnostics_theta5.pdf
#
# Interpretation boundary:
#   Gap statistic and silhouette are diagnostic metrics, not automatic decision
#   rules. Do not write "gap/silhouette determined k=4".

.libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

base_dir <- "05_恶性细胞分亚群与Neftel对照"
tab_dir <- file.path(base_dir, "tables/02c_sweep_theta5")
fig_dir <- file.path(base_dir, "figures")
source_dir <- file.path(fig_dir, "source_data")
dir.create(source_dir, showWarnings = FALSE, recursive = TRUE)

gap_path <- file.path(tab_dir, "gap_curve.csv")
sil_path <- file.path(tab_dir, "silhouette.csv")
risk_path <- file.path(tab_dir, "subtype_risk_per_k.csv")
out_pdf <- file.path(fig_dir, "S1_k_selection_diagnostics_theta5.pdf")

gap_df <- read.csv(gap_path, check.names = FALSE)
sil_df <- read.csv(sil_path, check.names = FALSE)
risk_df <- read.csv(risk_path, check.names = FALSE)

stopifnot(all(c("k", "gap", "SE.sim") %in% colnames(gap_df)))
stopifnot(all(c("k", "avg_sil") %in% colnames(sil_df)))
stopifnot(all(c("k", "subtype", "top_patient_pct", "cycling_pct") %in% colnames(risk_df)))
stopifnot(any(gap_df$k == 4), any(sil_df$k == 4), any(risk_df$k == 4))
stopifnot(!anyNA(gap_df[, c("k", "gap", "SE.sim")]))
stopifnot(!anyNA(sil_df[, c("k", "avg_sil")]))
stopifnot(!anyNA(risk_df[, c("k", "top_patient_pct", "cycling_pct")]))

highlight_col <- "#BC3C29"
line_col <- "grey25"
bar_col <- "grey55"

theme_supp <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      axis.title = element_text(size = base_size, color = "black"),
      axis.text = element_text(size = base_size - 1, color = "black"),
      strip.background = element_rect(fill = "grey92", color = "grey55", linewidth = 0.3),
      strip.text = element_text(size = base_size - 1, face = "bold", color = "black"),
      legend.position = "none",
      plot.margin = margin(4, 4, 4, 4, "pt")
    )
}

p_gap <- ggplot(gap_df, aes(x = k, y = gap)) +
  geom_errorbar(aes(ymin = gap - SE.sim, ymax = gap + SE.sim), width = 0.15, linewidth = 0.35, color = "grey55") +
  geom_line(linewidth = 0.45, color = line_col) +
  geom_point(size = 1.7, color = line_col) +
  geom_vline(xintercept = 4, linetype = "dashed", linewidth = 0.35, color = highlight_col) +
  geom_point(data = gap_df %>% filter(k == 4), size = 2.5, color = highlight_col) +
  scale_x_continuous(breaks = gap_df$k) +
  labs(title = "A. Gap statistic", x = "Number of subtypes (k)", y = "Gap statistic") +
  theme_supp()

p_sil <- ggplot(sil_df, aes(x = k, y = avg_sil)) +
  geom_line(linewidth = 0.45, color = line_col) +
  geom_point(size = 1.7, color = line_col) +
  geom_vline(xintercept = 4, linetype = "dashed", linewidth = 0.35, color = highlight_col) +
  geom_point(data = sil_df %>% filter(k == 4), size = 2.5, color = highlight_col) +
  scale_x_continuous(breaks = sil_df$k) +
  labs(title = "B. Average silhouette", x = "Number of subtypes (k)", y = "Average silhouette") +
  theme_supp()

risk_plot_df <- risk_df %>%
  mutate(
    k = factor(k, levels = sort(unique(k))),
    subtype = factor(subtype, levels = paste0("Subtype", 1:5))
  )

p_patient <- ggplot(risk_plot_df, aes(x = subtype, y = top_patient_pct)) +
  geom_col(fill = bar_col, width = 0.72) +
  geom_hline(yintercept = 70, linetype = "dotted", color = highlight_col, linewidth = 0.35) +
  facet_grid(. ~ k, scales = "free_x", space = "free_x", labeller = label_both) +
  scale_y_continuous(labels = label_number(suffix = "%"), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "C. Patient dominance", x = NULL, y = "Top patient fraction") +
  theme_supp() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

p_cycle <- ggplot(risk_plot_df, aes(x = subtype, y = cycling_pct)) +
  geom_col(fill = bar_col, width = 0.72) +
  facet_grid(. ~ k, scales = "free_x", space = "free_x", labeller = label_both) +
  scale_y_continuous(labels = label_number(suffix = "%"), expand = expansion(mult = c(0, 0.08))) +
  labs(title = "D. Cycling fraction", x = NULL, y = "Cycling fraction") +
  theme_supp() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

fig <- (p_gap | p_sil) / (p_patient | p_cycle) +
  plot_annotation(
    title = "Diagnostics supporting a coarse-grained k=4 subtype framework",
    caption = "Subtype labels within each k panel represent independent partitions and are not directly comparable across k values.",
    theme = theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.02),
      plot.caption = element_text(size = 7.5, color = "grey25", hjust = 0.02)
    )
  )

ggsave(out_pdf, fig, width = 8.0, height = 6.6, device = cairo_pdf)

write.csv(gap_df, file.path(source_dir, "S1_k_selection_gap_curve.csv"), row.names = FALSE)
write.csv(sil_df, file.path(source_dir, "S1_k_selection_silhouette.csv"), row.names = FALSE)
write.csv(risk_df, file.path(source_dir, "S1_k_selection_subtype_risk.csv"), row.names = FALSE)

cat("Output PDF:", out_pdf, "\n")
cat("Source data directory:", source_dir, "\n")
cat("k=4 gap:", gap_df$gap[gap_df$k == 4], "\n")
cat("k=4 silhouette:", sil_df$avg_sil[sil_df$k == 4], "\n")
cat("Sanity checks passed.\n")
