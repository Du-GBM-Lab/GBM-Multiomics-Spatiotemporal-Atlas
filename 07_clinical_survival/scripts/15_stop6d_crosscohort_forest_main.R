# STOP6d: cross-cohort MES_V per-SD HR forest for the main R4 figure.
# Uses saved STOP5b cross-cohort summary; visualisation only.

suppressPackageStartupMessages({
  library(ggplot2)
})

root <- getwd()
fig_dir <- file.path(root, "figures")
tab_dir <- file.path(root, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(file.path(tab_dir, "STOP5b_crosscohort_MESV_direction_summary.csv"), stringsAsFactors = FALSE)

label_map <- c(
  CGGA_HGG_504 = "CGGA HGG discovery",
  CGGA_HGG_primary = "CGGA HGG primary",
  TCGA_HGG = "TCGA HGG validation",
  CGGA325_HGG = "CGGA-325 HGG validation",
  CGGA325_HGG_primary = "CGGA-325 primary"
)
df$label <- label_map[df$cohort]
df$label <- ifelse(is.na(df$label), df$cohort, df$label)
df$n_events <- paste0("n=", df$n, ", events=", df$events)
df$label_full <- paste0(df$label, "\n", df$n_events)

order_levels <- c("CGGA HGG discovery", "CGGA HGG primary", "TCGA HGG validation",
                  "CGGA-325 HGG validation", "CGGA-325 primary")
df$label <- factor(df$label, levels = rev(order_levels))
df$label_full <- factor(df$label_full, levels = rev(df$label_full[match(order_levels, as.character(df$label))]))

df$p_label <- ifelse(df$MESV_p < 1e-4, "p < 1e-4", paste0("p = ", signif(df$MESV_p, 2)))
df$hr_label <- sprintf("%.2f (%.2f-%.2f)", df$MESV_HR, df$lo, df$hi)
df$analysis_type <- ifelse(grepl("primary", df$cohort), "Primary-only sensitivity", "HGG landscape")

plot_df <- df
plot_df$y <- seq_len(nrow(plot_df))
plot_df <- plot_df[order(plot_df$label), , drop = FALSE]

p <- ggplot(plot_df, aes(MESV_HR, label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.5) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.16, orientation = "y", color = "grey25", linewidth = 0.7) +
  geom_point(aes(shape = analysis_type), size = 2.8, color = "#20854E") +
  geom_text(aes(x = hi + 0.08, label = hr_label), hjust = 0, size = 3.1, color = "grey20") +
  scale_shape_manual(values = c("HGG landscape" = 16, "Primary-only sensitivity" = 17), name = NULL) +
  scale_x_continuous(limits = c(0.8, 2.85), breaks = seq(0.8, 2.8, 0.4), expand = c(0, 0)) +
  labs(x = "Hazard ratio per 1 SD MES_V score (95% CI)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(6, 44, 6, 6)
  )

ggsave(file.path(fig_dir, "Fig_R4_crosscohort_MESV_forest.pdf"), p, width = 6.8, height = 3.8)

write.csv(plot_df[, c("cohort", "label", "n", "events", "MESV_HR", "lo", "hi", "MESV_p",
                      "MESV_zph_p", "argmax_p", "medianKM_p", "analysis_type", "hr_label", "p_label")],
          file.path(tab_dir, "STOP6d_crosscohort_MESV_forest_source.csv"), row.names = FALSE)

cat("\n########## STOP6d REPORT ##########\n")
cat("1. Forest written to figures/Fig_R4_crosscohort_MESV_forest.pdf.\n")
cat("2. Source written to tables/STOP6d_crosscohort_MESV_forest_source.csv.\n")
cat("3. HR range:", sprintf("%.2f-%.2f", min(plot_df$MESV_HR), max(plot_df$MESV_HR)), "\n")
cat("=> This is the cross-cohort direction-consistency panel; descriptive/prognostic association only.\n")
