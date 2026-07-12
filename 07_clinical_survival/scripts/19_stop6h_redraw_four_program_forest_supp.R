# STOP6h: redraw four-program per-SD Cox forest as supplementary figure.
# Shows all program directions; not the main MES_V validation panel.

suppressPackageStartupMessages({
  library(ggplot2)
})

root <- getwd()
fig_dir <- file.path(root, "figures")
tab_dir <- file.path(root, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

PROGS <- c("NPC_P", "OPC_M", "MES_V", "MES_I")
LABELS <- c(NPC_P = "NPC-P", OPC_M = "OPC-M", MES_V = "MES-V", MES_I = "MES-I")

cox <- read.csv(file.path(tab_dir, "STOP4_perSD_Cox_results.csv"), stringsAsFactors = FALSE)
df <- cox[cox$model %in% c("four_program_idh_grade", "four_program_age_sex") &
            cox$score %in% PROGS, , drop = FALSE]

df$program <- factor(LABELS[df$score], levels = rev(LABELS[PROGS]))
df$analysis <- ifelse(df$set == "all504",
                      "CGGA HGG\nadjusted for IDH/grade",
                      "Primary IDH-wt WHO IV GBM\nadjusted for age/sex")
df$analysis <- factor(df$analysis, levels = c("CGGA HGG\nadjusted for IDH/grade",
                                              "Primary IDH-wt WHO IV GBM\nadjusted for age/sex"))
df$hr_label <- sprintf("%.2f (%.2f-%.2f)", df$HR, df$lo, df$hi)
df$p_label <- ifelse(df$p < 1e-4, "p < 1e-4", paste0("p = ", signif(df$p, 2)))

write.csv(df, file.path(tab_dir, "STOP6h_four_program_forest_source.csv"), row.names = FALSE)

p <- ggplot(df, aes(HR, program)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.5) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.16, orientation = "y", color = "grey25", linewidth = 0.7) +
  geom_point(size = 2.6, color = "#0072B5") +
  facet_wrap(~ analysis, nrow = 1, scales = "free_x") +
  labs(x = "Hazard ratio per 1 SD score (95% CI)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 10),
    panel.spacing = grid::unit(8, "mm"),
    plot.margin = margin(6, 6, 6, 6)
  )

ggsave(file.path(fig_dir, "FigS_R4_four_program_forest_clean.pdf"), p, width = 7.2, height = 3.6)

cat("\n########## STOP6h REPORT ##########\n")
cat("1. Clean four-program forest written to figures/FigS_R4_four_program_forest_clean.pdf.\n")
cat("2. Source written to tables/STOP6h_four_program_forest_source.csv.\n")
print(df[, c("analysis", "program", "n", "events", "HR", "lo", "hi", "p", "zph_p", "covars")])
cat("=> Supplementary defense figure; do not overstate primary IDH-wt GBM sensitivity.\n")
