# =============================================================================
# R2 · scWGCNA —— STAGE 8: specificity plots + PLAUR-in-MES-V visualization
# -----------------------------------------------------------------------------
# Read-only plotting: read STAGE 6 object and existing csv files, do not
# recalculate network, modify object, or save a new object.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(hdWGCNA)
  library(tidyverse)
  library(qs2)
  library(patchwork)
})

OUT_DIR <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/R2_scWGCNA"
OBJ_S6 <- file.path(OUT_DIR, "obj_malignant_R2_stage6_ME_subtype.qs2")
fig_dir <- file.path(OUT_DIR, "figures")
tab_dir <- file.path(OUT_DIR, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

SUBTYPE_LEVELS <- c("NPC-P", "OPC-M", "MES-V", "MES-I")
SUBTYPE_COL <- c(
  "NPC-P" = "#0072B5",
  "OPC-M" = "#E18727",
  "MES-V" = "#20854E",
  "MES-I" = "#BC3C29"
)
MESV_MOD <- "green"

obj <- qs2::qs_read(OBJ_S6)
obj$subtype_bio <- factor(obj$subtype_bio, levels = SUBTYPE_LEVELS)

tau <- read.csv(file.path(tab_dir, "module_specificity_tau.csv"), stringsAsFactors = FALSE)
tau$top_subtype <- factor(tau$top_subtype, levels = SUBTYPE_LEVELS)
tau$module <- factor(tau$module, levels = tau$module[order(tau$tau)])

p_tau <- ggplot(tau, aes(tau, module, fill = top_subtype)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0.7, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = SUBTYPE_COL, name = "Assigned subtype") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(c(0, 0.02))) +
  labs(x = "Module specificity (tau)", y = NULL, title = "Module-subtype specificity") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")
ggsave(file.path(fig_dir, "S_module_tau_barplot.pdf"), p_tau, width = 6, height = 5)

hMEs <- GetMEs(obj, harmonized = TRUE)
stopifnot(MESV_MOD %in% colnames(hMEs))
df <- data.frame(subtype = obj$subtype_bio, hME = hMEs[[MESV_MOD]])

p_me <- ggplot(df, aes(subtype, hME, fill = subtype)) +
  geom_violin(scale = "width", trim = TRUE, color = NA, alpha = 0.85) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = SUBTYPE_COL, guide = "none") +
  labs(x = NULL, y = paste0(MESV_MOD, " module eigengene (hME)"),
       title = "MES-V mesenchymal module activity") +
  theme_classic(base_size = 12)
ggsave(file.path(fig_dir, "S_green_module_hME_by_subtype.pdf"), p_me, width = 5, height = 4)

DefaultAssay(obj) <- "RNA"
p_plaur <- VlnPlot(
  obj,
  features = "PLAUR",
  group.by = "subtype_bio",
  cols = SUBTYPE_COL[SUBTYPE_LEVELS],
  pt.size = 0,
  layer = "data"
) +
  labs(x = NULL, title = "PLAUR expression") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "S_PLAUR_expr_by_subtype.pdf"), p_plaur, width = 5, height = 4)

write.csv(df, file.path(tab_dir, "stage8_green_hME_by_cell.csv"), row.names = FALSE)

plaur_expr <- FetchData(obj, "PLAUR", layer = "data")[, 1]
plaur_summary <- data.frame(
  subtype = SUBTYPE_LEVELS,
  mean_PLAUR = as.numeric(tapply(plaur_expr, obj$subtype_bio, mean)[SUBTYPE_LEVELS]),
  pct_detected = as.numeric(tapply(plaur_expr > 0, obj$subtype_bio, mean)[SUBTYPE_LEVELS]) * 100,
  mean_green_hME = as.numeric(tapply(df$hME, df$subtype, mean)[SUBTYPE_LEVELS])
)
write.csv(plaur_summary, file.path(tab_dir, "stage8_PLAUR_green_summary.csv"), row.names = FALSE)
print(plaur_summary)

cat("\nSTAGE 8 done. Plots and source data written; no object saved.\n")
