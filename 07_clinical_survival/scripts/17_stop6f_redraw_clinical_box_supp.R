# STOP6f: redraw supplementary MES_V by clinical annotations boxplots.
# Visual cleanup only; same CGGA HGG scores/clinical annotations as STOP6a.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

root <- getwd()
fig_dir <- file.path(root, "figures")
tab_dir <- file.path(root, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

PROGS <- c("NPC_P", "OPC_M", "MES_V", "MES_I")

es <- as.data.frame(readRDS(file.path(root, "data", "processed", "ssgsea_scores_harmonized.rds")))
if (all(PROGS %in% rownames(es))) es <- as.data.frame(t(as.matrix(es)))
hgg <- readRDS(file.path(root, "data", "raw", "CGGA_old_manuscript", "hgg_data.rds"))
clin <- hgg$clinical

common_ids <- intersect(rownames(es), clin$CGGA_ID)
es <- es[common_ids, , drop = FALSE]
clin <- clin[match(common_ids, clin$CGGA_ID), , drop = FALSE]
stopifnot(all(rownames(es) == clin$CGGA_ID))

d <- data.frame(
  CGGA_ID = rownames(es),
  MES_V = es$MES_V,
  IDH = clin$IDH_mutation_status,
  Grade = clin$Grade,
  codel = clin$X1p19q_codeletion_status,
  PRS = clin$PRS_type,
  stringsAsFactors = FALSE
)

make_long <- function(var, label, order_levels) {
  dd <- d[!is.na(d[[var]]) & d[[var]] != "", c("CGGA_ID", "MES_V", var), drop = FALSE]
  names(dd)[3] <- "group"
  dd$variable <- label
  dd$group <- factor(dd$group, levels = order_levels)
  dd
}

long <- rbind(
  make_long("IDH", "IDH status", c("Mutant", "Wildtype")),
  make_long("Grade", "WHO grade", c("WHO III", "WHO IV")),
  make_long("codel", "1p/19q status", c("Codel", "Non-codel")),
  make_long("PRS", "Sample type", c("Primary", "Recurrent"))
)
long <- long[!is.na(long$group), , drop = FALSE]
long$variable <- factor(long$variable, levels = c("IDH status", "WHO grade", "1p/19q status", "Sample type"))

test_p <- function(y, g) {
  g <- factor(g)
  if (nlevels(g) == 2) wilcox.test(y ~ g)$p.value else kruskal.test(y, g)$p.value
}
stats <- long %>%
  group_by(variable) %>%
  summarise(p = test_p(MES_V, group), .groups = "drop") %>%
  mutate(p_label = ifelse(p < 1e-4, "p < 1e-4", paste0("p = ", signif(p, 2))))

source_df <- long %>% left_join(stats, by = "variable")
write.csv(source_df, file.path(tab_dir, "STOP6f_MESV_by_clinical_box_source.csv"), row.names = FALSE)
write.csv(stats, file.path(tab_dir, "STOP6f_MESV_by_clinical_box_stats.csv"), row.names = FALSE)

p <- ggplot(source_df, aes(group, MES_V)) +
  geom_boxplot(outlier.size = 0.4, width = 0.58, fill = "grey92", color = "grey25", linewidth = 0.55) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 0.35, alpha = 0.35, color = "grey35") +
  facet_wrap(~ variable, scales = "free_x", nrow = 2) +
  geom_text(data = stats, aes(x = 1.5, y = 0.77, label = p_label), inherit.aes = FALSE, size = 3.1) +
  scale_y_continuous(limits = c(0.32, 0.78), breaks = seq(0.35, 0.75, 0.1), expand = c(0.01, 0)) +
  labs(x = NULL, y = "MES-V ssGSEA score") +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 9),
    panel.spacing = grid::unit(5, "mm"),
    plot.margin = margin(6, 6, 6, 6)
  )

ggsave(file.path(fig_dir, "FigS_R4_MESV_by_clinical_box_clean.pdf"), p, width = 7.2, height = 5.2)

cat("\n########## STOP6f REPORT ##########\n")
cat("1. Clean clinical boxplot written to figures/FigS_R4_MESV_by_clinical_box_clean.pdf.\n")
cat("2. Source/stat tables written to STOP6f_MESV_by_clinical_box_*.csv.\n")
print(stats)
cat("=> Supplementary figure only; shows clinical gradient, not independence or causality.\n")
