# STOP6a: clinical association supplementary plots for R4.
# Descriptive only: associated/differ across; no causality, no PLAUR/FOSL.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
})

root <- getwd()
tab_dir <- file.path(root, "tables")
fig_dir <- file.path(root, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

PROGS <- c("NPC_P", "OPC_M", "MES_V", "MES_I")

es <- as.data.frame(readRDS(file.path(root, "data", "processed", "ssgsea_scores_harmonized.rds")))
if (all(PROGS %in% rownames(es))) es <- as.data.frame(t(as.matrix(es)))
hgg <- readRDS(file.path(root, "data", "raw", "CGGA_old_manuscript", "hgg_data.rds"))
common_ids <- intersect(rownames(es), hgg$clinical$CGGA_ID)
es <- es[common_ids, , drop = FALSE]
clin <- hgg$clinical[match(common_ids, hgg$clinical$CGGA_ID), , drop = FALSE]
stopifnot(all(rownames(es) == clin$CGGA_ID))

d <- data.frame(
  es,
  IDH = clin$IDH_mutation_status,
  Grade = clin$Grade,
  codel = clin$X1p19q_codeletion_status,
  PRS = clin$PRS_type,
  age = suppressWarnings(as.numeric(clin$Age)),
  sex = clin$Gender,
  stringsAsFactors = FALSE
)
d$dom <- factor(PROGS[max.col(scale(as.matrix(d[, PROGS])), ties.method = "first")], levels = PROGS)

test_p <- function(y, g) {
  ok <- !is.na(y) & !is.na(g) & g != ""
  y <- y[ok]
  g <- factor(g[ok])
  if (nlevels(g) < 2) return(NA_real_)
  if (nlevels(g) == 2) wilcox.test(y ~ g)$p.value else kruskal.test(y, g)$p.value
}

box_one <- function(var) {
  dd <- d[!is.na(d[[var]]) & d[[var]] != "", , drop = FALSE]
  g <- factor(dd[[var]])
  p <- test_p(dd$MES_V, dd[[var]])
  ggplot(dd, aes(g, MES_V)) +
    geom_boxplot(outlier.size = 0.5, fill = "grey90", color = "grey25") +
    labs(x = var, y = "MES_V ssGSEA",
         subtitle = sprintf("p = %.2g (%s)", p, ifelse(nlevels(g) == 2, "Wilcoxon", "Kruskal-Wallis"))) +
    theme_bw(base_size = 11)
}

box_plots <- lapply(c("IDH", "Grade", "codel", "PRS"), box_one)
ggsave(file.path(fig_dir, "FigS_MESV_by_clinical_box.pdf"),
       gridExtra::marrangeGrob(box_plots, nrow = 2, ncol = 2),
       width = 9, height = 8)

stack_one <- function(var) {
  dd <- d[!is.na(d[[var]]) & d[[var]] != "", , drop = FALSE]
  ggplot(dd, aes(dom, fill = factor(.data[[var]]))) +
    geom_bar(position = "fill", color = "white", linewidth = 0.1) +
    labs(x = "dominant program", y = "proportion", fill = var) +
    theme_bw(base_size = 11)
}
ggsave(file.path(fig_dir, "FigS_argmax_clinical_composition.pdf"),
       gridExtra::marrangeGrob(lapply(c("IDH", "Grade", "PRS"), stack_one), nrow = 1, ncol = 3),
       width = 12, height = 4)

rho <- cor(d$age, d$MES_V, use = "complete.obs", method = "spearman")
p_age <- suppressWarnings(cor.test(d$age, d$MES_V, method = "spearman", exact = FALSE)$p.value)
da <- d[is.finite(d$age) & is.finite(d$MES_V), , drop = FALSE]
ggsave(file.path(fig_dir, "FigS_age_vs_MESV.pdf"),
       ggplot(da, aes(age, MES_V)) +
         geom_point(size = 0.6, alpha = 0.5) +
         geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.6) +
         labs(x = "Age", y = "MES_V ssGSEA", subtitle = sprintf("Spearman rho = %.2f, p = %.2g", rho, p_age)) +
         theme_bw(base_size = 11),
       width = 5, height = 4)

ds <- d[!is.na(d$sex) & d$sex != "", , drop = FALSE]
p_sex <- test_p(ds$MES_V, ds$sex)
ggsave(file.path(fig_dir, "FigS_sex_vs_MESV.pdf"),
       ggplot(ds, aes(sex, MES_V)) +
         geom_boxplot(outlier.size = 0.5, fill = "grey90", color = "grey25") +
         labs(x = "Sex", y = "MES_V ssGSEA", subtitle = sprintf("p = %.2g (Wilcoxon)", p_sex)) +
         theme_bw(base_size = 11),
       width = 4, height = 4)

stat <- data.frame(
  test = c("MES_V_by_IDH", "MES_V_by_Grade", "MES_V_by_1p19q", "MES_V_by_PRS", "age_by_MES_V_spearman_rho", "age_by_MES_V_spearman_p", "MES_V_by_sex"),
  value = c(test_p(d$MES_V, d$IDH),
            test_p(d$MES_V, d$Grade),
            test_p(d$MES_V, d$codel),
            test_p(d$MES_V, d$PRS),
            rho, p_age, p_sex)
)
write.csv(stat, file.path(tab_dir, "STOP6a_clinical_association_stats.csv"), row.names = FALSE)

cat("\n########## STOP6a REPORT ##########\n")
cat("1. MES_V by IDH/grade/1p19q/PRS p-values:\n")
print(stat[stat$test %in% c("MES_V_by_IDH", "MES_V_by_Grade", "MES_V_by_1p19q", "MES_V_by_PRS"), ])
cat("2. argmax x clinical composition source is FigS_argmax_clinical_composition.pdf.\n")
cat(sprintf("3. age rho = %.3f, age p = %.3g; sex p = %.3g.\n", rho, p_age, p_sex))
cat("=> Supplementary only; do not use as causal or independent evidence.\n")
