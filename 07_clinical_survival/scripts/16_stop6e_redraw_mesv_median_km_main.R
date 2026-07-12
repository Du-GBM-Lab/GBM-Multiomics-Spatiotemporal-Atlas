# STOP6e: redraw MES_V median-split KM for the main R4 figure.
# Same grouping/statistics as STOP5_KM_MESV_median_CGGA_HGG_504; visual cleanup only.

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
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
  MES_V = es$MES_V,
  OS.time = suppressWarnings(as.numeric(clin$OS)),
  OS.event = suppressWarnings(as.numeric(clin$Censor..alive.0..dead.1.)),
  stringsAsFactors = FALSE
)
d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS.event) & d$OS.event %in% c(0, 1), , drop = FALSE]
cut <- median(d$MES_V, na.rm = TRUE)
d$MESV_group <- factor(ifelse(d$MES_V > cut, "MES_V-high", "MES_V-low"),
                       levels = c("MES_V-low", "MES_V-high"))

fit <- survival::survfit(Surv(OS.time, OS.event) ~ MESV_group, data = d)
lr <- survival::survdiff(Surv(OS.time, OS.event) ~ MESV_group, data = d)
pval <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)

s <- summary(fit)
km <- data.frame(
  time_years = s$time / 365,
  survival = s$surv,
  strata = sub("^MESV_group=", "", as.character(s$strata)),
  stringsAsFactors = FALSE
)
km$group <- factor(km$strata, levels = c("MES_V-low", "MES_V-high"))
counts <- table(d$MESV_group)
legend_labels <- c(
  "MES_V-low" = paste0("MES-V low (n=", counts["MES_V-low"], ")"),
  "MES_V-high" = paste0("MES-V high (n=", counts["MES_V-high"], ")")
)

p <- ggplot(km, aes(time_years, survival, color = group)) +
  geom_step(linewidth = 1.05) +
  scale_color_manual(values = c("MES_V-low" = "#2C7BB6", "MES_V-high" = "#BC3C29"),
                     breaks = c("MES_V-low", "MES_V-high"),
                     labels = legend_labels,
                     name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0.01, 0)) +
  scale_x_continuous(limits = c(0, 8.5), breaks = seq(0, 8, 2), expand = c(0.01, 0)) +
  annotate("text", x = 0.35, y = 0.08,
           label = paste0("log-rank p = ", format.pval(pval, digits = 2, eps = 1e-300)),
           hjust = 0, size = 3.4) +
  labs(x = "Time (years)", y = "Overall survival") +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 9),
    axis.line = element_line(color = "black", linewidth = 0.45),
    axis.ticks = element_line(color = "black", linewidth = 0.45),
    plot.margin = margin(6, 8, 6, 6)
  )

ggsave(file.path(fig_dir, "Fig_R4_MESV_median_KM_CGGA_HGG_504_clean.pdf"), p, width = 5.2, height = 4.0)

write.csv(km, file.path(tab_dir, "STOP6e_MESV_median_KM_CGGA_HGG_504_curve_data.csv"), row.names = FALSE)
write.csv(data.frame(group = names(counts), n = as.integer(counts), cutoff = cut),
          file.path(tab_dir, "STOP6e_MESV_median_KM_CGGA_HGG_504_counts.csv"), row.names = FALSE)
write.csv(data.frame(chisq = lr$chisq, df = length(lr$n) - 1, p = pval, cutoff = cut),
          file.path(tab_dir, "STOP6e_MESV_median_KM_CGGA_HGG_504_logrank.csv"), row.names = FALSE)

cat("\n########## STOP6e REPORT ##########\n")
cat("1. Clean KM written to figures/Fig_R4_MESV_median_KM_CGGA_HGG_504_clean.pdf.\n")
cat("2. Group counts:\n")
print(data.frame(group = names(counts), n = as.integer(counts), cutoff = cut))
cat(sprintf("3. log-rank p = %.3g; median cutoff = %.4f.\n", pval, cut))
cat("=> Visual cleanup only; median split grouping and statistics unchanged.\n")
