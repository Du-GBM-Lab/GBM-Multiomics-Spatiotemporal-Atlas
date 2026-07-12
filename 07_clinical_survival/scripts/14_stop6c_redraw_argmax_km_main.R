# STOP6c: redraw main dominant-program KM panel with publication-facing labels.
# Same data/statistics as STOP5_KM_argmax_CGGA_HGG_504; visual cleanup only.

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
LABELS <- c(NPC_P = "NPC-P", OPC_M = "OPC-M", MES_V = "MES-V", MES_I = "MES-I")
COLORS <- c(NPC_P = "#0072B5", OPC_M = "#E18727", MES_V = "#20854E", MES_I = "#BC3C29")

es <- as.data.frame(readRDS(file.path(root, "data", "processed", "ssgsea_scores_harmonized.rds")))
if (all(PROGS %in% rownames(es))) es <- as.data.frame(t(as.matrix(es)))
hgg <- readRDS(file.path(root, "data", "raw", "CGGA_old_manuscript", "hgg_data.rds"))
clin <- hgg$clinical

common_ids <- intersect(rownames(es), clin$CGGA_ID)
es <- es[common_ids, , drop = FALSE]
clin <- clin[match(common_ids, clin$CGGA_ID), , drop = FALSE]
stopifnot(all(rownames(es) == clin$CGGA_ID))

d <- data.frame(
  es[, PROGS, drop = FALSE],
  OS.time = suppressWarnings(as.numeric(clin$OS)),
  OS.event = suppressWarnings(as.numeric(clin$Censor..alive.0..dead.1.)),
  stringsAsFactors = FALSE
)
d <- d[!is.na(d$OS.time) & d$OS.time > 0 & !is.na(d$OS.event) & d$OS.event %in% c(0, 1), , drop = FALSE]
zall <- scale(as.matrix(d[, PROGS]))
d$dominant_program <- factor(PROGS[max.col(zall, ties.method = "first")], levels = PROGS)

fit <- survival::survfit(Surv(OS.time, OS.event) ~ dominant_program, data = d)
lr <- survival::survdiff(Surv(OS.time, OS.event) ~ dominant_program, data = d)
pval <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)

s <- summary(fit)
km <- data.frame(
  time_years = s$time / 365,
  survival = s$surv,
  strata = sub("^dominant_program=", "", as.character(s$strata)),
  stringsAsFactors = FALSE
)
km$program <- factor(km$strata, levels = PROGS)
counts <- table(d$dominant_program)
legend_labels <- paste0(LABELS[names(counts)], " (n=", as.integer(counts), ")")
names(legend_labels) <- names(counts)

p <- ggplot(km, aes(time_years, survival, color = program)) +
  geom_step(linewidth = 1.05) +
  scale_color_manual(values = COLORS, breaks = PROGS, labels = legend_labels, name = "Dominant program") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0.01, 0)) +
  scale_x_continuous(limits = c(0, 8.5), breaks = seq(0, 8, 2), expand = c(0.01, 0)) +
  annotate("text", x = 0.35, y = 0.08,
           label = paste0("log-rank p = ", format.pval(pval, digits = 2, eps = 1e-300)),
           hjust = 0, size = 3.4) +
  labs(
    x = "Time (years)",
    y = "Overall survival"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.line = element_line(color = "black", linewidth = 0.45),
    axis.ticks = element_line(color = "black", linewidth = 0.45),
    plot.margin = margin(6, 8, 6, 6)
  )

ggsave(file.path(fig_dir, "Fig_R4_argmax_KM_CGGA_HGG_504_clean.pdf"), p, width = 5.6, height = 4.0)

write.csv(data.frame(program = names(counts), label = LABELS[names(counts)], n = as.integer(counts)),
          file.path(tab_dir, "STOP6c_argmax_KM_CGGA_HGG_504_counts.csv"), row.names = FALSE)
write.csv(km, file.path(tab_dir, "STOP6c_argmax_KM_CGGA_HGG_504_curve_data.csv"), row.names = FALSE)
write.csv(data.frame(chisq = lr$chisq, df = length(lr$n) - 1, p = pval),
          file.path(tab_dir, "STOP6c_argmax_KM_CGGA_HGG_504_logrank.csv"), row.names = FALSE)

cat("\n########## STOP6c REPORT ##########\n")
cat("1. Clean KM written to figures/Fig_R4_argmax_KM_CGGA_HGG_504_clean.pdf.\n")
cat("2. Sample counts:\n")
print(data.frame(program = names(counts), label = LABELS[names(counts)], n = as.integer(counts)))
cat(sprintf("3. log-rank p = %.3g.\n", pval))
cat("=> Visual cleanup only; dominant-program grouping and statistics unchanged.\n")
