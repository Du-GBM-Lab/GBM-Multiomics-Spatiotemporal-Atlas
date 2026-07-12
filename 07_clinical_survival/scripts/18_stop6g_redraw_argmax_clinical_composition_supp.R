# STOP6g: redraw supplementary dominant-program clinical composition stacked bars.
# Shows IDH/grade/PRS composition across dominant program arms; descriptive only.

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
LABELS <- c(NPC_P = "NPC-P", OPC_M = "OPC-M", MES_V = "MES-V", MES_I = "MES-I")

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
  es[, PROGS, drop = FALSE],
  IDH = clin$IDH_mutation_status,
  Grade = clin$Grade,
  PRS = clin$PRS_type,
  stringsAsFactors = FALSE
)
d$dominant_program <- PROGS[max.col(scale(as.matrix(d[, PROGS])), ties.method = "first")]
d$dominant_program <- factor(d$dominant_program, levels = PROGS, labels = LABELS[PROGS])

make_long <- function(var, label, levels_keep, labels_keep = levels_keep) {
  dd <- d[!is.na(d[[var]]) & d[[var]] != "" & d[[var]] %in% levels_keep,
          c("CGGA_ID", "dominant_program", var), drop = FALSE]
  names(dd)[3] <- "category"
  dd$variable <- label
  dd$category <- factor(dd$category, levels = levels_keep, labels = labels_keep)
  dd
}

long <- rbind(
  make_long("IDH", "IDH status", c("Mutant", "Wildtype")),
  make_long("Grade", "WHO grade", c("WHO III", "WHO IV")),
  make_long("PRS", "Sample type", c("Primary", "Secondary", "Recurrent"))
)
long$variable <- factor(long$variable, levels = c("IDH status", "WHO grade", "Sample type"))

comp <- long %>%
  count(variable, dominant_program, category, name = "n") %>%
  group_by(variable, dominant_program) %>%
  mutate(total = sum(n), proportion = n / total) %>%
  ungroup()
write.csv(comp, file.path(tab_dir, "STOP6g_argmax_clinical_composition_source.csv"), row.names = FALSE)

colors <- c(
  Mutant = "#377EB8",
  Wildtype = "#E41A1C",
  "WHO III" = "#FDAE61",
  "WHO IV" = "#A50026",
  Primary = "#984EA3",
  Secondary = "grey55",
  Recurrent = "#FF7F00"
)

p <- ggplot(comp, aes(dominant_program, proportion, fill = category)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  facet_wrap(~ variable, nrow = 1) +
  scale_fill_manual(values = colors, name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = "Dominant program", y = "Proportion") +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 9, angle = 0),
    legend.position = "right",
    legend.text = element_text(size = 9),
    panel.spacing = grid::unit(5, "mm"),
    plot.margin = margin(6, 6, 6, 6)
  )

ggsave(file.path(fig_dir, "FigS_R4_argmax_clinical_composition_clean.pdf"), p, width = 8.2, height = 3.8)

cat("\n########## STOP6g REPORT ##########\n")
cat("1. Clean composition plot written to figures/FigS_R4_argmax_clinical_composition_clean.pdf.\n")
cat("2. Source written to tables/STOP6g_argmax_clinical_composition_source.csv.\n")
cat("3. Composition summary:\n")
print(comp)
cat("=> Supplementary only; shows clinical-gradient composition of dominant-program arms.\n")
