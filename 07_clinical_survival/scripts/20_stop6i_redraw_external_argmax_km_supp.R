# STOP6i: redraw external dominant-program KM panels with publication-facing labels.
# Visual cleanup only; dominant-program assignment and log-rank statistics follow STOP5b.

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
})

root <- getwd()
fig_dir <- file.path(root, "figures")
tab_dir <- file.path(root, "tables")
proc_dir <- file.path(root, "data", "processed")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

PROGS <- c("NPC_P", "OPC_M", "MES_V", "MES_I")
LABELS <- c(NPC_P = "NPC-P", OPC_M = "OPC-M", MES_V = "MES-V", MES_I = "MES-I")
COLORS <- c(NPC_P = "#0072B5", OPC_M = "#E18727", MES_V = "#20854E", MES_I = "#BC3C29")

clean_tcga_clin <- function(project) {
  cache <- file.path(proc_dir, paste0(project, "_clinical_survival.rds"))
  if (file.exists(cache)) return(readRDS(cache))
  x <- TCGAbiolinks::GDCquery_clinic(project, "clinical")
  for (nm in c("submitter_id", "vital_status", "days_to_death", "days_to_last_follow_up")) {
    if (!nm %in% colnames(x)) x[[nm]] <- NA
  }
  out <- data.frame(
    id = x$submitter_id,
    OS.event = ifelse(x$vital_status == "Dead", 1L, 0L),
    OS.time = ifelse(x$vital_status == "Dead",
                     suppressWarnings(as.numeric(x$days_to_death)),
                     suppressWarnings(as.numeric(x$days_to_last_follow_up))),
    stringsAsFactors = FALSE
  )
  saveRDS(out, cache)
  out
}

build_km <- function(es, surv, tag, display_name, x_max = NULL) {
  es <- as.data.frame(es)
  if (all(PROGS %in% rownames(es))) es <- as.data.frame(t(as.matrix(es)))
  stopifnot(all(PROGS %in% colnames(es)))

  surv <- surv[surv$id %in% rownames(es), , drop = FALSE]
  d <- data.frame(
    es[surv$id, PROGS, drop = FALSE],
    OS.time = suppressWarnings(as.numeric(surv$OS.time)),
    OS.event = suppressWarnings(as.numeric(surv$OS.event)),
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
  count_df <- data.frame(program = names(counts), label = LABELS[names(counts)], n = as.integer(counts))
  legend_labels <- paste0(LABELS[names(counts)], " (n=", as.integer(counts), ")")
  names(legend_labels) <- names(counts)

  if (is.null(x_max)) {
    x_max <- ceiling(max(km$time_years, na.rm = TRUE))
    x_max <- max(6, min(x_max, 12))
  }

  p <- ggplot(km, aes(time_years, survival, color = program)) +
    geom_step(linewidth = 1.0) +
    scale_color_manual(values = COLORS, breaks = PROGS, labels = legend_labels, name = "Dominant program") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0.01, 0)) +
    scale_x_continuous(limits = c(0, x_max), breaks = seq(0, x_max, by = 2), expand = c(0.01, 0)) +
    annotate("text", x = 0.35, y = 0.08,
             label = paste0("log-rank p = ", format.pval(pval, digits = 2, eps = 1e-300)),
             hjust = 0, size = 3.3) +
    labs(x = "Time (years)", y = "Overall survival", title = display_name) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 11, face = "bold"),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      axis.line = element_line(color = "black", linewidth = 0.45),
      axis.ticks = element_line(color = "black", linewidth = 0.45),
      plot.margin = margin(6, 8, 6, 6)
    )

  ggsave(file.path(fig_dir, paste0("FigS_R4_argmax_KM_", tag, "_clean.pdf")), p, width = 5.8, height = 4.0)
  write.csv(km, file.path(tab_dir, paste0("STOP6i_argmax_KM_", tag, "_curve_data.csv")), row.names = FALSE)
  write.csv(count_df, file.path(tab_dir, paste0("STOP6i_argmax_KM_", tag, "_counts.csv")), row.names = FALSE)
  write.csv(data.frame(tag = tag, chisq = lr$chisq, df = length(lr$n) - 1, p = pval),
            file.path(tab_dir, paste0("STOP6i_argmax_KM_", tag, "_logrank.csv")), row.names = FALSE)

  data.frame(tag = tag, display_name = display_name, n = nrow(d), events = sum(d$OS.event == 1),
             logrank_p = pval, stringsAsFactors = FALSE)
}

# TCGA HGG validation.
suppressPackageStartupMessages(library(TCGAbiolinks))
es_tcga <- readRDS(file.path(proc_dir, "ssgsea_scores_TCGA_HGG.rds"))
clin_t <- rbind(clean_tcga_clin("TCGA-GBM"), clean_tcga_clin("TCGA-LGG"))
tcga_res <- build_km(es_tcga, clin_t, "TCGA_HGG", "TCGA HGG validation", x_max = 12)

# CGGA-325 validation and primary-only sensitivity.
es_325 <- readRDS(file.path(proc_dir, "ssgsea_scores_CGGA325_HGG.rds"))
clin_file <- file.path(root, "data", "raw", "CGGA_325", "CGGA.mRNAseq_325_clinical.20200506.txt")
clin <- read.delim(clin_file, check.names = FALSE, stringsAsFactors = FALSE)
censor_col <- "Censor (alive=0; dead=1)"
surv_325 <- data.frame(
  id = clin$CGGA_ID,
  OS.time = suppressWarnings(as.numeric(clin$OS)),
  OS.event = suppressWarnings(as.numeric(clin[[censor_col]])),
  PRS_type = clin$PRS_type,
  stringsAsFactors = FALSE
)
surv_325 <- surv_325[surv_325$id %in% rownames(es_325), , drop = FALSE]
cgga325_res <- build_km(es_325, surv_325, "CGGA325_HGG", "CGGA-325 HGG validation", x_max = 10)

primary_ids <- surv_325$id[grepl("Primary", surv_325$PRS_type, ignore.case = TRUE)]
cgga325_primary_res <- build_km(es_325[primary_ids, , drop = FALSE],
                                surv_325[surv_325$id %in% primary_ids, , drop = FALSE],
                                "CGGA325_primary", "CGGA-325 primary sensitivity", x_max = 10)

summary_df <- rbind(tcga_res, cgga325_res, cgga325_primary_res)
write.csv(summary_df, file.path(tab_dir, "STOP6i_external_argmax_KM_summary.csv"), row.names = FALSE)

cat("\n########## STOP6i REPORT ##########\n")
print(summary_df)
cat("Clean external dominant-program KM panels written to figures/FigS_R4_argmax_KM_*_clean.pdf.\n")
cat("Visual cleanup only; dominant-program grouping and log-rank tests follow STOP5b.\n")
