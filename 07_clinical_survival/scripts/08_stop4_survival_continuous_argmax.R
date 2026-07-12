# R4 08 / STOP4: survival analysis after STOP3 route change.
# Main logic:
# 1) Four program scores are continuous bulk programs, not discrete clusters.
# 2) Per-SD Cox estimates each program's prognostic direction.
# 3) Argmax dominant-program KM asks which program each tumor is most biased toward.
# No PLAUR/FOSL1/FOSL2 in R4.

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(dplyr)
})

root <- getwd()
raw_dir <- file.path(root, "data", "raw", "CGGA_old_manuscript")
proc_dir <- file.path(root, "data", "processed")
tab_dir <- file.path(root, "tables")
fig_dir <- file.path(root, "figures")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

PROGS <- c("NPC_P", "OPC_M", "MES_V", "MES_I")
FORBIDDEN <- c("PLAUR", "FOSL1", "FOSL2")

score_path <- file.path(proc_dir, "ssgsea_scores_harmonized.rds")
hgg_path <- file.path(raw_dir, "hgg_data.rds")
gene_set_path <- file.path(proc_dir, "sc_subtype_markers.rds")
stopifnot(file.exists(score_path), file.exists(hgg_path), file.exists(gene_set_path))
gene_sets <- readRDS(gene_set_path)
stopifnot(length(intersect(toupper(unlist(gene_sets)), FORBIDDEN)) == 0)

es <- as.matrix(readRDS(score_path))
if (all(PROGS %in% rownames(es))) {
  score <- as.data.frame(t(es[PROGS, , drop = FALSE]))
} else if (all(PROGS %in% colnames(es))) {
  score <- as.data.frame(es[, PROGS, drop = FALSE])
} else {
  stop("Cannot identify ssGSEA score orientation.")
}
score$CGGA_ID <- rownames(score)

hgg <- readRDS(hgg_path)
clin <- hgg$clinical
clin <- clin[match(score$CGGA_ID, as.character(clin$CGGA_ID)), , drop = FALSE]
stopifnot(all(score$CGGA_ID == as.character(clin$CGGA_ID)))

norm_string <- function(x) {
  trimws(as.character(x))
}
df <- cbind(score[, c("CGGA_ID", PROGS)], data.frame(
  OS.time = suppressWarnings(as.numeric(clin$OS)),
  OS.event = suppressWarnings(as.numeric(clin$Censor..alive.0..dead.1.)),
  age = suppressWarnings(as.numeric(clin$Age)),
  sex = norm_string(clin$Gender),
  idh = norm_string(clin$IDH_mutation_status),
  grade = norm_string(clin$Grade),
  prs = norm_string(clin$PRS_type),
  mgmt = norm_string(clin$MGMTp_methylation_status),
  codel = norm_string(clin$X1p19q_codeletion_status),
  histology = norm_string(clin$Histology),
  stringsAsFactors = FALSE
))

df <- df %>%
  filter(!is.na(OS.time), OS.time > 0, !is.na(OS.event), OS.event %in% c(0, 1))

df$sex <- factor(ifelse(grepl("^male$", df$sex, ignore.case = TRUE), "Male",
                        ifelse(grepl("^female$", df$sex, ignore.case = TRUE), "Female", NA)))
df$idh <- factor(ifelse(grepl("wild", df$idh, ignore.case = TRUE), "Wildtype",
                        ifelse(grepl("mut", df$idh, ignore.case = TRUE), "Mutant", NA)),
                 levels = c("Mutant", "Wildtype"))
df$grade <- factor(df$grade)
df$prs <- factor(ifelse(grepl("primary", df$prs, ignore.case = TRUE), "Primary",
                        ifelse(grepl("recurrent", df$prs, ignore.case = TRUE), "Recurrent", NA)),
                 levels = c("Primary", "Recurrent"))
df$mgmt <- factor(ifelse(grepl("^methylated$", df$mgmt, ignore.case = TRUE), "methylated",
                         ifelse(grepl("un-methylated|unmethylated", df$mgmt, ignore.case = TRUE),
                                "un-methylated", NA)))

sets <- list(
  all504 = df,
  primary = df %>% filter(prs == "Primary"),
  idhwt_gbm_primary = df %>% filter(prs == "Primary", idh == "Wildtype", grade == "WHO IV")
)

set_summary <- data.frame(
  set = names(sets),
  n = sapply(sets, nrow),
  events = sapply(sets, function(d) sum(d$OS.event == 1, na.rm = TRUE)),
  median_OS = sapply(sets, function(d) median(d$OS.time, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
write.csv(set_summary, file.path(tab_dir, "STOP4_analysis_set_summary.csv"), row.names = FALSE)
cat("== analysis set n / events ==\n")
print(set_summary)

cox_per_sd <- function(d, score_col, covars = character()) {
  keep <- c("OS.time", "OS.event", score_col, covars)
  dd <- d[, keep, drop = FALSE]
  dd <- dd[complete.cases(dd), , drop = FALSE]
  dd$z_score <- as.numeric(scale(dd[[score_col]]))
  rhs <- c("z_score", covars)
  f <- as.formula(paste("Surv(OS.time, OS.event) ~", paste(rhs, collapse = " + ")))
  m <- survival::coxph(f, data = dd)
  sm <- summary(m)
  zph <- tryCatch(survival::cox.zph(m)$table["z_score", "p"], error = function(e) NA_real_)
  data.frame(
    score = score_col,
    n = m$n,
    events = sum(dd$OS.event == 1),
    HR = unname(sm$coef["z_score", "exp(coef)"]),
    lo = unname(sm$conf.int["z_score", "lower .95"]),
    hi = unname(sm$conf.int["z_score", "upper .95"]),
    p = unname(sm$coef["z_score", "Pr(>|z|)"]),
    zph_p = zph,
    covars = paste(covars, collapse = "+"),
    stringsAsFactors = FALSE
  )
}

cox_rows <- data.frame()
cox_rows <- rbind(cox_rows, transform(cox_per_sd(sets$all504, "MES_V"),
                                      set = "all504", model = "MES_V_uni"))
cox_rows <- rbind(cox_rows, transform(cox_per_sd(sets$all504, "MES_V",
                                                 c("age", "sex", "idh", "grade", "mgmt")),
                                      set = "all504", model = "MES_V_adj_age_sex_idh_grade_mgmt"))
cox_rows <- rbind(cox_rows, transform(cox_per_sd(sets$idhwt_gbm_primary, "MES_V"),
                                      set = "idhwt_gbm_primary", model = "MES_V_uni"))
cox_rows <- rbind(cox_rows, transform(cox_per_sd(sets$idhwt_gbm_primary, "MES_V",
                                                 c("age", "sex", "mgmt")),
                                      set = "idhwt_gbm_primary", model = "MES_V_adj_age_sex_mgmt"))

forest_rows <- data.frame()
for (set_name in c("all504", "idhwt_gbm_primary")) {
  covs <- if (set_name == "all504") c("idh", "grade") else c("age", "sex")
  tmp <- do.call(rbind, lapply(PROGS, function(p) cox_per_sd(sets[[set_name]], p, covs)))
  tmp$set <- set_name
  tmp$model <- paste0("four_program_", paste(covs, collapse = "_"))
  forest_rows <- rbind(forest_rows, tmp)
}
cox_all <- rbind(cox_rows, forest_rows)
cox_all <- cox_all[, c("set", "model", "score", "n", "events", "HR", "lo", "hi", "p", "zph_p", "covars")]
write.csv(cox_all, file.path(tab_dir, "STOP4_perSD_Cox_results.csv"), row.names = FALSE)

cat("\n== MES_V per-SD Cox ==\n")
print(cox_rows[, c("set", "model", "n", "events", "HR", "lo", "hi", "p", "zph_p", "covars")])
cat("\n== four-program per-SD Cox forest rows ==\n")
print(forest_rows[, c("set", "score", "n", "events", "HR", "lo", "hi", "p", "zph_p", "covars")])

forest_plot <- forest_rows
forest_plot$score <- factor(forest_plot$score, levels = rev(PROGS))
forest_plot$set <- factor(forest_plot$set, levels = c("all504", "idhwt_gbm_primary"))
p_forest <- ggplot(forest_plot, aes(HR, score)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey45") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.18, orientation = "y", color = "grey25") +
  geom_point(size = 2.3, color = "#0072B5") +
  facet_wrap(~ set, scales = "free_x") +
  scale_x_log10() +
  labs(x = "HR per 1 SD (95% CI)", y = NULL) +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "STOP4_forest_4program_perSD.pdf"), p_forest, width = 8, height = 3.6)

make_km_data <- function(fit) {
  s <- summary(fit)
  strata <- if (is.null(s$strata)) rep("all", length(s$time)) else as.character(s$strata)
  data.frame(time = s$time, surv = s$surv, lower = s$lower, upper = s$upper, strata = strata)
}

plot_km <- function(d, group_col, tag, title) {
  f <- as.formula(paste("Surv(OS.time, OS.event) ~", group_col))
  fit <- survfit(f, data = d)
  lr <- survdiff(f, data = d)
  pval <- pchisq(lr$chisq, length(lr$n) - 1, lower.tail = FALSE)
  km <- make_km_data(fit)
  p <- ggplot(km, aes(time / 365, surv, color = strata)) +
    geom_step(linewidth = 0.8) +
    labs(x = "Time (years)", y = "Overall survival", color = NULL,
         title = title, subtitle = paste0("log-rank p = ", signif(pval, 3))) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 11))
  ggsave(file.path(fig_dir, paste0("STOP4_KM_", tag, ".pdf")), p, width = 6, height = 4.3)
  data.frame(tag = tag, chisq = lr$chisq, df = length(lr$n) - 1, p = pval)
}

km_rows <- data.frame()
for (set_name in c("primary", "idhwt_gbm_primary")) {
  d <- sets[[set_name]]
  d$MES_V_median_group <- factor(
    ifelse(d$MES_V > median(d$MES_V, na.rm = TRUE), "MES_V-high", "MES_V-low"),
    levels = c("MES_V-low", "MES_V-high")
  )
  tab <- table(d$MES_V_median_group)
  write.csv(data.frame(group = names(tab), n = as.integer(tab)),
            file.path(tab_dir, paste0("STOP4_KM_MESV_median_", set_name, "_group_counts.csv")),
            row.names = FALSE)
  km_rows <- rbind(km_rows, plot_km(d, "MES_V_median_group",
                                    paste0("MESV_median_", set_name),
                                    paste0("MES_V median split: ", set_name)))
}

assign_argmax <- function(d) {
  z <- scale(as.matrix(d[, PROGS]))
  factor(PROGS[max.col(z, ties.method = "first")], levels = PROGS)
}

argmax_rows <- data.frame()
for (set_name in c("all504", "primary")) {
  d <- sets[[set_name]]
  d$dominant_program <- assign_argmax(d)
  dom_tab <- table(d$dominant_program, useNA = "ifany")
  dom_idh <- table(d$dominant_program, d$idh, useNA = "ifany")
  dom_grade <- table(d$dominant_program, d$grade, useNA = "ifany")
  write.csv(data.frame(program = names(dom_tab), n = as.integer(dom_tab)),
            file.path(tab_dir, paste0("STOP4_argmax_", set_name, "_dominant_program_counts.csv")),
            row.names = FALSE)
  write.csv(as.data.frame(dom_idh),
            file.path(tab_dir, paste0("STOP4_argmax_", set_name, "_dominant_program_x_IDH.csv")),
            row.names = FALSE)
  write.csv(as.data.frame(dom_grade),
            file.path(tab_dir, paste0("STOP4_argmax_", set_name, "_dominant_program_x_grade.csv")),
            row.names = FALSE)
  cat("\n== argmax dominant program: ", set_name, " ==\n", sep = "")
  print(dom_tab)
  cat("dom x IDH:\n")
  print(dom_idh)
  cat("dom x grade:\n")
  print(dom_grade)
  km_rows <- rbind(km_rows, plot_km(d, "dominant_program",
                                    paste0("argmax_", set_name),
                                    paste0("Dominant program: ", set_name)))
  argmax_rows <- rbind(argmax_rows, data.frame(
    set = set_name,
    program = names(dom_tab),
    n = as.integer(dom_tab),
    stringsAsFactors = FALSE
  ))
}

write.csv(km_rows, file.path(tab_dir, "STOP4_KM_logrank_results.csv"), row.names = FALSE)
write.csv(argmax_rows, file.path(tab_dir, "STOP4_argmax_dominant_program_counts.csv"), row.names = FALSE)

cat("\n== KM log-rank results ==\n")
print(km_rows)

cat("\n\n########## STOP4 REPORT ##########\n")
cat("1. MES_V per-SD HR(CI,p): see STOP4_perSD_Cox_results.csv and printed MES_V rows; cox.zph p also reported.\n")
cat("2. all504 adjusted MES_V result indicates whether signal remains after age/sex/IDH/grade/MGMT.\n")
cat("3. idhwt_gbm_primary MES_V result reports significant/trend/negative within primary IDH-wt WHO IV GBM.\n")
cat("4. forest rows report all four programs in all504 and idhwt_gbm_primary.\n")
cat("5. argmax four-arm KM includes dominant program counts plus dom x IDH/grade for all504 and primary.\n")
cat("6. MES_V median KM generated for primary and idhwt_gbm_primary.\n")
cat("=> STOP4 complete. STOP5 remains optimal-cut sensitivity + ECMcore robustness if needed.\n")
