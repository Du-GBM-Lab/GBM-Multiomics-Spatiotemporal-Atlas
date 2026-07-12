# R4 09 / STOP5a: discovery cohort finalization for CGGA HGG.
# Descriptive survival association only. No PLAUR/FOSL1/FOSL2.
# Includes ECMcore robustness and maxstat optimal-cut sensitivity.

suppressPackageStartupMessages({
  library(survival)
  library(maxstat)
})

root <- getwd()
source(file.path(root, "scripts", "R4_helpers.R"))

raw_dir <- file.path(root, "data", "raw", "CGGA_old_manuscript")
proc_dir <- file.path(root, "data", "processed")
tab_dir <- file.path(root, "tables")
fig_dir <- file.path(root, "figures")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

FORBIDDEN <- c("PLAUR", "FOSL1", "FOSL2")
score_path <- file.path(proc_dir, "ssgsea_scores_harmonized.rds")
hgg_path <- file.path(raw_dir, "hgg_data.rds")
gene_set_path <- file.path(proc_dir, "sc_subtype_markers.rds")
stopifnot(file.exists(score_path), file.exists(hgg_path), file.exists(gene_set_path))

gene_sets <- readRDS(gene_set_path)
stopifnot(length(intersect(toupper(unlist(gene_sets)), FORBIDDEN)) == 0)

es0 <- as.matrix(readRDS(score_path))
if (all(PROGS %in% rownames(es0))) {
  es <- t(es0[PROGS, , drop = FALSE])
} else if (all(PROGS %in% colnames(es0))) {
  es <- es0[, PROGS, drop = FALSE]
} else {
  stop("Cannot identify score orientation.")
}

hgg <- readRDS(hgg_path)
clin <- hgg$clinical
clin <- clin[match(rownames(es), as.character(clin$CGGA_ID)), , drop = FALSE]
stopifnot(all(rownames(es) == as.character(clin$CGGA_ID)))

surv <- data.frame(
  id = as.character(clin$CGGA_ID),
  OS.time = suppressWarnings(as.numeric(clin$OS)),
  OS.event = suppressWarnings(as.numeric(clin$Censor..alive.0..dead.1.)),
  prs = as.character(clin$PRS_type),
  stringsAsFactors = FALSE
)
primary_ids <- surv$id[grepl("primary", surv$prs, ignore.case = TRUE)]

cat("== STOP5a discovery landscape ==\n")
disc <- rbind(
  mes_landscape(es, surv, "CGGA_HGG_504", fig_dir),
  mes_landscape(es[primary_ids, , drop = FALSE],
                surv[surv$id %in% primary_ids, , drop = FALSE],
                "CGGA_HGG_primary", fig_dir)
)

write.csv(disc, file.path(tab_dir, "STOP5a_discovery_summary.csv"), row.names = FALSE)
saveRDS(disc, file.path(proc_dir, "discovery_summary.rds"))

cat("\n== ECMcore robustness ==\n")
expr_df <- hgg$expression
stopifnot("Gene_Name" %in% colnames(expr_df))
expr <- as.matrix(expr_df[, setdiff(colnames(expr_df), "Gene_Name")])
mode(expr) <- "numeric"
rownames(expr) <- as.character(expr_df$Gene_Name)

ECM_anchor <- c(
  "ACTA2", "TAGLN", "FN1", "COL1A2", "COL4A1", "COL4A2", "COL8A1",
  "CCN2", "CTGF", "FSTL1", "IGFBP7", "LUM", "MMP14", "MYL9",
  "PCOLCE", "SPARCL1", "TFPI", "TGFB1I1", "TGM2", "TPM1", "TPM2",
  "ANXA2", "VIM", "COL6A1", "COL1A1", "COL3A1"
)
MESV_current <- to_current(gene_sets$MES_V)
anchor_current <- to_current(ECM_anchor)
ECM_core <- intersect(MESV_current[!is.na(MESV_current)], anchor_current[!is.na(anchor_current)])
ECM_core <- unique(ECM_core)
cat("ECMcore n=", length(ECM_core), "\n", sep = "")
print(sort(ECM_core))
write.csv(data.frame(gene = sort(ECM_core)),
          file.path(tab_dir, "STOP5a_MESV_ECMcore_genes.csv"), row.names = FALSE)

es_ecm <- score_cohort(expr, list(MES_V_ECMcore = ECM_core), tag = "CGGA_HGG_ECMcore", out_dir = tab_dir)
surv_ecm <- surv
d_ecm <- data.frame(id = rownames(es_ecm), MES_V_ECMcore = es_ecm[, 1])
d_ecm <- merge(d_ecm, surv_ecm[, c("id", "OS.time", "OS.event")], by = "id", all.x = TRUE)
d_ecm <- d_ecm[!is.na(d_ecm$OS.time) & d_ecm$OS.time > 0 & !is.na(d_ecm$OS.event), ]
d_ecm$z <- as.numeric(scale(d_ecm$MES_V_ECMcore))
m_ecm <- survival::coxph(Surv(OS.time, OS.event) ~ z, data = d_ecm)
s_ecm <- summary(m_ecm)
ecm_res <- data.frame(
  cohort = "CGGA_HGG_504",
  score = "MES_V_ECMcore",
  n = m_ecm$n,
  events = sum(d_ecm$OS.event == 1),
  HR = s_ecm$coef["z", "exp(coef)"],
  lo = s_ecm$conf.int["z", "lower .95"],
  hi = s_ecm$conf.int["z", "upper .95"],
  p = s_ecm$coef["z", "Pr(>|z|)"],
  zph_p = tryCatch(survival::cox.zph(m_ecm)$table["z", "p"], error = function(e) NA_real_)
)
write.csv(ecm_res, file.path(tab_dir, "STOP5a_ECMcore_Cox.csv"), row.names = FALSE)
cat(sprintf("[ECMcore all504] perSD HR=%.3f (%.3f-%.3f) p=%.2e\n",
            ecm_res$HR, ecm_res$lo, ecm_res$hi, ecm_res$p))

cat("\n== optimal-cut sensitivity with maxstat adjusted p ==\n")
d_opt <- data.frame(id = rownames(es), MES_V = es[, "MES_V"])
d_opt <- merge(d_opt, surv[, c("id", "OS.time", "OS.event")], by = "id", all.x = TRUE)
d_opt <- d_opt[!is.na(d_opt$OS.time) & d_opt$OS.time > 0 & !is.na(d_opt$OS.event), ]
set.seed(1)
mx <- maxstat::maxstat.test(
  Surv(OS.time, OS.event) ~ MES_V,
  data = d_opt,
  smethod = "LogRank",
  pmethod = "condMC",
  B = 9999
)
opt_res <- data.frame(
  cohort = "CGGA_HGG_504",
  score = "MES_V",
  cut = as.numeric(mx$estimate),
  maxstat = as.numeric(mx$statistic),
  adjusted_p = as.numeric(mx$p.value)
)
write.csv(opt_res, file.path(tab_dir, "STOP5a_MESV_optimal_cut_maxstat.csv"), row.names = FALSE)
cat(sprintf("[optimal-cut all504] cut=%.3f, maxstat adjusted p=%.4g\n",
            opt_res$cut, opt_res$adjusted_p))

cat("\n\n########## STOP5a REPORT ##########\n")
cat("1. CGGA 504 and primary MES_V perSD HR/p, argmax p, median-KM p are in STOP5a_discovery_summary.csv.\n")
cat("2. ECMcore HR/p is in STOP5a_ECMcore_Cox.csv; compare direction/scale with full MES_V.\n")
cat("3. optimal-cut maxstat adjusted p is in STOP5a_MESV_optimal_cut_maxstat.csv.\n")
cat("=> STOP5a complete. Discovery set can be judged before external validation.\n")
