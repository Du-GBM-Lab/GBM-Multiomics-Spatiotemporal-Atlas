# STOP5b CGGA mRNAseq_325 validation
# Same descriptive MES landscape口径 as discovery/TCGA:
# independent ssGSEA per cohort, HGG = WHO III/IV, no PLAUR/FOSL genes in signatures.

root <- getwd()
source(file.path(root, "scripts", "R4_helpers.R"))

raw_dir <- file.path(root, "data", "raw", "CGGA_325")
proc_dir <- file.path(root, "data", "processed")
tab_dir <- file.path(root, "tables")
fig_dir <- file.path(root, "figures")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

expr_file <- file.path(raw_dir, "CGGA.mRNAseq_325.RSEM-genes.20200506.txt")
clin_file <- file.path(raw_dir, "CGGA.mRNAseq_325_clinical.20200506.txt")
gene_sets <- readRDS(file.path(proc_dir, "sc_subtype_markers.rds"))

expr_df <- read.delim(expr_file, check.names = FALSE, stringsAsFactors = FALSE)
gene_col <- "Gene_Name"
stopifnot(gene_col %in% colnames(expr_df))
expr <- as.matrix(expr_df[, setdiff(colnames(expr_df), gene_col), drop = FALSE])
mode(expr) <- "numeric"
rownames(expr) <- expr_df[[gene_col]]

clin <- read.delim(clin_file, check.names = FALSE, stringsAsFactors = FALSE)
censor_col <- "Censor (alive=0; dead=1)"
stopifnot(all(c("CGGA_ID", "Grade", "PRS_type", "OS", censor_col) %in% colnames(clin)))

hgg_ids <- clin$CGGA_ID[clin$Grade %in% c("WHO III", "WHO IV")]
hgg_ids <- intersect(hgg_ids, colnames(expr))
expr_hgg <- expr[, hgg_ids, drop = FALSE]

surv <- data.frame(
  id = clin$CGGA_ID,
  OS.time = suppressWarnings(as.numeric(clin$OS)),
  OS.event = suppressWarnings(as.numeric(clin[[censor_col]])),
  PRS_type = clin$PRS_type,
  Grade = clin$Grade,
  IDH = clin$IDH_mutation_status,
  stringsAsFactors = FALSE
)
surv_hgg <- surv[surv$id %in% hgg_ids, , drop = FALSE]
surv_hgg <- surv_hgg[match(colnames(expr_hgg), surv_hgg$id), , drop = FALSE]
stopifnot(all(colnames(expr_hgg) == surv_hgg$id))

write.csv(data.frame(metric = c("genes", "patients", "OS_available", "events"),
                     value = c(nrow(expr_hgg), ncol(expr_hgg),
                               sum(!is.na(surv_hgg$OS.time) & surv_hgg$OS.time > 0),
                               sum(surv_hgg$OS.event == 1, na.rm = TRUE))),
          file.path(tab_dir, "STOP5b_CGGA325_HGG_expression_survival_audit.csv"),
          row.names = FALSE)

es_325 <- score_cohort(expr_hgg, gene_sets, tag = "CGGA325_HGG", out_dir = tab_dir)
saveRDS(es_325, file.path(proc_dir, "ssgsea_scores_CGGA325_HGG.rds"))
write.csv(data.frame(CGGA_ID = rownames(es_325), es_325, check.names = FALSE),
          file.path(tab_dir, "STOP5b_CGGA325_HGG_ssGSEA_scores.csv"), row.names = FALSE)

val_325_all <- mes_landscape(es_325, surv_hgg, "CGGA325_HGG", fig_dir)

primary_ids <- surv_hgg$id[grepl("Primary", surv_hgg$PRS_type, ignore.case = TRUE)]
val_325_primary <- mes_landscape(es_325[primary_ids, , drop = FALSE],
                                 surv_hgg[surv_hgg$id %in% primary_ids, , drop = FALSE],
                                 "CGGA325_HGG_primary", fig_dir)

val_325 <- rbind(val_325_all, val_325_primary)
write.csv(val_325, file.path(tab_dir, "STOP5b_CGGA325_validation_summary.csv"), row.names = FALSE)
saveRDS(val_325, file.path(proc_dir, "CGGA325_validation_summary.rds"))

cat("\n== val_325 ==\n")
print(val_325)

cat("\n########## STOP5b CGGA325 REPORT ##########\n")
cat("1. CGGA325 HGG = WHO III + WHO IV, expression/clinical IDs exactly matched.\n")
cat("2. score_cohort coverage printed above and saved to STOP5_CGGA325_HGG_score_coverage.csv.\n")
cat("3. val_325 printed above and saved to STOP5b_CGGA325_validation_summary.csv.\n")
