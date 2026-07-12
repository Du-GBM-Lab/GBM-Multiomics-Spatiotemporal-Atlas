# R4 10 / STOP5b: external validation in TCGA-GBMLGG HGG.
# Main validation: HGG landscape = TCGA-GBM G4 + TCGA-LGG G3.
# Uses one primary sample per patient before 12-digit barcode truncation.
# CGGA-325 validation remains a later placeholder.

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})
options(timeout = 600)

root <- getwd()
source(file.path(root, "scripts", "R4_helpers.R"))

raw_dir <- file.path(root, "data", "raw", "TCGA")
proc_dir <- file.path(root, "data", "processed")
tab_dir <- file.path(root, "tables")
fig_dir <- file.path(root, "figures")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

gene_sets <- readRDS(file.path(proc_dir, "sc_subtype_markers.rds"))

dedup_primary <- function(m) {
  st <- suppressWarnings(as.integer(substr(colnames(m), 14, 15)))
  st[is.na(st)] <- 99
  ord <- order(st)
  m <- m[, ord, drop = FALSE]
  pt <- substr(colnames(m), 1, 12)
  keep <- !duplicated(pt)
  m <- m[, keep, drop = FALSE]
  colnames(m) <- pt[keep]
  m
}

get_tcga_expr <- function(project) {
  cache_file <- file.path(proc_dir, paste0(project, "_STAR_TPM_primary_dedup_symbol.rds"))
  if (file.exists(cache_file)) {
    cat("Using cached expression:", cache_file, "\n")
    return(readRDS(cache_file))
  }

  cat("Querying/downloading:", project, "\n")
  q <- GDCquery(
    project = project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts",
    sample.type = "Primary Tumor"
  )
  GDCdownload(q, directory = raw_dir, files.per.chunk = 25)
  se <- GDCprepare(q, directory = raw_dir)

  assay_names <- names(SummarizedExperiment::assays(se))
  cat(project, "assays:", paste(assay_names, collapse = ", "), "\n")
  assay_name <- if ("tpm_unstrand" %in% assay_names) {
    "tpm_unstrand"
  } else if ("fpkm_unstrand" %in% assay_names) {
    "fpkm_unstrand"
  } else if ("unstranded" %in% assay_names) {
    "unstranded"
  } else {
    assay_names[1]
  }

  m <- SummarizedExperiment::assay(se, assay_name)
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  gene_col <- if ("gene_name" %in% colnames(rd)) "gene_name" else if ("external_gene_name" %in% colnames(rd)) "external_gene_name" else NA
  if (is.na(gene_col)) stop("Cannot find gene_name in rowData.")
  rownames(m) <- as.character(rd[[gene_col]])
  m <- m[!is.na(rownames(m)) & rownames(m) != "", , drop = FALSE]
  m <- collapse_duplicate_rows(m)
  m <- dedup_primary(m)
  saveRDS(m, cache_file)
  write.csv(data.frame(project = project, assay = assay_name, genes = nrow(m), patients = ncol(m)),
            file.path(tab_dir, paste0("STOP5b_", project, "_expression_audit.csv")),
            row.names = FALSE)
  m
}

sub_gbm <- TCGAquery_subtype("gbm")
sub_lgg <- TCGAquery_subtype("lgg")
write.csv(sub_gbm, file.path(tab_dir, "STOP5b_TCGAquery_subtype_gbm.csv"), row.names = FALSE)
write.csv(sub_lgg, file.path(tab_dir, "STOP5b_TCGAquery_subtype_lgg.csv"), row.names = FALSE)

sub <- rbind(
  data.frame(project = "TCGA-GBM", patient = sub_gbm$patient, Grade = sub_gbm$Grade,
             IDH = sub_gbm$IDH.status, codel = sub_gbm$X1p.19q.codeletion,
             stringsAsFactors = FALSE),
  data.frame(project = "TCGA-LGG", patient = sub_lgg$patient, Grade = sub_lgg$Grade,
             IDH = sub_lgg$IDH.status, codel = sub_lgg$X1p.19q.codeletion,
             stringsAsFactors = FALSE)
)
sub$HGG_landscape <- sub$Grade %in% c("G3", "G4")
sub$WHO2021_strict_GBM <- sub$Grade == "G4" & sub$IDH == "WT"
write.csv(sub, file.path(tab_dir, "STOP5b_TCGA_subtype_HGG_filter_annotation.csv"), row.names = FALSE)

cat("== TCGA subtype HGG counts ==\n")
print(table(sub$project, sub$Grade, useNA = "ifany"))
print(table(sub$project, sub$IDH, useNA = "ifany"))
cat("HGG landscape patients:", sum(sub$HGG_landscape, na.rm = TRUE), "\n")
cat("WHO2021 strict GBM patients:", sum(sub$WHO2021_strict_GBM, na.rm = TRUE), "\n")

expr_gbm <- get_tcga_expr("TCGA-GBM")
expr_lgg <- get_tcga_expr("TCGA-LGG")
common_genes <- intersect(rownames(expr_gbm), rownames(expr_lgg))
expr_hgg <- cbind(expr_gbm[common_genes, , drop = FALSE],
                  expr_lgg[common_genes, , drop = FALSE])
hgg_pt <- sub$patient[sub$HGG_landscape]
expr_hgg <- expr_hgg[, colnames(expr_hgg) %in% hgg_pt, drop = FALSE]

cat("TCGA HGG expression genes:", nrow(expr_hgg), "patients:", ncol(expr_hgg), "\n")
write.csv(data.frame(metric = c("genes", "patients"),
                     value = c(nrow(expr_hgg), ncol(expr_hgg))),
          file.path(tab_dir, "STOP5b_TCGA_HGG_expression_after_filter.csv"),
          row.names = FALSE)

clean_tcga_clin <- function(project) {
  x <- NULL
  for (i in 1:3) {
    x <- tryCatch(GDCquery_clinic(project, "clinical"), error = function(e) e)
    if (!inherits(x, "error")) break
    if (i == 3) stop(x)
    Sys.sleep(10)
  }
  needed <- c("submitter_id", "vital_status", "days_to_death", "days_to_last_follow_up")
  for (nm in setdiff(needed, colnames(x))) x[[nm]] <- NA
  data.frame(
    id = x$submitter_id,
    OS.event = ifelse(x$vital_status == "Dead", 1L, 0L),
    OS.time = ifelse(x$vital_status == "Dead",
                     suppressWarnings(as.numeric(x$days_to_death)),
                     suppressWarnings(as.numeric(x$days_to_last_follow_up))),
    stringsAsFactors = FALSE
  )
}
clin_t <- rbind(clean_tcga_clin("TCGA-GBM"), clean_tcga_clin("TCGA-LGG"))
surv_t <- clin_t[clin_t$id %in% colnames(expr_hgg), c("id", "OS.time", "OS.event")]
surv_t <- surv_t[!duplicated(surv_t$id), , drop = FALSE]

es_tcga <- score_cohort(expr_hgg, gene_sets, tag = "TCGA_HGG", out_dir = tab_dir)
saveRDS(es_tcga, file.path(proc_dir, "ssgsea_scores_TCGA_HGG.rds"))
write.csv(data.frame(patient = rownames(es_tcga), es_tcga, check.names = FALSE),
          file.path(tab_dir, "STOP5b_TCGA_HGG_ssGSEA_scores.csv"), row.names = FALSE)

val_tcga <- mes_landscape(es_tcga, surv_t, "TCGA_HGG", fig_dir)
write.csv(val_tcga, file.path(tab_dir, "STOP5b_TCGA_HGG_validation_summary.csv"), row.names = FALSE)
saveRDS(val_tcga, file.path(proc_dir, "TCGA_HGG_validation_summary.rds"))

cat("\n== val_tcga ==\n")
print(val_tcga)

cat("\n########## STOP5b TCGA REPORT ##########\n")
cat("1. TCGA HGG = GBM G4 + LGG G3, grade-only landscape filter; strict GBM sensitivity not run here.\n")
cat("2. score_cohort coverage printed above and saved to STOP5_TCGA_HGG_score_coverage.csv.\n")
cat("3. val_tcga printed above and saved to STOP5b_TCGA_HGG_validation_summary.csv.\n")
cat("=> CGGA-325 validation remains next.\n")
