# =============================================================================
# R4 · 临床生存分析 —— STAGE 1：TCGA 下载与结构审计准备
# -----------------------------------------------------------------------------
# 默认不下载。先设置 RUN_DOWNLOAD <- TRUE 后才会调用 GDCdownload。
# 当前 TCGAbiolinks RNA-seq workflow.type 使用 "STAR - Counts"。
# STOP 1：回报 query 样本数、clinical 字段、可用 OS/IDH/MGMT 字段，再决定建模口径。
# =============================================================================

## ===== CONFIG =====
ROOT <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/R4_临床生存分析"
RAW_DIR <- file.path(ROOT, "data", "raw", "TCGA_GBM")
PROC_DIR <- file.path(ROOT, "data", "processed")
TAB_DIR <- file.path(ROOT, "tables")
dir.create(RAW_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PROC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

RUN_DOWNLOAD <- FALSE
## ==================

needed <- c("TCGAbiolinks", "SummarizedExperiment", "dplyr", "tibble", "qs2")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
})

query <- TCGAbiolinks::GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)

query_results <- TCGAbiolinks::getResults(query)
write.csv(query_results, file.path(TAB_DIR, "stage1_tcga_gbm_star_counts_query_results.csv"), row.names = FALSE)
cat("TCGA-GBM STAR - Counts primary tumor query rows:", nrow(query_results), "\n")

clinical <- TCGAbiolinks::GDCquery_clinic(project = "TCGA-GBM", type = "clinical")
write.csv(clinical, file.path(TAB_DIR, "stage1_tcga_gbm_clinical_raw.csv"), row.names = FALSE)
cat("Clinical rows:", nrow(clinical), " columns:", ncol(clinical), "\n")
cat("Clinical columns:\n")
print(colnames(clinical))

if (RUN_DOWNLOAD) {
  TCGAbiolinks::GDCdownload(query, directory = RAW_DIR, method = "api", files.per.chunk = 50)
  se <- TCGAbiolinks::GDCprepare(query, directory = RAW_DIR)
  qs2::qs_save(se, file.path(PROC_DIR, "tcga_gbm_star_counts_primary_tumor_se.qs2"))
  cat("Downloaded and saved SummarizedExperiment.\n")
} else {
  cat("RUN_DOWNLOAD is FALSE. No GDCdownload/GDCprepare executed.\n")
}

cat("\nSTOP 1 complete. Report query count, clinical columns, and whether IDH/MGMT/OS fields are available.\n")

