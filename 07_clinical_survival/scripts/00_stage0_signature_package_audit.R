# =============================================================================
# R4 · 临床生存分析 —— STAGE 0：包与 signature 审计
# -----------------------------------------------------------------------------
# 目的：在下载/建模前确认 R4 的工具链与 signature 定义。
# 边界：R4 只做 MES / mesenchymal program rationale，不提 PLAUR，不含 FOSL1/FOSL2。
# 输出：signature 表、包版本表、禁用基因检查表。
# =============================================================================

## ===== CONFIG =====
ROOT <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/R4_临床生存分析"
NEFTEL_XLSX <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/05_恶性细胞分亚群与Neftel对照/data/Neftel2019_TableS2.xlsx"
OUT_TABLES <- file.path(ROOT, "tables")
dir.create(OUT_TABLES, recursive = TRUE, showWarnings = FALSE)

FORBIDDEN_GENES <- c("PLAUR", "FOSL1", "FOSL2")

MESV_GREEN_CORE <- c(
  "VIM", "FN1", "TIMP1", "ANXA2", "LGALS1", "S100A6", "S100A11", "COL6A1"
)
## ==================

pkg_ver <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    as.character(utils::packageVersion(pkg))
  } else {
    "NOT_INSTALLED"
  }
}

packages <- c(
  "TCGAbiolinks", "SummarizedExperiment", "GSVA", "survival", "survminer",
  "maxstat", "readxl", "dplyr", "tidyr", "tibble", "ggplot2", "qs2"
)
pkg_tbl <- data.frame(package = packages, version = vapply(packages, pkg_ver, character(1)))
write.csv(pkg_tbl, file.path(OUT_TABLES, "stage0_package_versions.csv"), row.names = FALSE)
print(pkg_tbl)

stopifnot(file.exists(NEFTEL_XLSX))
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("readxl is required to extract Neftel MES genes. Install readxl and rerun STAGE 0.")
}

# Neftel Table S2 structure:
# rows 1-4 are notes; row 5 contains module headers; below rows are genes.
neftel <- readxl::read_excel(NEFTEL_XLSX, sheet = "Table S2", skip = 4)
neftel_cols <- colnames(neftel)
stopifnot(all(c("MES1", "MES2") %in% neftel_cols))

neftel_mes <- unique(c(neftel$MES1, neftel$MES2))
neftel_mes <- toupper(trimws(as.character(neftel_mes)))
neftel_mes <- neftel_mes[!is.na(neftel_mes) & neftel_mes != ""]

signature_tbl <- data.frame(
  signature = c(rep("MESV_green_core", length(MESV_GREEN_CORE)), rep("Neftel_MES", length(neftel_mes))),
  gene = c(MESV_GREEN_CORE, neftel_mes)
)
signature_tbl$gene <- toupper(signature_tbl$gene)

forbidden_hits <- signature_tbl[signature_tbl$gene %in% FORBIDDEN_GENES, , drop = FALSE]
write.csv(signature_tbl, file.path(OUT_TABLES, "stage0_signatures_no_PLAUR.csv"), row.names = FALSE)
write.csv(forbidden_hits, file.path(OUT_TABLES, "stage0_forbidden_gene_hits.csv"), row.names = FALSE)

cat("\nSignature sizes:\n")
print(table(signature_tbl$signature))
cat("\nForbidden gene hits (must be 0 rows):\n")
print(forbidden_hits)

if (nrow(forbidden_hits) > 0) {
  stop("Forbidden genes detected in R4 signatures. Remove PLAUR/FOSL1/FOSL2 before proceeding.")
}

cat("\nSTOP 0 complete. Report package table + signature sizes before proceeding.\n")

