# =============================================================================
# R4 · 临床队列生存分析 —— STOP1：CGGA 队列结构审计
# -----------------------------------------------------------------------------
# 只审计，不打分、不聚类、不读 gene set。
# 目的：确认 CGGA mRNAseq_693 / 325 的文件结构、字段名、样本口径、
#       grade 构成、OS 编码、关键协变量缺失率、表达矩阵与临床 ID 对齐。
#
# R4 边界：此脚本不涉及 PLAUR / FOSL1 / FOSL2，也不做任何机制 claim。
# =============================================================================

## ---- CONFIG：只改这里 --------------------------------------------------------
PROJECT_DIR <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/R4_临床生存分析"
CGGA_DIR    <- file.path(PROJECT_DIR, "data", "raw", "CGGA")
OUT         <- file.path(PROJECT_DIR, "tables")
dir.create(CGGA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# 若自动识别错文件，手动填这两行；留 NA 则自动按文件名关键词找。
EXPR_FILE_OVERRIDE <- NA_character_
CLIN_FILE_OVERRIDE <- NA_character_
## ----------------------------------------------------------------------------

has_dt <- requireNamespace("data.table", quietly = TRUE)
read_tab <- function(f, nrows = Inf) {
  if (has_dt) {
    data.table::fread(
      f,
      nrows = if (is.finite(nrows)) nrows else -1,
      data.table = FALSE,
      check.names = FALSE
    )
  } else {
    read.delim(
      f,
      check.names = FALSE,
      nrows = if (is.finite(nrows)) nrows else -1
    )
  }
}

clean_missing <- function(x) is.na(x) | trimws(as.character(x)) %in% c("", "NA", "N/A", "null", "NULL", "--")
miss_rate <- function(x) round(mean(clean_missing(x)), 3)

find_col <- function(df, pats) {
  for (p in pats) {
    h <- grep(p, colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(h)) return(h[1])
  }
  NA_character_
}

cat("PROJECT_DIR:", PROJECT_DIR, "\n")
cat("CGGA_DIR   :", CGGA_DIR, "\n\n")

files <- list.files(CGGA_DIR, full.names = TRUE)
cat("== CGGA_DIR files ==\n")
print(basename(files))

if (!length(files)) {
  stop("CGGA_DIR is empty. Put CGGA expression/clinical files here or set CGGA_DIR/override paths.")
}

if (is.na(EXPR_FILE_OVERRIDE)) {
  expr_file <- grep("expr|RSEM|FPKM|TPM|mRNAseq.*\\.txt|mRNAseq.*\\.csv", files, ignore.case = TRUE, value = TRUE)
  clin_like <- grep("clinic|clinical", files, ignore.case = TRUE, value = TRUE)
  expr_file <- setdiff(expr_file, clin_like)[1]
} else {
  expr_file <- EXPR_FILE_OVERRIDE
}

if (is.na(CLIN_FILE_OVERRIDE)) {
  clin_file <- grep("clinic|clinical", files, ignore.case = TRUE, value = TRUE)[1]
} else {
  clin_file <- CLIN_FILE_OVERRIDE
}

cat("\nExpression file:", expr_file, "\n")
cat("Clinical file  :", clin_file, "\n")
if (is.na(expr_file) || !file.exists(expr_file)) stop("Expression file not found.")
if (is.na(clin_file) || !file.exists(clin_file)) stop("Clinical file not found.")

## ---- 1. Clinical audit -------------------------------------------------------
clin <- read_tab(clin_file)
cat("\n== Clinical dimensions ==\n")
print(dim(clin))
cat("\n== Clinical column names ==\n")
print(colnames(clin))
cat("\n== Clinical head ==\n")
print(utils::head(clin, 3))

map <- list(
  id    = find_col(clin, c("^CGGA", "CGGA_ID", "sample", "^ID$", "case")),
  prs   = find_col(clin, c("PRS", "primary.*recur", "tumor.*type", "sample.*type")),
  grade = find_col(clin, c("WHO.*grade", "grade", "WHO")),
  hist  = find_col(clin, c("histolog")),
  idh   = find_col(clin, c("IDH.*mut", "IDH.?status", "^IDH")),
  codel = find_col(clin, c("1p.?19q", "codel")),
  mgmt  = find_col(clin, c("MGMT")),
  os    = find_col(clin, c("^OS$", "OS.*day", "overall.*surv", "surv.*time")),
  event = find_col(clin, c("censor", "vital", "dead", "event", "status$")),
  age   = find_col(clin, c("^age", "age")),
  sex   = find_col(clin, c("gender", "^sex"))
)
cat("\n== Auto field mapping (verify manually) ==\n")
print(unlist(map))
write.csv(data.frame(field = names(map), column = unlist(map)),
          file.path(OUT, "stage1_cgga_field_mapping.csv"), row.names = FALSE)

cat("\n== Key variable values + missing rates ==\n")
for (k in c("prs", "grade", "idh", "codel", "mgmt", "sex", "hist")) {
  col <- map[[k]]
  if (is.na(col)) {
    cat(sprintf("[%s] no matched column\n", k))
    next
  }
  cat(sprintf("\n[%s] -> %s  missing=%.3f\n", k, col, miss_rate(clin[[col]])))
  print(table(clin[[col]], useNA = "ifany"))
}

for (k in c("os", "age")) {
  col <- map[[k]]
  if (is.na(col)) next
  cat(sprintf("\n[%s] -> %s  missing=%.3f  numeric summary:\n", k, col, miss_rate(clin[[col]])))
  print(summary(suppressWarnings(as.numeric(clin[[col]]))))
}

## ---- 2. Grade / PRS / IDH cross tabs ----------------------------------------
if (!is.na(map$grade) && !is.na(map$prs)) {
  cat("\n== grade x PRS cross-tab ==\n")
  print(table(clin[[map$grade]], clin[[map$prs]], useNA = "ifany"))
}

if (!is.na(map$grade) && !is.na(map$idh)) {
  cat("\n== grade x IDH cross-tab ==\n")
  print(table(clin[[map$grade]], clin[[map$idh]], useNA = "ifany"))
}

## ---- 3. OS coding ------------------------------------------------------------
if (!is.na(map$event)) {
  cat("\n== event/censor values (confirm dead=1) ==\n")
  print(table(clin[[map$event]], useNA = "ifany"))
}

if (!is.na(map$os) && !is.na(map$event)) {
  os <- suppressWarnings(as.numeric(clin[[map$os]]))
  ev <- suppressWarnings(as.numeric(clin[[map$event]]))
  cat(sprintf(
    "\nEvents (=1)=%d  censored (=0)=%d  median OS=%.0f  unit guess=%s\n",
    sum(ev == 1, na.rm = TRUE),
    sum(ev == 0, na.rm = TRUE),
    median(os, na.rm = TRUE),
    ifelse(median(os, na.rm = TRUE) > 100, "days", "not obviously days")
  ))
}

## ---- 4. Expression matrix peek ----------------------------------------------
peek <- read_tab(expr_file, nrows = 5)
cat("\n== Expression peek dimensions (first 5 rows) ==\n")
print(dim(peek))
cat("\n== Expression column examples ==\n")
print(utils::head(colnames(peek), 10))

g1 <- as.character(peek[[1]])
cat("Gene ID examples:", paste(utils::head(g1, 3), collapse = ", "),
    "->", ifelse(any(grepl("^ENSG", g1)), "Ensembl", "Symbol suspected"), "\n")
vals <- suppressWarnings(as.numeric(unlist(peek[, -1])))
cat(sprintf(
  "Value range=[%.2f, %.2f], median=%.2f -> scale guess=%s\n",
  min(vals, na.rm = TRUE),
  max(vals, na.rm = TRUE),
  median(vals, na.rm = TRUE),
  ifelse(max(vals, na.rm = TRUE) > 1000, "raw count or unlogged abundance", "FPKM/TPM or already log-like")
))

# Full read only if needed to get exact dimensions and sample names.
expr_header <- read_tab(expr_file, nrows = 1)
expr_samples <- colnames(expr_header)[-1]
gene_ids <- read_tab(expr_file)[[1]]
cat("Genes (rows) =", length(gene_ids), " Samples (columns) =", length(expr_samples), "\n")

## ---- 5. Expression-column vs clinical-ID alignment ---------------------------
if (!is.na(map$id)) {
  clin_ids <- as.character(clin[[map$id]])
  ov_direct <- length(intersect(expr_samples, clin_ids))
  ov_upper <- length(intersect(toupper(expr_samples), toupper(clin_ids)))
  cat(sprintf(
    "\n== Sample ID alignment: expression=%d, clinical=%d, direct overlap=%d, upper-case overlap=%d ==\n",
    length(expr_samples), nrow(clin), ov_direct, ov_upper
  ))
  if (ov_direct == 0 && ov_upper == 0) {
    cat("!! IDs do not match directly. Check prefixes/suffixes/versioned sample IDs.\n")
  }
}

## ---- 6. Minimal output tables ------------------------------------------------
write.csv(clin, file.path(OUT, "stage1_cgga_clinical_raw_peeked.csv"), row.names = FALSE)
write.csv(data.frame(expr_sample = expr_samples),
          file.path(OUT, "stage1_cgga_expression_sample_ids.csv"), row.names = FALSE)

cat("\n\n########## STOP1 REPORT ##########\n")
cat("1. HGG口径: grade构成见上；拟定 HGG = ?；原发/复发筛 = ?\n")
cat("2. 筛 HGG(原发) 后 n = ?；其中 IDH-wt / IDH-mut = ?\n")
cat("3. 缺失率: IDH / 1p19q / MGMT / age / sex = ?；是否足够进 cluster 注释 + Block4 Cox 协变量?\n")
cat("4. OS: 事件数/删失、单位(days?)、event 编码(dead=1?) 见上。\n")
cat("5. 表达: 量纲(count/FPKM/TPM/log?)、基因ID(symbol/ensembl)、维度、与临床 ID 交集见上。\n")
cat("=> STOP1 complete. Do not proceed to ssGSEA/clustering until reviewed.\n")

