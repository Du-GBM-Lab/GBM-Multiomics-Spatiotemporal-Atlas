# R4 STOP1 audit for the old manuscript CGGA HGG RDS.
# Read-only: no ssGSEA, no clustering, no survival modeling.

root <- getwd()
raw_dir <- file.path(root, "data", "raw", "CGGA_old_manuscript")
out_dir <- file.path(root, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

hgg_path <- file.path(raw_dir, "hgg_data.rds")
gbm_no_c4_path <- file.path(raw_dir, "gbm_data_without_c4.rds")
marker_path <- file.path(raw_dir, "markers_subtypes_corl.csv")
old_script_path <- file.path(raw_dir, "代码.R")

cat("PROJECT_ROOT:", root, "\n")
cat("RAW_DIR:", raw_dir, "\n")
cat("hgg_data.rds exists:", file.exists(hgg_path), "\n")
cat("gbm_data_without_c4.rds exists:", file.exists(gbm_no_c4_path), "\n")
cat("markers_subtypes_corl.csv exists:", file.exists(marker_path), "\n")
cat("old script exists:", file.exists(old_script_path), "\n\n")

stopifnot(file.exists(hgg_path))
hgg <- readRDS(hgg_path)
stopifnot(is.list(hgg), all(c("clinical", "expression") %in% names(hgg)))

clin <- hgg$clinical
expr <- hgg$expression

cat("========== hgg_data.rds structure ==========\n")
cat("object class:", paste(class(hgg), collapse = ", "), "\n")
cat("names:", paste(names(hgg), collapse = ", "), "\n")
cat("clinical dim:", paste(dim(clin), collapse = " x "), "\n")
cat("expression dim:", paste(dim(expr), collapse = " x "), "\n")
cat("clinical columns:\n")
print(colnames(clin))
cat("expression first columns:\n")
print(head(colnames(expr), 10))

id_col <- "CGGA_ID"
grade_col <- "Grade"
prs_col <- "PRS_type"
hist_col <- "Histology"
idh_col <- "IDH_mutation_status"
codel_col <- "X1p19q_codeletion_status"
mgmt_col <- "MGMTp_methylation_status"
os_col <- "OS"
event_col <- "Censor..alive.0..dead.1."
age_col <- "Age"
sex_col <- "Gender"
radio_col <- "Radio_status..treated.1.un.treated.0."
chemo_col <- "Chemo_status..TMZ.treated.1.un.treated.0."

required <- c(id_col, grade_col, prs_col, hist_col, idh_col, codel_col, mgmt_col,
              os_col, event_col, age_col, sex_col, radio_col, chemo_col)
missing_required <- setdiff(required, colnames(clin))
cat("\nmissing required clinical columns:", paste(missing_required, collapse = ", "), "\n")
stopifnot(length(missing_required) == 0)

clean_missing <- function(x) {
  is.na(x) | trimws(as.character(x)) %in% c("", "NA", "N/A", "null", "NULL", "--")
}
miss_rate <- function(x) round(mean(clean_missing(x)), 3)

cat("\n========== clinical tables ==========\n")
for (col in c(prs_col, hist_col, grade_col, sex_col, idh_col, codel_col, mgmt_col,
              radio_col, chemo_col, event_col)) {
  cat("\n[", col, "] missing=", miss_rate(clin[[col]]), "\n", sep = "")
  print(table(clin[[col]], useNA = "ifany"))
}

cat("\n========== numeric summaries ==========\n")
for (col in c(age_col, os_col)) {
  cat("\n[", col, "] missing=", miss_rate(clin[[col]]), "\n", sep = "")
  print(summary(suppressWarnings(as.numeric(clin[[col]]))))
}

cat("\n========== cross tabs ==========\n")
cat("\nGrade x PRS_type\n")
print(table(clin[[grade_col]], clin[[prs_col]], useNA = "ifany"))
cat("\nGrade x IDH\n")
print(table(clin[[grade_col]], clin[[idh_col]], useNA = "ifany"))
cat("\nGrade x 1p19q\n")
print(table(clin[[grade_col]], clin[[codel_col]], useNA = "ifany"))
cat("\nIDH x 1p19q\n")
print(table(clin[[idh_col]], clin[[codel_col]], useNA = "ifany"))
cat("\nIDH x MGMT\n")
print(table(clin[[idh_col]], clin[[mgmt_col]], useNA = "ifany"))

os <- suppressWarnings(as.numeric(clin[[os_col]]))
ev <- suppressWarnings(as.numeric(clin[[event_col]]))
cat("\n========== OS coding ==========\n")
cat("events (=1):", sum(ev == 1, na.rm = TRUE), "\n")
cat("censored (=0):", sum(ev == 0, na.rm = TRUE), "\n")
cat("median OS:", median(os, na.rm = TRUE), "\n")
cat("OS unit guess:", ifelse(median(os, na.rm = TRUE) > 100, "days", "not obviously days"), "\n")

cat("\n========== expression audit ==========\n")
has_gene_name <- "Gene_Name" %in% colnames(expr)
cat("Gene_Name column exists:", has_gene_name, "\n")
expr_samples <- if (has_gene_name) setdiff(colnames(expr), "Gene_Name") else colnames(expr)
cat("expression samples:", length(expr_samples), "\n")
cat("clinical samples:", nrow(clin), "\n")
cat("direct sample overlap:", length(intersect(expr_samples, as.character(clin[[id_col]]))), "\n")
cat("first genes:", paste(head(if (has_gene_name) expr$Gene_Name else rownames(expr), 8), collapse = ", "), "\n")
vals <- suppressWarnings(as.numeric(unlist(expr[seq_len(min(20, nrow(expr))), setdiff(seq_len(ncol(expr)), which(colnames(expr) == "Gene_Name"))[seq_len(min(20, length(expr_samples)))]])))
cat("expression value range peek:", paste(round(range(vals, na.rm = TRUE), 3), collapse = " - "), "\n")
cat("expression value median peek:", round(median(vals, na.rm = TRUE), 3), "\n")

if (file.exists(marker_path)) {
  markers <- read.csv(marker_path, stringsAsFactors = FALSE)
  cat("\n========== old marker file ==========\n")
  cat("marker dim:", paste(dim(markers), collapse = " x "), "\n")
  cat("marker columns:", paste(colnames(markers), collapse = " | "), "\n")
  if (all(c("gene", "cluster") %in% colnames(markers))) {
    cat("marker genes:", length(unique(markers$gene)), "\n")
    cat("marker clusters:\n")
    print(table(markers$cluster, useNA = "ifany"))
    available <- intersect(unique(markers$gene), as.character(expr$Gene_Name))
    cat("old marker genes available in expression:", length(available), "\n")
  }
}

summary_rows <- data.frame(
  item = c("clinical_n", "expression_sample_n", "gene_n", "WHO_III", "WHO_IV",
           "IDH_mutant", "IDH_wildtype", "IDH_missing", "OS_events", "OS_censored",
           "median_OS", "expr_clin_overlap"),
  value = c(
    nrow(clin),
    length(expr_samples),
    nrow(expr),
    sum(clin[[grade_col]] == "WHO III", na.rm = TRUE),
    sum(clin[[grade_col]] == "WHO IV", na.rm = TRUE),
    sum(clin[[idh_col]] == "Mutant", na.rm = TRUE),
    sum(clin[[idh_col]] == "Wildtype", na.rm = TRUE),
    sum(clean_missing(clin[[idh_col]])),
    sum(ev == 1, na.rm = TRUE),
    sum(ev == 0, na.rm = TRUE),
    median(os, na.rm = TRUE),
    length(intersect(expr_samples, as.character(clin[[id_col]])))
  )
)
write.csv(summary_rows, file.path(out_dir, "STOP1_旧原稿CGGA_HGG审计_summary.csv"), row.names = FALSE)
write.csv(data.frame(column = colnames(clin),
                     missing_rate = sapply(clin, miss_rate)),
          file.path(out_dir, "STOP1_旧原稿CGGA_HGG临床字段缺失率.csv"), row.names = FALSE)

cat("\n########## STOP1 REPORT ##########\n")
cat("1. Old manuscript data source: hgg_data.rds, processed CGGA HGG-like cohort.\n")
cat("2. Cohort n=", nrow(clin), "; Grade WHO III=", sum(clin[[grade_col]] == "WHO III", na.rm = TRUE),
    ", WHO IV=", sum(clin[[grade_col]] == "WHO IV", na.rm = TRUE), ".\n", sep = "")
cat("3. IDH: Mutant=", sum(clin[[idh_col]] == "Mutant", na.rm = TRUE),
    ", Wildtype=", sum(clin[[idh_col]] == "Wildtype", na.rm = TRUE),
    ", missing=", sum(clean_missing(clin[[idh_col]])), ".\n", sep = "")
cat("4. OS: events=", sum(ev == 1, na.rm = TRUE),
    ", censored=", sum(ev == 0, na.rm = TRUE),
    ", median=", median(os, na.rm = TRUE), " days; event coding dead=1 if source label is correct.\n", sep = "")
cat("5. Expression: ", nrow(expr), " genes x ", length(expr_samples),
    " samples; gene ID=Symbol/Gene_Name; overlap with clinical IDs=",
    length(intersect(expr_samples, as.character(clin[[id_col]]))), ".\n", sep = "")
cat("=> STOP1 complete. Do not proceed to ssGSEA/clustering until reviewed.\n")
