suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(readr)
  library(tidyr)
})

step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir) {
  setwd(file.path(getwd(), step_dir))
}

obj_path <- "outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"
out_cols <- "tables/F7_metadata_columns_audit.csv"
out_candidates <- "tables/F7_clinical_metadata_candidates.csv"
out_patient <- "tables/F7_patient_level_metadata_preview.csv"
out_subtype_comp <- "tables/F7_patient_subtype_composition_preview.csv"
out_session <- "tables/F7_metadata_audit_session_info.txt"

obj <- qs2::qs_read(obj_path)
md <- obj@meta.data
md$cell_id <- rownames(md)

subtype_col <- "subtype_k4"
patient_candidates <- c("Pt_number", "patient", "Patient", "patient_id", "Patient_ID", "case", "Case", "sample", "Sample", "orig.ident")
patient_col <- patient_candidates[patient_candidates %in% colnames(md)][1]
if (is.na(patient_col)) stop("No patient/sample column found among: ", paste(patient_candidates, collapse = ", "))
if (!subtype_col %in% colnames(md)) stop("Missing subtype_k4")

candidate_regex <- "(age|sex|gender|grade|who|idh|mgmt|tert|egfr|chr|1p|19q|surv|os|pfs|death|status|recur|therapy|treat|radio|chemo|clinical|hist|path|diagnos|tumor|location|sample|patient|pt_|pt|barcode|病|年龄|性别|生存|复发|治疗|分级|级别|诊断|部位|临床)"

col_audit <- lapply(colnames(md), function(col) {
  x <- md[[col]]
  x_chr <- as.character(x)
  n_unique <- dplyr::n_distinct(x_chr, na.rm = TRUE)
  n_missing <- sum(is.na(x) | x_chr == "")
  examples <- unique(x_chr[!(is.na(x) | x_chr == "")])
  examples <- head(examples, 8)
  data.frame(
    column = col,
    class = paste(class(x), collapse = ";"),
    n_unique = n_unique,
    n_missing = n_missing,
    missing_pct = round(n_missing / nrow(md) * 100, 3),
    example_values = paste(examples, collapse = " | "),
    candidate_clinical = grepl(candidate_regex, col, ignore.case = TRUE) || (n_unique > 1 && n_unique <= 40 && !col %in% c("cell_id")),
    stringsAsFactors = FALSE
  )
}) %>% bind_rows()
readr::write_csv(col_audit, out_cols)

candidate_cols <- col_audit %>%
  filter(candidate_clinical) %>%
  arrange(desc(grepl(candidate_regex, column, ignore.case = TRUE)), n_missing, n_unique)
readr::write_csv(candidate_cols, out_candidates)

patient_level_cols <- setdiff(
  candidate_cols$column[candidate_cols$column %in% colnames(md)],
  c(patient_col, subtype_col, "cell_id")
)

patient_level <- md %>%
  group_by(.data[[patient_col]]) %>%
  summarise(
    n_cells = n(),
    n_subtypes = n_distinct(.data[[subtype_col]]),
    across(
      .cols = all_of(patient_level_cols),
      .fns = ~ {
        vals <- unique(as.character(.x[!(is.na(.x) | as.character(.x) == "")]))
        paste(head(vals, 5), collapse = " | ")
      },
      .names = "{.col}"
    ),
    .groups = "drop"
  ) %>%
  rename(patient_id = 1)
readr::write_csv(patient_level, out_patient)

subtype_comp <- md %>%
  count(.data[[patient_col]], .data[[subtype_col]], name = "n_cells") %>%
  group_by(.data[[patient_col]]) %>%
  mutate(patient_total_cells = sum(n_cells), fraction = n_cells / patient_total_cells) %>%
  ungroup() %>%
  rename(patient_id = 1, subtype_k4 = 2) %>%
  arrange(patient_id, subtype_k4)
readr::write_csv(subtype_comp, out_subtype_comp)

writeLines(capture.output(sessionInfo()), out_session)

cat("F7 metadata audit completed.\n")
cat("Patient column:", patient_col, "\n")
cat("Cells:", nrow(md), "\n")
cat("Patients:", n_distinct(md[[patient_col]]), "\n")
cat("Metadata columns:", ncol(md), "\n")
cat("Candidate clinical columns:", nrow(candidate_cols), "\n")
cat("Top candidate columns:\n")
print(head(candidate_cols, 30))
cat("Patient-level preview:\n")
print(patient_level)
cat("Subtype composition preview:\n")
print(head(subtype_comp, 20))
