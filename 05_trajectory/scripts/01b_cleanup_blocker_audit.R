suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
})

if (!dir.exists("06_恶性细胞拟时序") && basename(getwd()) != "06_恶性细胞拟时序") {
  stop("Run from 重新分析 root or 06_恶性细胞拟时序.")
}
if (dir.exists("06_恶性细胞拟时序")) {
  setwd("06_恶性细胞拟时序")
}

source(file.path("scripts", "_naming.R"))

obj_path <- file.path(
  "..",
  "05_恶性细胞分亚群与Neftel对照",
  "outputs",
  "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"
)
obj <- qs2::qs_read(obj_path)

pt <- as.character(obj@meta.data$Pt_number)
patient_ids <- sort(unique(pt))
patient_audit <- tibble(
  Pt_number_raw = patient_ids,
  Pt_number_trim = trimws(patient_ids),
  has_leading_or_trailing_space = Pt_number_raw != Pt_number_trim,
  has_internal_space = str_detect(Pt_number_trim, "\\s"),
  has_leading_zero = str_detect(Pt_number_trim, "^Pt0[0-9]"),
  normalized_upper = toupper(Pt_number_trim)
)
patient_summary <- tibble(
  n_unique_Pt_number = length(patient_ids),
  n_unique_after_trim = length(unique(trimws(pt))),
  any_space = any(patient_audit$has_leading_or_trailing_space | patient_audit$has_internal_space),
  any_leading_zero = any(patient_audit$has_leading_zero),
  patient_ids = paste(patient_ids, collapse = "; ")
)
write_csv(patient_audit, file.path("tables", "patient_id_format_audit.csv"))
write_csv(patient_summary, file.path("tables", "patient_count_summary.csv"))

blacklisted <- read_csv(file.path("tables", "trajectory_driver_gene_candidates_blacklisted.csv"), show_col_types = FALSE) |>
  left_join(
    subtype_naming_mapping |>
      select(subtype_k4, subtype_label_final, abbreviation),
    by = "subtype_k4"
  ) |>
  mutate(
    direction = case_when(
      avg_log2FC > 0 ~ paste0(as.character(subtype_k4), "_up"),
      avg_log2FC < 0 ~ paste0(as.character(subtype_k4), "_down"),
      TRUE ~ "zero"
    )
  ) |>
  select(
    subtype_k4,
    subtype_label_final,
    abbreviation,
    gene,
    direction,
    avg_log2FC,
    BH_q,
    pct.1,
    pct.2,
    is_blacklisted
  )
write_csv(blacklisted, file.path("tables", "blacklisted_DE_rows.csv"))

patterns <- readLines(file.path("data", "gene_blacklist.txt"), warn = FALSE)
patterns <- patterns[!str_detect(patterns, "^\\s*(#|$)")]
blacklist_summary <- tibble(
  n_blacklist_regex_patterns = length(patterns),
  regex_patterns = paste(patterns, collapse = "; "),
  n_blacklisted_rows_in_current_DE = nrow(blacklisted),
  n_blacklisted_unique_genes_in_current_DE = n_distinct(blacklisted$gene),
  blacklisted_genes_in_current_DE = paste(sort(unique(blacklisted$gene)), collapse = "; ")
)
write_csv(blacklist_summary, file.path("tables", "gene_blacklist_summary.csv"))

print(patient_summary)
print(patient_audit, n = Inf)
print(blacklisted, n = Inf)
print(blacklist_summary)
