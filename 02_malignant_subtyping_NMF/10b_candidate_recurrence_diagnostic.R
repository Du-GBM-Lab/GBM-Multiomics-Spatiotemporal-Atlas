# Diagnostic only: recurrent metaprogram count across candidate cluster counts.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(dplyr)
  library(readr)
})

base <- "05_恶性细胞分亚群与Neftel对照"
hc <- readRDS(file.path(base, "outputs", "nmf", "10b_program_dendrogram.rds"))
program_meta <- read_csv(file.path(base, "tables", "10b_all_programs_top50_genes.csv"), show_col_types = FALSE) |>
  distinct(program_id, patient)

diagnostic <- lapply(6:20, function(k) {
  assignment <- cutree(hc, k = k)
  df <- tibble(
    program_id = names(assignment),
    cluster = sprintf("MP%02d", assignment)
  ) |>
    left_join(program_meta, by = "program_id")

  summary <- df |>
    group_by(cluster) |>
    summarise(
      n_programs = n(),
      n_patients_unique = n_distinct(patient),
      .groups = "drop"
    )

  tibble(
    candidate_k = k,
    n_recurrent_ge6 = sum(summary$n_patients_unique >= 6),
    n_nonrecurrent_lt6 = sum(summary$n_patients_unique < 6),
    min_patients = min(summary$n_patients_unique),
    max_patients = max(summary$n_patients_unique)
  )
}) |>
  bind_rows()

write_csv(diagnostic, file.path(base, "tables", "10b_candidate_recurrence_diagnostic.csv"))
print(diagnostic)
