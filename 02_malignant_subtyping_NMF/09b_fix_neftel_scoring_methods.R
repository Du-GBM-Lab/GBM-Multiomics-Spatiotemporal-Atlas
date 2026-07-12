# 05_恶性细胞分亚群与Neftel对照/09b_fix_neftel_scoring_methods.R
# Fix Neftel scoring-derived labels before NMF/marker heatmap.
#
# Source of signatures:
# Neftel et al., Cell 2019, "An Integrative Model of Cellular States,
# Plasticity, and Genetics for Glioblastoma", doi:10.1016/j.cell.2019.06.024,
# Supplementary Table S2. The source Excel file is copied into ./data for
# project-local provenance and future reproducibility.
#
# Key fixes:
# 1) 4-state dominance is assigned by 6-way submodule which.max first, then
#    collapsing MES1/MES2 -> MES and NPC1/NPC2 -> NPC.
# 2) 2D Neftel placement uses y = max(OPC, NPC) - max(AC, MES), following
#    Neftel STAR Methods. State labels from this placement are audit labels.
# 3) Cycling is pmax(G1S, G2M), not G1S + G2M. AMS threshold uses > 1.
#    UCell threshold is the 95th percentile of random G1S/G2M-size gene sets.
# 4) MES1/MES2 binary is only assigned inside MES-dominant cells.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(qs2)
  library(Seurat)
  library(UCell)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
})

set.seed(42)

params <- list(
  project_dir = "05_恶性细胞分亚群与Neftel对照",
  input_object = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.submodule_labeled.qs2"
  ),
  output_object = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.v2.qs2"
  ),
  source_signature_xlsx = "<DATA_ROOT>/zetora/storage/7SGY8BYX/亚型marker.xlsx",
  local_signature_xlsx = file.path(
    "05_恶性细胞分亚群与Neftel对照",
    "data",
    "Neftel2019_TableS2.xlsx"
  ),
  tables_dir = file.path("05_恶性细胞分亚群与Neftel对照", "tables"),
  subtype_col = "subtype_k4",
  sample_col = "Pt_number",
  subtype_levels = paste0("Subtype", 1:4),
  ams_cycling_threshold = 1.0,
  ucell_random_pairs = 50,
  ucell_random_threshold_quantile = 0.95,
  mes_tie_tolerance = 0.05
)

dir.create(dirname(params$local_signature_xlsx), showWarnings = FALSE, recursive = TRUE)
dir.create(params$tables_dir, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

write_session_info <- function(path) {
  con <- file(path, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  sink(con)
  print(sessionInfo())
  sink()
}

read_neftel_signatures <- function(path) {
  raw <- readxl::read_excel(path, sheet = "Table S2", skip = 3, col_types = "text")
  names(raw) <- gsub("\\s+", "", names(raw))
  names(raw) <- gsub("/", "", names(raw))
  expected <- c("MES2", "MES1", "AC", "OPC", "NPC1", "NPC2", "G1S", "G2M")
  missing <- setdiff(expected, names(raw))
  if (length(missing) > 0) {
    stop("Missing signature columns in Excel: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  sigs <- lapply(expected, function(module) {
    genes <- trimws(as.character(raw[[module]]))
    unique(genes[!is.na(genes) & genes != ""])
  })
  names(sigs) <- expected
  sigs[c("MES1", "MES2", "AC", "OPC", "NPC1", "NPC2", "G1S", "G2M")]
}

six_way_assign <- function(score_df, method) {
  modules <- c("MES1", "MES2", "AC", "OPC", "NPC1", "NPC2")
  cols <- if (method == "AMS") paste0("AMS_", modules) else paste0(modules, "_UCell")
  mat <- as.matrix(score_df[, cols])
  max_idx <- max.col(mat, ties.method = "first")
  submodule <- modules[max_idx]
  state <- dplyr::case_when(
    submodule %in% c("MES1", "MES2") ~ "MES",
    submodule %in% c("NPC1", "NPC2") ~ "NPC",
    TRUE ~ submodule
  )
  list(submodule = submodule, state = state, max_score = mat[cbind(seq_len(nrow(mat)), max_idx)])
}

calc_2d <- function(score_df, method) {
  if (method == "AMS") {
    mes1 <- score_df$AMS_MES1
    mes2 <- score_df$AMS_MES2
    ac <- score_df$AMS_AC
    opc <- score_df$AMS_OPC
    npc1 <- score_df$AMS_NPC1
    npc2 <- score_df$AMS_NPC2
  } else {
    mes1 <- score_df$MES1_UCell
    mes2 <- score_df$MES2_UCell
    ac <- score_df$AC_UCell
    opc <- score_df$OPC_UCell
    npc1 <- score_df$NPC1_UCell
    npc2 <- score_df$NPC2_UCell
  }
  mes <- pmax(mes1, mes2)
  npc <- pmax(npc1, npc2)
  y <- pmax(opc, npc) - pmax(ac, mes)
  x <- ifelse(
    y > 0,
    ifelse(opc >= npc, -log2(abs(opc - npc) + 1), log2(abs(opc - npc) + 1)),
    ifelse(ac >= mes, -log2(abs(ac - mes) + 1), log2(abs(ac - mes) + 1))
  )
  state <- ifelse(
    y > 0,
    ifelse(opc >= npc, "OPC", "NPC"),
    ifelse(y < 0, ifelse(ac >= mes, "AC", "MES"), "Boundary")
  )
  tibble(x = x, y = y, state = state, mes = mes, npc = npc)
}

cliffs_delta <- function(x, y) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  ranks <- rank(c(x, y), ties.method = "average")
  rx <- sum(ranks[seq_along(x)])
  n_x <- length(x)
  n_y <- length(y)
  u <- rx - n_x * (n_x + 1) / 2
  (2 * u / (n_x * n_y)) - 1
}

score_test_one <- function(md, score_col) {
  s3 <- md[[score_col]][md$subtype_k4 == "Subtype3"]
  s4 <- md[[score_col]][md$subtype_k4 == "Subtype4"]
  wt <- suppressWarnings(wilcox.test(s3, s4))
  tibble(
    score = score_col,
    n_S3 = sum(!is.na(s3)),
    n_S4 = sum(!is.na(s4)),
    mean_S3 = mean(s3, na.rm = TRUE),
    mean_S4 = mean(s4, na.rm = TRUE),
    median_S3 = median(s3, na.rm = TRUE),
    median_S4 = median(s4, na.rm = TRUE),
    delta_mean_S3_minus_S4 = mean_S3 - mean_S4,
    wilcox_p = wt$p.value,
    cliffs_delta_S3_vs_S4 = cliffs_delta(s3, s4)
  )
}

copy_signature_file <- function() {
  if (!file.exists(params$local_signature_xlsx)) {
    if (!file.exists(params$source_signature_xlsx)) {
      stop("Neither local nor source signature xlsx exists.", call. = FALSE)
    }
    file.copy(params$source_signature_xlsx, params$local_signature_xlsx, overwrite = FALSE)
  }
}

msg("Copying/checking local Neftel signature file.")
copy_signature_file()
signature_info <- file.info(params$local_signature_xlsx)
signature_sha <- unname(tools::sha256sum(params$local_signature_xlsx))
neftel_sigs_raw <- read_neftel_signatures(params$local_signature_xlsx)

msg("Loading object:", params$input_object)
stopifnot(file.exists(params$input_object))
obj <- qs2::qs_read(params$input_object)
md <- obj@meta.data
stopifnot(identical(rownames(obj@meta.data), rownames(md)))

required_cols <- c(
  params$sample_col,
  params$subtype_col,
  "AMS_MES1", "AMS_MES2", "AMS_AC", "AMS_OPC", "AMS_NPC1", "AMS_NPC2", "AMS_G1S", "AMS_G2M",
  "MES1_UCell", "MES2_UCell", "AC_UCell", "OPC_UCell", "NPC1_UCell", "NPC2_UCell", "G1S_UCell", "G2M_UCell"
)
missing_cols <- setdiff(required_cols, colnames(md))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

md$cell_id <- rownames(md)
md$subtype_k4 <- factor(as.character(md[[params$subtype_col]]), levels = params$subtype_levels)

gene_recovery <- tibble(
  module = names(neftel_sigs_raw),
  n_genes_in_excel = vapply(neftel_sigs_raw, length, integer(1)),
  n_genes_used = vapply(neftel_sigs_raw, function(g) sum(g %in% rownames(obj)), integer(1)),
  recovery_pct = round(100 * n_genes_used / n_genes_in_excel, 2)
)

six_ams <- six_way_assign(md, "AMS")
six_ucell <- six_way_assign(md, "UCell")
two_ams <- calc_2d(md, "AMS")
two_ucell <- calc_2d(md, "UCell")

md$neftel_submodule_6way_AMS <- six_ams$submodule
md$neftel_state_4way_AMS <- six_ams$state
md$neftel_submodule_6way_UCell <- six_ucell$submodule
md$neftel_state_4way_UCell <- six_ucell$state

md$neftel_2D_x_AMS <- two_ams$x
md$neftel_2D_y_AMS <- two_ams$y
md$neftel_state_2D_AMS <- two_ams$state
md$neftel_2D_x_UCell <- two_ucell$x
md$neftel_2D_y_UCell <- two_ucell$y
md$neftel_state_2D_UCell <- two_ucell$state

md$cycling_score_AMS <- pmax(md$AMS_G1S, md$AMS_G2M)
md$cycling_score_UCell <- pmax(md$G1S_UCell, md$G2M_UCell)
md$is_cycling_AMS <- md$cycling_score_AMS > params$ams_cycling_threshold

set.seed(42)
all_genes <- rownames(obj)
g1s_size <- sum(neftel_sigs_raw$G1S %in% all_genes)
g2m_size <- sum(neftel_sigs_raw$G2M %in% all_genes)
random_sigs <- list()
for (i in seq_len(params$ucell_random_pairs)) {
  random_sigs[[paste0("random_G1S_", i)]] <- sample(all_genes, g1s_size)
  random_sigs[[paste0("random_G2M_", i)]] <- sample(all_genes, g2m_size)
}

msg("Running UCell random gene-set threshold with ", params$ucell_random_pairs, " random pairs.")
ucell_ncores <- if (.Platform$OS.type == "windows") 1 else 4
obj_random <- UCell::AddModuleScore_UCell(obj, features = random_sigs, ncores = ucell_ncores)
random_md <- obj_random@meta.data
random_pair_scores <- lapply(seq_len(params$ucell_random_pairs), function(i) {
  pmax(random_md[[paste0("random_G1S_", i, "_UCell")]], random_md[[paste0("random_G2M_", i, "_UCell")]])
})
random_scores <- unlist(random_pair_scores, use.names = FALSE)
ucell_threshold <- as.numeric(quantile(random_scores, probs = params$ucell_random_threshold_quantile, na.rm = TRUE))
md$is_cycling_UCell <- md$cycling_score_UCell > ucell_threshold

md$mes_submodule_within_mes_AMS <- ifelse(
  md$neftel_state_4way_AMS == "MES",
  ifelse(md$AMS_MES1 >= md$AMS_MES2, "MES1", "MES2"),
  NA_character_
)
md$mes_submodule_within_mes_UCell <- ifelse(
  md$neftel_state_4way_UCell == "MES",
  ifelse(md$MES1_UCell >= md$MES2_UCell, "MES1", "MES2"),
  NA_character_
)
md$MES1_MES2_absdiff_AMS <- abs(md$AMS_MES1 - md$AMS_MES2)
md$MES1_MES2_absdiff_UCell <- abs(md$MES1_UCell - md$MES2_UCell)
md$MES1_MES2_tie_like_AMS <- md$MES1_MES2_absdiff_AMS < params$mes_tie_tolerance
md$MES1_MES2_tie_like_UCell <- md$MES1_MES2_absdiff_UCell < params$mes_tie_tolerance

state_cols <- c(
  "neftel_submodule_6way_AMS", "neftel_state_4way_AMS",
  "neftel_submodule_6way_UCell", "neftel_state_4way_UCell",
  "neftel_2D_x_AMS", "neftel_2D_y_AMS", "neftel_state_2D_AMS",
  "neftel_2D_x_UCell", "neftel_2D_y_UCell", "neftel_state_2D_UCell",
  "cycling_score_AMS", "cycling_score_UCell", "is_cycling_AMS", "is_cycling_UCell",
  "mes_submodule_within_mes_AMS", "mes_submodule_within_mes_UCell",
  "MES1_MES2_absdiff_AMS", "MES1_MES2_absdiff_UCell",
  "MES1_MES2_tie_like_AMS", "MES1_MES2_tie_like_UCell"
)

stopifnot(identical(rownames(obj@meta.data), rownames(md)))
for (col in state_cols) {
  obj@meta.data[[col]] <- md[[col]]
}
obj@meta.data[[params$subtype_col]] <- md$subtype_k4

per_cell_cols <- intersect(
  c(
    "cell_id", params$sample_col, params$subtype_col, "Phase",
    "AMS_MES1", "AMS_MES2", "AMS_AC", "AMS_OPC", "AMS_NPC1", "AMS_NPC2", "AMS_G1S", "AMS_G2M",
    "MES1_UCell", "MES2_UCell", "AC_UCell", "OPC_UCell", "NPC1_UCell", "NPC2_UCell", "G1S_UCell", "G2M_UCell",
    state_cols
  ),
  colnames(md)
)
per_cell_scores <- md[, per_cell_cols] |>
  as_tibble()

score_cols_test <- c(
  "AMS_MES1", "AMS_MES2", "AMS_AC", "AMS_OPC", "AMS_NPC1", "AMS_NPC2", "AMS_G1S", "AMS_G2M",
  "MES1_UCell", "MES2_UCell", "AC_UCell", "OPC_UCell", "NPC1_UCell", "NPC2_UCell", "G1S_UCell", "G2M_UCell",
  "cycling_score_AMS", "cycling_score_UCell", "neftel_2D_x_AMS", "neftel_2D_y_AMS", "neftel_2D_x_UCell", "neftel_2D_y_UCell"
)
s3_s4_tests <- bind_rows(lapply(score_cols_test, function(col) score_test_one(md, col))) |>
  mutate(wilcox_p_adj_BH = p.adjust(wilcox_p, method = "BH")) |>
  arrange(desc(abs(cliffs_delta_S3_vs_S4)))

score_by_sample_subtype <- md |>
  group_by(sample = .data[[params$sample_col]], subtype_k4) |>
  summarise(
    n_cells = n(),
    across(all_of(score_cols_test), \(x) mean(x, na.rm = TRUE), .names = "{.col}_mean"),
    .groups = "drop"
  ) |>
  arrange(sample, subtype_k4)

within_patient_s3_s4_delta <- score_by_sample_subtype |>
  filter(subtype_k4 %in% c("Subtype3", "Subtype4")) |>
  select(sample, subtype_k4, n_cells, ends_with("_mean")) |>
  pivot_wider(
    names_from = subtype_k4,
    values_from = c(n_cells, ends_with("_mean")),
    names_sep = "__"
  ) |>
  filter(!is.na(n_cells__Subtype3), !is.na(n_cells__Subtype4)) |>
  mutate(
    n_cells_min_S3_S4 = pmin(n_cells__Subtype3, n_cells__Subtype4)
  )

for (col in score_cols_test) {
  s3_col <- paste0(col, "_mean__Subtype3")
  s4_col <- paste0(col, "_mean__Subtype4")
  if (all(c(s3_col, s4_col) %in% colnames(within_patient_s3_s4_delta))) {
    within_patient_s3_s4_delta[[paste0(col, "_delta_S3_minus_S4")]] <-
      within_patient_s3_s4_delta[[s3_col]] - within_patient_s3_s4_delta[[s4_col]]
  }
}

within_delta_summary <- within_patient_s3_s4_delta |>
  select(sample, starts_with("n_cells"), ends_with("_delta_S3_minus_S4")) |>
  pivot_longer(
    cols = ends_with("_delta_S3_minus_S4"),
    names_to = "score_delta",
    values_to = "delta_S3_minus_S4"
  ) |>
  group_by(score_delta) |>
  summarise(
    n_patients_with_both = n(),
    n_positive = sum(delta_S3_minus_S4 > 0, na.rm = TRUE),
    n_negative = sum(delta_S3_minus_S4 < 0, na.rm = TRUE),
    median_delta = median(delta_S3_minus_S4, na.rm = TRUE),
    mean_delta = mean(delta_S3_minus_S4, na.rm = TRUE),
    .groups = "drop"
  )

state_summary <- bind_rows(
  lapply(c("neftel_state_4way_AMS", "neftel_state_4way_UCell", "neftel_state_2D_AMS", "neftel_state_2D_UCell"), function(col) {
    md |>
      count(subtype_k4, state = .data[[col]], name = "n_cells") |>
      group_by(subtype_k4) |>
      mutate(row_pct = round(100 * n_cells / sum(n_cells), 2), state_column = col) |>
      ungroup() |>
      select(state_column, subtype_k4, state, n_cells, row_pct)
  })
)

old_new_state_comparison <- bind_rows(
  lapply(
    c(
      old_AMS = "neftel_state_AMS",
      old_UCell = "neftel_state_UCell",
      new_4way_AMS = "neftel_state_4way_AMS",
      new_4way_UCell = "neftel_state_4way_UCell",
      new_2D_AMS = "neftel_state_2D_AMS",
      new_2D_UCell = "neftel_state_2D_UCell"
    ),
    function(col) {
      if (!col %in% colnames(md)) return(NULL)
      md |>
        count(subtype_k4, state = .data[[col]], name = "n_cells") |>
        group_by(subtype_k4) |>
        mutate(row_pct = round(100 * n_cells / sum(n_cells), 2), state_column = col) |>
        ungroup()
    }
  ),
  .id = "source"
)

submodule_spearman <- bind_rows(lapply(c("MES1", "MES2", "AC", "OPC", "NPC1", "NPC2", "G1S", "G2M"), function(module) {
  ams_col <- paste0("AMS_", module)
  ucell_col <- paste0(module, "_UCell")
  tibble(
    module = module,
    spearman_AMS_UCell = cor(md[[ams_col]], md[[ucell_col]], method = "spearman", use = "complete.obs")
  )
}))

mes_binary_summary <- bind_rows(
  lapply(c("AMS", "UCell"), function(method) {
    state_col <- paste0("neftel_state_4way_", method)
    mes_col <- paste0("mes_submodule_within_mes_", method)
    tie_col <- paste0("MES1_MES2_tie_like_", method)
    md |>
      group_by(subtype_k4) |>
      summarise(
        method = method,
        n_cells = n(),
        n_mes_dominant = sum(.data[[state_col]] == "MES", na.rm = TRUE),
        mes_dominant_pct = round(100 * n_mes_dominant / n_cells, 2),
        MES1_within_mes_pct = round(100 * mean(.data[[mes_col]] == "MES1", na.rm = TRUE), 2),
        MES2_within_mes_pct = round(100 * mean(.data[[mes_col]] == "MES2", na.rm = TRUE), 2),
        tie_like_all_cells_pct = round(100 * mean(.data[[tie_col]], na.rm = TRUE), 2),
        .groups = "drop"
      )
  })
)

sanity_checks <- tibble(
  check = c(
    "n_cells",
    "n_features",
    "n_samples",
    "subtype_na",
    "score_na_total",
    "new_state_na_total",
    "per_cell_rows",
    "ucell_random_threshold",
    "ucell_random_pairs",
    "signature_file_exists",
    "signature_file_size_bytes"
  ),
  value = c(
    ncol(obj),
    nrow(obj),
    n_distinct(md[[params$sample_col]]),
    sum(is.na(md$subtype_k4)),
    sum(is.na(md[, required_cols[required_cols %in% colnames(md)]])),
    sum(is.na(md[, state_cols[!grepl("^mes_submodule_within_mes", state_cols)]])),
    nrow(per_cell_scores),
    ucell_threshold,
    params$ucell_random_pairs,
    file.exists(params$local_signature_xlsx),
    signature_info$size
  )
)

provenance <- tibble(
  item = c(
    "source",
    "doi",
    "local_file",
    "sha256",
    "file_size_bytes",
    "sheet",
    "skip_rows",
    "n_random_pairs_for_UCell_threshold",
    "UCell_random_95pct_threshold"
  ),
  value = c(
    "Neftel et al., Cell 2019, Supplementary Table S2",
    "10.1016/j.cell.2019.06.024",
    params$local_signature_xlsx,
    signature_sha,
    signature_info$size,
    "Table S2",
    "3",
    params$ucell_random_pairs,
    ucell_threshold
  )
)

msg("Saving v2 object:", params$output_object)
qs2::qs_save(obj, params$output_object)

write_tsv(per_cell_scores, file.path(params$tables_dir, "09b_per_cell_scores.tsv"))
write_csv(s3_s4_tests, file.path(params$tables_dir, "09b_subtype3_vs_subtype4_score_test.csv"))
write_csv(score_by_sample_subtype, file.path(params$tables_dir, "09b_score_by_sample_and_subtype.csv"))
write_csv(within_patient_s3_s4_delta, file.path(params$tables_dir, "09b_within_patient_S3_S4_score_delta.csv"))
write_csv(within_delta_summary, file.path(params$tables_dir, "09b_within_patient_S3_S4_delta_summary.csv"))
write_csv(state_summary, file.path(params$tables_dir, "09b_subtype_x_neftel_state_v2.csv"))
write_csv(old_new_state_comparison, file.path(params$tables_dir, "09b_old_vs_v2_state_comparison.csv"))
write_csv(submodule_spearman, file.path(params$tables_dir, "09b_AMS_UCell_submodule_spearman.csv"))
write_csv(mes_binary_summary, file.path(params$tables_dir, "09b_MES1_MES2_within_MES_summary.csv"))
write_csv(gene_recovery, file.path(params$tables_dir, "09b_neftel_module_gene_recovery.csv"))
write_csv(provenance, file.path(params$tables_dir, "09b_neftel_signature_provenance.csv"))
write_csv(sanity_checks, file.path(params$tables_dir, "09b_sanity_checks.csv"))
write_session_info(file.path(params$tables_dir, "09b_session_info.txt"))

cat("\nSanity checks:\n")
print(sanity_checks)
cat("\nGene recovery:\n")
print(gene_recovery)
cat("\nSubtype x Neftel state v2 row %:\n")
print(state_summary)
cat("\nTop S3 vs S4 effect sizes:\n")
print(head(s3_s4_tests, 8))
cat("\nPatients with both S3 and S4:", nrow(within_patient_s3_s4_delta), "\n")
cat("\nWithin-patient delta summary, top rows:\n")
print(head(within_delta_summary, 8))
msg("Done.")
