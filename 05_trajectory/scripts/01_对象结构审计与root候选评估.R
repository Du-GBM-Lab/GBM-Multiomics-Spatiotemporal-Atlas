suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

# Step 1 read-only audit and root-candidate evaluation.
# Does not modify or save the Seurat object.

started_at <- Sys.time()
cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(cwd) == "scripts") {
  step_dir <- normalizePath(file.path(cwd, ".."), winslash = "/", mustWork = TRUE)
} else if (file.exists(file.path(cwd, "scripts", "_naming.R"))) {
  step_dir <- cwd
} else {
  candidate <- file.path(cwd, "06_恶性细胞拟时序")
  if (!dir.exists(candidate)) candidate <- file.path(cwd, "06_鎭舵€х粏鑳炴嫙鏃跺簭")
  step_dir <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
}
setwd(step_dir)

dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/source_data", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)

source(file.path("scripts", "_naming.R"))
write_subtype_naming_mapping(file.path("tables", "subtype_naming_mapping.csv"))

log_file <- file.path("logs", "01_audit.log")
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = " "))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}
if (file.exists(log_file)) file.remove(log_file)

expected_cells <- 28213L
short_colors <- setNames(subtype_naming_mapping$color, subtype_naming_mapping$abbreviation)
short_levels <- subtype_naming_mapping$abbreviation

obj_path <- file.path(
  "..",
  "05_恶性细胞分亚群与Neftel对照",
  "outputs",
  "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"
)
if (!file.exists(obj_path)) {
  obj_path <- file.path(
    "..",
    "05_鎭舵€х粏鑳炲垎浜氱兢涓嶯eftel瀵圭収",
    "outputs",
    "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"
  )
}
obj_path <- normalizePath(obj_path, winslash = "/", mustWork = TRUE)

log_msg("Input object:", obj_path)
log_msg("Naming source of truth:", file.path("scripts", "_naming.R"))
log_msg("Gene blacklist:", file.path("data", "gene_blacklist.txt"))

load_start <- Sys.time()
obj <- qs2::qs_read(obj_path)
load_end <- Sys.time()
log_msg("Object loaded in seconds:", round(as.numeric(difftime(load_end, load_start, units = "secs")), 2))

n_cells <- ncol(obj)
n_genes <- nrow(obj)
log_msg("Object dimensions genes x cells:", n_genes, "x", n_cells)
if (n_cells != expected_cells) {
  log_msg("WARNING: expected cell count", expected_cells, "but observed", n_cells)
} else {
  log_msg("Cell count sanity passed:", expected_cells)
}

md <- obj@meta.data |>
  tibble::rownames_to_column("cell_id")

required_cols <- c("subtype_k4", "subtype_label_final")
missing_required <- setdiff(required_cols, colnames(md))
if (length(missing_required) > 0) {
  stop("Missing required metadata columns: ", paste(missing_required, collapse = ", "))
}
md <- recode_subtype(md)

reduction_names <- names(obj@reductions)
reduction_audit <- lapply(reduction_names, function(red) {
  emb <- Seurat::Embeddings(obj, reduction = red)
  tibble::tibble(
    section = "reduction",
    item = red,
    field = c("n_cells", "n_dims"),
    value = as.character(c(nrow(emb), ncol(emb)))
  )
}) |>
  bind_rows()

metadata_audit <- tibble::tibble(
  section = "metadata_colname",
  item = colnames(obj@meta.data),
  field = "present",
  value = "TRUE"
)

object_audit <- tibble::tibble(
  section = "object",
  item = c("n_genes", "n_cells", "default_assay"),
  field = "value",
  value = c(as.character(n_genes), as.character(n_cells), Seurat::DefaultAssay(obj))
)

subtype_cross <- md |>
  count(subtype_k4, subtype_label_source_metadata, subtype_label_original, subtype_label_final, subtype_short, name = "n_cells") |>
  mutate(
    section = "subtype_cross_table",
    item = paste(subtype_k4, subtype_label_source_metadata, subtype_label_final, subtype_short, sep = " | "),
    field = "n_cells",
    value = as.character(n_cells)
  ) |>
  select(section, item, field, value)

key_patterns <- c("^MP0[1-6]$", "^MP0[1-6]_UCell$", "^AMS_", "_UCell$", "^Phase$", "^Pt_number$")
key_cols <- unique(unlist(lapply(key_patterns, function(p) grep(p, colnames(md), value = TRUE))))
na_audit <- tibble::tibble(
  section = "key_field_NA",
  item = key_cols,
  field = "n_NA",
  value = as.character(vapply(md[key_cols], function(x) sum(is.na(x)), integer(1)))
)

audit <- bind_rows(object_audit, metadata_audit, reduction_audit, subtype_cross, na_audit)
readr::write_csv(audit, file.path("tables", "拟时序_metadata审计.csv"))
log_msg("Wrote metadata audit:", file.path("tables", "拟时序_metadata审计.csv"))

gini_coefficient <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0 || sum(x) == 0) return(NA_real_)
  x <- sort(x)
  n <- length(x)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}

if ("Pt_number" %in% colnames(md)) {
  patient_subtype <- md |>
    count(subtype_k4, subtype_label_source_metadata, subtype_label_original, subtype_label_final, subtype_short, Pt_number, name = "n_cells") |>
    group_by(subtype_k4, subtype_label_final, subtype_short) |>
    mutate(
      subtype_total_cells = sum(n_cells),
      patient_fraction_within_subtype = n_cells / subtype_total_cells,
      n_patients_with_cells = sum(n_cells > 0),
      patient_contribution_gini = gini_coefficient(n_cells)
    ) |>
    ungroup() |>
    arrange(subtype_k4, desc(n_cells), Pt_number)
  readr::write_csv(patient_subtype, file.path("tables", "拟时序_per_patient_subtype分布.csv"))
  log_msg("Wrote per-patient subtype distribution:", file.path("tables", "拟时序_per_patient_subtype分布.csv"))
}

get_assay_matrix <- function(object, assay = "RNA", slot = "data") {
  if (packageVersion("SeuratObject") >= "5.0.0") {
    return(SeuratObject::LayerData(object, assay = assay, layer = slot))
  }
  Seurat::GetAssayData(object, assay = assay, slot = slot)
}

rna_assay <- if ("RNA" %in% names(obj@assays)) "RNA" else Seurat::DefaultAssay(obj)
data_mat <- get_assay_matrix(obj, assay = rna_assay, slot = "data")
genes <- c("SOX2", "SOX4", "SOX9", "SOX11", "GAP43")
present_genes <- intersect(genes, rownames(data_mat))
missing_genes <- setdiff(genes, present_genes)
if (length(missing_genes) > 0) log_msg("WARNING: missing genes:", paste(missing_genes, collapse = ", "))

expr_df <- as.matrix(data_mat[present_genes, , drop = FALSE]) |>
  t() |>
  as.data.frame(check.names = FALSE) |>
  tibble::rownames_to_column("cell_id") |>
  pivot_longer(-cell_id, names_to = "metric_name", values_to = "value")

mp_raw_cols <- intersect(c(paste0("MP0", 1:6), paste0("MP0", 1:6, "_UCell")), colnames(md))
score_cols <- intersect(c(mp_raw_cols, "AMS_NPC1", "AMS_NPC2", "AMS_OPC", "AMS_MES1", "AMS_MES2", "AMS_AC"), colnames(md))
score_df <- md |>
  select(cell_id, all_of(score_cols)) |>
  pivot_longer(-cell_id, names_to = "metric_name", values_to = "value") |>
  mutate(metric_name = sub("_UCell$", "", metric_name))

phase_df <- if ("Phase" %in% colnames(md)) {
  tidyr::expand_grid(
    cell_id = md$cell_id,
    metric_name = c("Phase_G1", "Phase_S", "Phase_G2M")
  ) |>
    left_join(md |> select(cell_id, Phase), by = "cell_id") |>
    mutate(value = as.numeric(metric_name == paste0("Phase_", Phase))) |>
    select(cell_id, metric_name, value)
} else {
  log_msg("WARNING: Phase column missing.")
  tibble::tibble(cell_id = character(), metric_name = character(), value = numeric())
}

install_cytotrace2_if_needed <- function() {
  if (requireNamespace("CytoTRACE2", quietly = TRUE)) return(TRUE)
  local_src <- file.path("data", "cytotrace2_src", "cytotrace2-main", "cytotrace2_r")
  if (dir.exists(local_src)) {
    log_msg("CytoTRACE2 not installed. Installing from local source:", local_src)
    install.packages(local_src, repos = NULL, type = "source")
    if (requireNamespace("CytoTRACE2", quietly = TRUE)) return(TRUE)
  }
  log_msg("CytoTRACE2 not installed. Installing from digitalcytometry/cytotrace2, subdir cytotrace2_r.")
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
  remotes::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r", upgrade = "never")
  requireNamespace("CytoTRACE2", quietly = TRUE)
}

run_cytotrace2 <- function(object) {
  cache_paths <- c(
    file.path("outputs", "CytoTRACE2_per_cell.csv"),
    file.path("tables", "拟时序_CytoTRACE2_per_cell.csv")
  )
  for (cache_path in cache_paths) {
    if (file.exists(cache_path)) {
      log_msg("Using cached CytoTRACE2 per-cell scores:", cache_path)
      cached <- readr::read_csv(cache_path, show_col_types = FALSE)
      if (all(c("cell_id", "CytoTRACE2_score") %in% colnames(cached))) {
        readr::write_csv(cached, file.path("outputs", "CytoTRACE2_per_cell.csv"))
        readr::write_csv(cached, file.path("tables", "拟时序_CytoTRACE2_per_cell.csv"))
        return(cached)
      }
    }
  }
  ok <- install_cytotrace2_if_needed()
  if (!ok) stop("CytoTRACE2 installation failed.")
  cytotrace_start <- Sys.time()
  result <- CytoTRACE2::cytotrace2(
    object,
    species = "human",
    is_seurat = TRUE,
    slot_type = "counts",
    batch_size = 10000,
    smooth_batch_size = 1000,
    parallelize_models = TRUE,
    parallelize_smoothing = TRUE,
    ncores = max(1L, parallel::detectCores(logical = FALSE) - 1L),
    seed = 14
  )
  cytotrace_end <- Sys.time()
  log_msg("CytoTRACE2 run seconds:", round(as.numeric(difftime(cytotrace_end, cytotrace_start, units = "secs")), 2))
  ct_md <- result@meta.data |> tibble::rownames_to_column("cell_id")
  score_col <- intersect(c("CytoTRACE2_Score", "CytoTRACE2_score", "cytotrace2_score"), colnames(ct_md))
  if (length(score_col) == 0) {
    stop("CytoTRACE2 finished but no score column was found. Metadata columns: ", paste(colnames(ct_md), collapse = ", "))
  }
  out <- ct_md |> transmute(cell_id, CytoTRACE2_score = .data[[score_col[1]]])
  readr::write_csv(out, file.path("outputs", "CytoTRACE2_per_cell.csv"))
  readr::write_csv(out, file.path("tables", "拟时序_CytoTRACE2_per_cell.csv"))
  log_msg("Wrote CytoTRACE2 per-cell scores:", file.path("outputs", "CytoTRACE2_per_cell.csv"))
  out
}

cytotrace_df <- run_cytotrace2(obj) |>
  transmute(cell_id, metric_name = "CytoTRACE2_score", value = CytoTRACE2_score)

metric_df <- bind_rows(cytotrace_df, phase_df, expr_df, score_df) |>
  left_join(
    md |> select(cell_id, subtype_k4, subtype_label_source_metadata, subtype_label_original, subtype_label_final, subtype_short),
    by = "cell_id"
  ) |>
  rename(subtype_label_manuscript = subtype_label_final)

read_gene_blacklist <- function(path = file.path("data", "gene_blacklist.txt")) {
  patterns <- readLines(path, warn = FALSE)
  patterns[!grepl("^\\s*(#|$)", patterns)]
}

write_cleaned_driver_candidates <- function() {
  patterns <- read_gene_blacklist()
  de_path <- file.path("..", "05_恶性细胞分亚群与Neftel对照", "tables", "10e_one_vs_rest_FindAllMarkers_wilcoxon_downsampled.csv")
  if (!file.exists(de_path)) {
    de_path <- file.path("..", "05_鎭舵€х粏鑳炲垎浜氱兢涓嶯eftel瀵圭収", "tables", "10e_one_vs_rest_FindAllMarkers_wilcoxon_downsampled.csv")
  }
  if (!file.exists(de_path)) {
    log_msg("WARNING: raw DE table not found for cleaned candidate export.")
    return(invisible(NULL))
  }
  de <- readr::read_csv(de_path, show_col_types = FALSE)
  blacklisted <- Reduce(`|`, lapply(patterns, function(p) grepl(p, de$gene)))
  de_flagged <- de |> mutate(is_blacklisted = blacklisted)
  cleaned <- de_flagged |>
    filter(!is_blacklisted) |>
    recode_subtype()
  readr::write_csv(cleaned, file.path("tables", "trajectory_driver_gene_candidates_cleaned.csv"))
  readr::write_csv(de_flagged |> filter(is_blacklisted), file.path("tables", "trajectory_driver_gene_candidates_blacklisted.csv"))
  log_msg("Wrote cleaned driver gene candidates:", nrow(cleaned), "rows; blacklisted:", sum(blacklisted))
}
write_cleaned_driver_candidates()

root_summary <- metric_df |>
  group_by(subtype_k4, subtype_label_source_metadata, subtype_label_original, subtype_label_manuscript, subtype_short, metric_name) |>
  summarise(
    mean = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    pct_pos = mean(value > 0, na.rm = TRUE),
    n_cells = sum(!is.na(value)),
    .groups = "drop"
  ) |>
  arrange(subtype_k4, metric_name)
readr::write_csv(root_summary, file.path("tables", "拟时序_root候选评估.csv"))
log_msg("Wrote root candidate summary:", file.path("tables", "拟时序_root候选评估.csv"))

p_to_text <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 1e-4, formatC(p, format = "e", digits = 2), signif(p, 3)))
}
p_to_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

kw_one <- function(df, value_col, facet_col = NULL) {
  value_sym <- rlang::sym(value_col)
  if (is.null(facet_col)) {
    return(tibble::tibble(
      metric = "all",
      test = "Kruskal-Wallis",
      p_value = kruskal.test(df[[value_col]] ~ df$subtype_short)$p.value
    ))
  }
  df |>
    group_by(.data[[facet_col]]) |>
    summarise(
      test = "Kruskal-Wallis",
      p_value = kruskal.test(!!value_sym ~ subtype_short)$p.value,
      .groups = "drop"
    ) |>
    rename(metric = all_of(facet_col))
}

dunn_vs_root <- function(df, value_col, facet_col = NULL, root = "NPC-P") {
  if (!requireNamespace("dunn.test", quietly = TRUE)) {
    install.packages("dunn.test", repos = "https://cloud.r-project.org")
  }
  run_one <- function(x, metric = "all") {
    x <- x |> filter(!is.na(.data[[value_col]]), !is.na(subtype_short))
    dt <- dunn.test::dunn.test(x[[value_col]], x$subtype_short, method = "bh", kw = FALSE, list = TRUE)
    tibble::tibble(
      metric = metric,
      comparison = dt$comparisons,
      z = dt$Z,
      p_unadjusted = dt$P,
      p_adj = dt$P.adjusted
    ) |>
      filter(grepl(root, comparison)) |>
      mutate(
        compared_to_root = sub(paste0(root, " - | - ", root), "", comparison),
        star = p_to_star(p_adj)
      )
  }
  if (is.null(facet_col)) return(run_one(df))
  split(df, df[[facet_col]]) |>
    lapply(function(x) run_one(x, unique(x[[facet_col]])[1])) |>
    bind_rows()
}

plot_df <- metric_df |>
  mutate(subtype_short = factor(subtype_short, levels = short_levels))

cyto_source <- plot_df |>
  filter(metric_name == "CytoTRACE2_score") |>
  select(cell_id, subtype_k4, subtype_label_source_metadata, subtype_label_original, subtype_label_manuscript, subtype_short, CytoTRACE2_score = value)
readr::write_csv(cyto_source, file.path("figures/source_data", "01_CytoTRACE2_by_subtype.csv"))

phase_source <- plot_df |>
  filter(grepl("^Phase_", metric_name)) |>
  mutate(Phase = sub("^Phase_", "", metric_name)) |>
  group_by(subtype_k4, subtype_label_source_metadata, subtype_label_original, subtype_label_manuscript, subtype_short, Phase) |>
  summarise(n_cells = sum(value > 0, na.rm = TRUE), .groups = "drop_last") |>
  mutate(prop = n_cells / sum(n_cells)) |>
  ungroup()
readr::write_csv(phase_source, file.path("figures/source_data", "01_cell_cycle_phase_composition.csv"))

sox_source <- plot_df |>
  filter(metric_name %in% genes) |>
  select(cell_id, subtype_k4, subtype_label_source_metadata, subtype_label_original, subtype_label_manuscript, subtype_short, gene = metric_name, expression = value)
readr::write_csv(sox_source, file.path("figures/source_data", "01_SOX_GAP43_expression.csv"))
readr::write_csv(sox_source, file.path("figures/source_data", "01_SOX_family_expression.csv"))

neftel_source <- plot_df |>
  filter(metric_name %in% c("AMS_NPC1", "AMS_NPC2", "AMS_OPC", "AMS_MES1", "AMS_MES2", "AMS_AC")) |>
  select(cell_id, subtype_k4, subtype_label_source_metadata, subtype_label_original, subtype_label_manuscript, subtype_short, module = metric_name, score = value)
readr::write_csv(neftel_source, file.path("figures/source_data", "01_Neftel_6state_AMS.csv"))

cyto_kw <- kw_one(cyto_source, "CytoTRACE2_score")
cyto_dunn <- dunn_vs_root(cyto_source, "CytoTRACE2_score")
sox_kw <- kw_one(sox_source, "expression", "gene")
sox_dunn <- dunn_vs_root(sox_source, "expression", "gene")
neftel_kw <- kw_one(neftel_source, "score", "module")
neftel_dunn <- dunn_vs_root(neftel_source, "score", "module")
phase_chisq <- chisq.test(xtabs(n_cells ~ subtype_short + Phase, data = phase_source))

stats_tables <- list(
  CytoTRACE2_KW = cyto_kw,
  CytoTRACE2_Dunn_vs_NPC_P = cyto_dunn,
  SOX_GAP43_KW = sox_kw,
  SOX_GAP43_Dunn_vs_NPC_P = sox_dunn,
  Neftel_KW = neftel_kw,
  Neftel_Dunn_vs_NPC_P = neftel_dunn,
  Phase_chisq = tibble::tibble(test = "Chi-square", statistic = unname(phase_chisq$statistic), df = unname(phase_chisq$parameter), p_value = phase_chisq$p.value)
)
readr::write_csv(bind_rows(stats_tables, .id = "table"), file.path("tables", "拟时序_root候选统计检验.csv"))
log_msg("Wrote statistical tests:", file.path("tables", "拟时序_root候选统计检验.csv"))

make_star_ann <- function(source, value_col, stat_df, facet_col = NULL) {
  x_map <- tibble::tibble(compared_to_root = short_levels[-1], x = 2:4)
  if (is.null(facet_col)) {
    ymax <- max(source[[value_col]], na.rm = TRUE)
    ymin <- min(source[[value_col]], na.rm = TRUE)
    span <- max(ymax - ymin, 0.1)
    return(stat_df |>
      left_join(x_map, by = "compared_to_root") |>
      mutate(y = ymax + span * 0.08 * row_number(), label = star))
  }
  source_range <- source |>
    group_by(.data[[facet_col]]) |>
    summarise(ymax = max(.data[[value_col]], na.rm = TRUE), ymin = min(.data[[value_col]], na.rm = TRUE), .groups = "drop") |>
    rename(metric = all_of(facet_col)) |>
    mutate(span = pmax(ymax - ymin, 0.1))
  stat_df |>
    left_join(x_map, by = "compared_to_root") |>
    left_join(source_range, by = "metric") |>
    group_by(metric) |>
    mutate(y = ymax + span * 0.08 * row_number(), label = star) |>
    ungroup()
}

cyto_ann <- make_star_ann(cyto_source, "CytoTRACE2_score", cyto_dunn)
sox_ann <- make_star_ann(sox_source, "expression", sox_dunn, "gene")
neftel_ann <- make_star_ann(neftel_source, "score", neftel_dunn, "module")

sox_labels <- sox_kw |> mutate(gene_label = paste0(metric, "\nKW p=", p_to_text(p_value))) |> select(gene = metric, gene_label)
neftel_labels <- neftel_kw |> mutate(module_label = paste0(metric, "\nKW p=", p_to_text(p_value))) |> select(module = metric, module_label)
sox_source <- sox_source |> left_join(sox_labels, by = "gene")
sox_ann <- sox_ann |> left_join(sox_labels, by = c("metric" = "gene"))
neftel_source <- neftel_source |> left_join(neftel_labels, by = "module")
neftel_ann <- neftel_ann |> left_join(neftel_labels, by = c("metric" = "module"))

base_theme <- theme_bw(base_size = 7) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    strip.text = element_text(size = 6.5),
    plot.title = element_text(size = 8),
    plot.subtitle = element_text(size = 6.5)
  )

p1 <- ggplot(cyto_source, aes(x = subtype_short, y = CytoTRACE2_score, fill = subtype_short)) +
  geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.15, linewidth = 0.2) +
  geom_text(data = cyto_ann, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
  scale_fill_manual(values = short_colors, guide = "none") +
  labs(x = NULL, y = "CytoTRACE2 score", title = "A. CytoTRACE2", subtitle = paste0("KW p=", p_to_text(cyto_kw$p_value), "; stars: Dunn BH vs NPC-P")) +
  base_theme

p2 <- ggplot(phase_source, aes(x = subtype_short, y = prop, fill = Phase)) +
  geom_col(width = 0.75, linewidth = 0.2, color = "white") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Cell fraction", title = "B. Cell cycle phase", subtitle = paste0("Chi-square p=", p_to_text(phase_chisq$p.value))) +
  base_theme

p3 <- ggplot(sox_source, aes(x = subtype_short, y = expression, fill = subtype_short)) +
  geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.15, linewidth = 0.2) +
  geom_text(data = sox_ann, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 1.8) +
  facet_wrap(~gene_label, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = short_colors, guide = "none") +
  labs(x = NULL, y = "RNA normalized expression", title = "C. SOX family and GAP43", subtitle = "Stars: Dunn BH vs NPC-P") +
  base_theme

p4 <- ggplot(neftel_source, aes(x = subtype_short, y = score, fill = subtype_short)) +
  geom_violin(scale = "width", linewidth = 0.2, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.size = 0.15, linewidth = 0.2) +
  geom_text(data = neftel_ann, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 1.8) +
  facet_wrap(~module_label, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = short_colors, guide = "none") +
  labs(x = NULL, y = "AMS score", title = "D. Neftel 6-state AMS", subtitle = "Stars: Dunn BH vs NPC-P") +
  base_theme

pdf(file.path("figures", "01_root候选评估.pdf"), width = 8.4, height = 7.6, useDingbats = FALSE)
if (requireNamespace("patchwork", quietly = TRUE)) {
  print((p1 | p2) / (p3 / p4))
} else if (requireNamespace("gridExtra", quietly = TRUE)) {
  gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
} else {
  print(p1); print(p2); print(p3); print(p4)
}
dev.off()
log_msg("Wrote figure:", file.path("figures", "01_root候选评估.pdf"))

session_info <- capture.output(sessionInfo())
cat(session_info, sep = "\n", file = file.path("logs", "01_session_info.txt"))
log_msg("Session info:", file.path("logs", "01_session_info.txt"))
log_msg("Total elapsed seconds:", round(as.numeric(difftime(Sys.time(), started_at, units = "secs")), 2))
log_msg("STOP: Cleanup Step 1 finished. Do not proceed to Slingshot or Monocle3 automatically.")
