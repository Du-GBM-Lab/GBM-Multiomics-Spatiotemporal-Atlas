suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(qs2)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

if (!requireNamespace("tidyr", quietly = TRUE)) {
  stop("Required package not installed: tidyr", call. = FALSE)
}

config <- list(
  input_rds = normalizePath(file.path("..", "GBM.RNA.integrated.24.rds"), winslash = "\\", mustWork = TRUE),
  out_dir = normalizePath(".", winslash = "\\", mustWork = TRUE),
  sample_col = "Pt_number",
  annotation_col = "anno_ident",
  assay = "RNA",
  mt_pattern = "^MT-",
  mt_hard_cap = 20,
  nmads = 3,
  run_doubletfinder = TRUE,
  run_pk_sweep = FALSE,
  default_pk = 0.09,
  pN = 0.25,
  pcs = 1:30,
  sample_resolution = 0.5,
  min_cells_for_doubletfinder = 500,
  keep_doubletfinder_skipped = TRUE,
  max_doublet_rate = 0.20,
  random_seed = 20260506
)

dir.create(file.path(config$out_dir, "outputs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

set.seed(config$random_seed)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Required package not installed: ", pkg, "\n",
      "For DoubletFinder, install it before running this script, for example:\n",
      "  remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')\n",
      call. = FALSE
    )
  }
}

get_df_function <- function(candidates) {
  ns <- asNamespace("DoubletFinder")
  for (fn in candidates) {
    if (exists(fn, envir = ns, inherits = FALSE)) {
      return(get(fn, envir = ns))
    }
  }
  stop("Cannot find compatible DoubletFinder function: ", paste(candidates, collapse = " / "), call. = FALSE)
}

join_rna_layers_if_needed <- function(obj, assay = "RNA") {
  obj <- UpdateSeuratObject(obj)
  DefaultAssay(obj) <- assay
  if (inherits(obj[[assay]], "Assay5")) {
    obj[[assay]] <- JoinLayers(obj[[assay]])
  }
  obj
}

estimate_10x_doublet_rate <- function(n_cells, max_rate = 0.20) {
  # Approximate 10x expected doublet rate: about 0.8% per 1,000 recovered cells.
  min(max_rate, 0.008 * n_cells / 1000)
}

robust_bounds <- function(x, nmads = 3, lower = TRUE, upper = TRUE) {
  med <- median(x, na.rm = TRUE)
  mad_value <- mad(x, constant = 1, na.rm = TRUE)
  if (is.na(mad_value) || mad_value == 0) {
    mad_value <- stats::IQR(x, na.rm = TRUE) / 1.349
  }
  if (is.na(mad_value) || mad_value == 0) {
    mad_value <- 0
  }
  list(
    median = med,
    mad = mad_value,
    lower = if (lower) med - nmads * mad_value else -Inf,
    upper = if (upper) med + nmads * mad_value else Inf
  )
}

add_samplewise_qc_flags <- function(obj, config) {
  md <- obj@meta.data
  stopifnot(config$sample_col %in% colnames(md))

  needed_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
  missing_metrics <- setdiff(needed_metrics, colnames(md))
  if (length(missing_metrics) > 0) {
    stop("Missing QC metrics: ", paste(missing_metrics, collapse = ", "), call. = FALSE)
  }

  md$qc_fail_nCount_low <- FALSE
  md$qc_fail_nCount_high <- FALSE
  md$qc_fail_nFeature_low <- FALSE
  md$qc_fail_nFeature_high <- FALSE
  md$qc_fail_percent_mt_mad <- FALSE
  md$qc_fail_percent_mt20 <- md$percent.mt >= config$mt_hard_cap

  threshold_rows <- list()
  samples <- sort(unique(as.character(md[[config$sample_col]])))

  for (sample_id in samples) {
    idx <- which(as.character(md[[config$sample_col]]) == sample_id)
    log_n_count <- log1p(md$nCount_RNA[idx])
    log_n_feature <- log1p(md$nFeature_RNA[idx])
    b_count <- robust_bounds(log_n_count, nmads = config$nmads, lower = TRUE, upper = TRUE)
    b_feature <- robust_bounds(log_n_feature, nmads = config$nmads, lower = TRUE, upper = TRUE)
    b_mt <- robust_bounds(md$percent.mt[idx], nmads = config$nmads, lower = FALSE, upper = TRUE)

    md$qc_fail_nCount_low[idx] <- log_n_count < b_count$lower
    md$qc_fail_nCount_high[idx] <- log_n_count > b_count$upper
    md$qc_fail_nFeature_low[idx] <- log_n_feature < b_feature$lower
    md$qc_fail_nFeature_high[idx] <- log_n_feature > b_feature$upper
    md$qc_fail_percent_mt_mad[idx] <- md$percent.mt[idx] > b_mt$upper

    threshold_rows[[sample_id]] <- data.frame(
      sample = sample_id,
      n_cells_before = length(idx),
      log1p_nCount_lower = b_count$lower,
      log1p_nCount_upper = b_count$upper,
      nCount_lower = max(0, expm1(b_count$lower)),
      nCount_upper = expm1(b_count$upper),
      log1p_nFeature_lower = b_feature$lower,
      log1p_nFeature_upper = b_feature$upper,
      nFeature_lower = max(0, expm1(b_feature$lower)),
      nFeature_upper = expm1(b_feature$upper),
      percent_mt_mad_upper = b_mt$upper,
      percent_mt_hard_upper = config$mt_hard_cap,
      stringsAsFactors = FALSE
    )
  }

  md$qc_pass_mad <- !(
    md$qc_fail_nCount_low |
      md$qc_fail_nCount_high |
      md$qc_fail_nFeature_low |
      md$qc_fail_nFeature_high |
      md$qc_fail_percent_mt_mad
  )
  md$qc_pass_mt20 <- !md$qc_fail_percent_mt20
  md$qc_pass_pre_doublet <- md$qc_pass_mad & md$qc_pass_mt20
  obj@meta.data <- md

  thresholds <- bind_rows(threshold_rows)
  list(obj = obj, thresholds = thresholds)
}

plot_qc_violin <- function(md, out_file) {
  plot_df <- md |>
    select(sample, qc_stage, nCount_RNA, nFeature_RNA, percent.mt) |>
    tidyr::pivot_longer(
      cols = c(nCount_RNA, nFeature_RNA, percent.mt),
      names_to = "metric",
      values_to = "value"
    )

  p <- ggplot(plot_df, aes(x = sample, y = value, fill = qc_stage)) +
    geom_violin(scale = "width", linewidth = 0.15, trim = TRUE) +
    facet_wrap(~metric, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c(Before_QC = "#8A8F98", After_QC = "#2C7FB8")) +
    labs(x = NULL, y = NULL, fill = NULL) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "top"
    )
  ggsave(out_file, p, width = 10, height = 8, units = "in", dpi = 300)
}

plot_bar <- function(df, x, y, fill, out_file, y_label = "Cells") {
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = .data[[fill]])) +
    geom_col(width = 0.78) +
    labs(x = NULL, y = y_label, fill = NULL) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  ggsave(out_file, p, width = 9, height = 5.5, units = "in", dpi = 300)
}

choose_pk <- function(sample_obj, pcs, sample_id, config) {
  if (!config$run_pk_sweep) {
    return(config$default_pk)
  }
  param_sweep <- get_df_function(c("paramSweep", "paramSweep_v3"))
  summarize_sweep <- get_df_function(c("summarizeSweep"))
  find_pk <- get_df_function(c("find.pK"))

  pk_table <- tryCatch({
    sweep_res <- param_sweep(sample_obj, PCs = pcs, sct = FALSE)
    sweep_stats <- summarize_sweep(sweep_res, GT = FALSE)
    find_pk(sweep_stats)
  }, error = function(e) {
    warning("DoubletFinder pK sweep failed for ", sample_id, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(pk_table) || nrow(pk_table) == 0) {
    return(config$default_pk)
  }
  write.csv(
    pk_table,
    file.path(config$out_dir, "tables", paste0("DoubletFinder_pK_sweep_", sample_id, ".csv")),
    row.names = FALSE
  )
  pk_table$pK_numeric <- as.numeric(as.character(pk_table$pK))
  selected_pk <- pk_table$pK_numeric[which.max(pk_table$BCmetric)]
  if (is.na(selected_pk)) {
    selected_pk <- config$default_pk
  }
  selected_pk
}

run_doubletfinder_one_sample <- function(obj, sample_id, config) {
  require_pkg("DoubletFinder")
  set.seed(config$random_seed + sum(utf8ToInt(sample_id)))

  cells <- colnames(obj)[as.character(obj@meta.data[[config$sample_col]]) == sample_id]
  sample_obj <- subset(obj, cells = cells)
  sample_obj <- subset(sample_obj, cells = colnames(sample_obj)[sample_obj$qc_pass_pre_doublet])

  n_cells <- ncol(sample_obj)
  if (n_cells < config$min_cells_for_doubletfinder) {
    return(data.frame(
      cell = colnames(sample_obj),
      sample = sample_id,
      doubletfinder_class = "Skipped_low_cells",
      doubletfinder_pANN = NA_real_,
      doubletfinder_pK = NA_real_,
      doubletfinder_nExp = 0,
      doubletfinder_nExp_adj = 0,
      doubletfinder_rate = 0,
      doubletfinder_note = paste0("Skipped: fewer than ", config$min_cells_for_doubletfinder, " cells"),
      stringsAsFactors = FALSE
    ))
  }

  DefaultAssay(sample_obj) <- config$assay
  sample_obj <- NormalizeData(sample_obj, verbose = FALSE)
  sample_obj <- FindVariableFeatures(sample_obj, nfeatures = 2000, verbose = FALSE)
  sample_obj <- ScaleData(sample_obj, features = VariableFeatures(sample_obj), verbose = FALSE)
  sample_obj <- RunPCA(sample_obj, features = VariableFeatures(sample_obj), npcs = max(config$pcs), verbose = FALSE)
  pcs_use <- config$pcs[config$pcs <= ncol(Embeddings(sample_obj, "pca"))]
  sample_obj <- FindNeighbors(sample_obj, dims = pcs_use, verbose = FALSE)
  sample_obj <- FindClusters(sample_obj, resolution = config$sample_resolution, verbose = FALSE)

  pK <- choose_pk(sample_obj, pcs_use, sample_id, config)
  doublet_rate <- estimate_10x_doublet_rate(n_cells, config$max_doublet_rate)
  nExp_poi <- round(doublet_rate * n_cells)

  model_homotypic <- get_df_function(c("modelHomotypic"))
  sample_clusters <- as.character(Idents(sample_obj))
  homotypic_prop <- model_homotypic(sample_clusters)
  nExp_poi_adj <- round(nExp_poi * (1 - homotypic_prop))

  doublet_finder <- get_df_function(c("doubletFinder", "doubletFinder_v3"))
  before_cols <- colnames(sample_obj@meta.data)
  sample_obj <- doublet_finder(
    sample_obj,
    PCs = pcs_use,
    pN = config$pN,
    pK = pK,
    nExp = nExp_poi_adj,
    reuse.pANN = NULL,
    sct = FALSE
  )

  after_cols <- colnames(sample_obj@meta.data)
  new_cols <- setdiff(after_cols, before_cols)
  class_cols <- grep("^DF\\.classifications", new_cols, value = TRUE)
  pann_cols <- grep("^pANN", new_cols, value = TRUE)
  if (length(class_cols) == 0) {
    class_cols <- grep("^DF\\.classifications", after_cols, value = TRUE)
  }
  if (length(pann_cols) == 0) {
    pann_cols <- grep("^pANN", after_cols, value = TRUE)
  }
  class_col <- tail(class_cols, 1)
  pann_col <- tail(pann_cols, 1)
  if (length(class_col) == 0 || length(pann_col) == 0) {
    stop("DoubletFinder completed but classification or pANN column was not found for sample: ", sample_id, call. = FALSE)
  }

  data.frame(
    cell = colnames(sample_obj),
    sample = sample_id,
    doubletfinder_class = as.character(sample_obj@meta.data[[class_col]]),
    doubletfinder_pANN = as.numeric(sample_obj@meta.data[[pann_col]]),
    doubletfinder_pK = pK,
    doubletfinder_nExp = nExp_poi,
    doubletfinder_nExp_adj = nExp_poi_adj,
    doubletfinder_rate = doublet_rate,
    doubletfinder_note = "DoubletFinder completed",
    stringsAsFactors = FALSE
  )
}

msg("Loading object:", config$input_rds)
obj <- readRDS(config$input_rds)
obj <- join_rna_layers_if_needed(obj, config$assay)

if (!"percent.mt" %in% colnames(obj@meta.data)) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = config$mt_pattern, assay = config$assay)
}

qc <- add_samplewise_qc_flags(obj, config)
obj <- qc$obj

write.csv(qc$thresholds, file.path(config$out_dir, "tables", "samplewise_MAD_QC_thresholds.csv"), row.names = FALSE)

md <- obj@meta.data
md$sample <- as.character(md[[config$sample_col]])

qc_counts <- md |>
  mutate(
    qc_status = case_when(
      qc_pass_pre_doublet ~ "Pass_pre_doublet",
      qc_fail_percent_mt20 ~ "Fail_MT20",
      !qc_pass_mad ~ "Fail_MAD",
      TRUE ~ "Fail_other"
    )
  ) |>
  count(sample, qc_status, name = "n_cells") |>
  group_by(sample) |>
  mutate(percent = n_cells / sum(n_cells) * 100) |>
  ungroup()
write.csv(qc_counts, file.path(config$out_dir, "tables", "QC_filtering_summary_by_sample.csv"), row.names = FALSE)

plot_bar(
  qc_counts,
  x = "sample",
  y = "n_cells",
  fill = "qc_status",
  out_file = file.path(config$out_dir, "figures", "QC_filtering_summary_by_sample.pdf")
)

qc_plot_md <- bind_rows(
  md |> mutate(qc_stage = "Before_QC"),
  md |> filter(qc_pass_pre_doublet) |> mutate(qc_stage = "After_QC")
)
plot_qc_violin(qc_plot_md, file.path(config$out_dir, "figures", "QC_metrics_before_after_samplewise.pdf"))

if (config$run_doubletfinder) {
  require_pkg("DoubletFinder")
  samples <- sort(unique(as.character(obj@meta.data[[config$sample_col]])))
  df_results <- list()
  for (sample_id in samples) {
    msg("Running DoubletFinder for", sample_id)
    df_results[[sample_id]] <- run_doubletfinder_one_sample(obj, sample_id, config)
  }
  doublet_calls <- bind_rows(df_results)
  write.csv(doublet_calls, file.path(config$out_dir, "tables", "DoubletFinder_calls_by_cell.csv"), row.names = FALSE)

  obj$doubletfinder_class <- NA_character_
  obj$doubletfinder_pANN <- NA_real_
  obj$doubletfinder_pK <- NA_real_
  obj$doubletfinder_nExp <- NA_integer_
  obj$doubletfinder_nExp_adj <- NA_integer_
  obj$doubletfinder_rate <- NA_real_

  rownames(doublet_calls) <- doublet_calls$cell
  shared_cells <- intersect(colnames(obj), doublet_calls$cell)
  obj@meta.data[shared_cells, "doubletfinder_class"] <- doublet_calls[shared_cells, "doubletfinder_class"]
  obj@meta.data[shared_cells, "doubletfinder_pANN"] <- doublet_calls[shared_cells, "doubletfinder_pANN"]
  obj@meta.data[shared_cells, "doubletfinder_pK"] <- doublet_calls[shared_cells, "doubletfinder_pK"]
  obj@meta.data[shared_cells, "doubletfinder_nExp"] <- doublet_calls[shared_cells, "doubletfinder_nExp"]
  obj@meta.data[shared_cells, "doubletfinder_nExp_adj"] <- doublet_calls[shared_cells, "doubletfinder_nExp_adj"]
  obj@meta.data[shared_cells, "doubletfinder_rate"] <- doublet_calls[shared_cells, "doubletfinder_rate"]

  is_singlet <- replace(obj$doubletfinder_class == "Singlet", is.na(obj$doubletfinder_class), FALSE)
  is_skipped <- replace(obj$doubletfinder_class == "Skipped_low_cells", is.na(obj$doubletfinder_class), FALSE)
  obj$qc_pass_final <- obj$qc_pass_pre_doublet & (is_singlet | (config$keep_doubletfinder_skipped & is_skipped))

  doublet_summary <- obj@meta.data |>
    mutate(sample = as.character(.data[[config$sample_col]])) |>
    filter(qc_pass_pre_doublet) |>
    count(sample, doubletfinder_class, name = "n_cells") |>
    group_by(sample) |>
    mutate(percent = n_cells / sum(n_cells) * 100) |>
    ungroup()
  write.csv(doublet_summary, file.path(config$out_dir, "tables", "DoubletFinder_summary_by_sample.csv"), row.names = FALSE)
  plot_bar(
    doublet_summary,
    x = "sample",
    y = "n_cells",
    fill = "doubletfinder_class",
    out_file = file.path(config$out_dir, "figures", "DoubletFinder_summary_by_sample.pdf")
  )
} else {
  obj$doubletfinder_class <- NA_character_
  obj$qc_pass_final <- obj$qc_pass_pre_doublet
}

final_counts <- obj@meta.data |>
  mutate(
    sample = as.character(.data[[config$sample_col]]),
    final_status = if_else(qc_pass_final, "Final_pass", "Filtered")
  ) |>
  count(sample, final_status, name = "n_cells") |>
  group_by(sample) |>
  mutate(percent = n_cells / sum(n_cells) * 100) |>
  ungroup()
write.csv(final_counts, file.path(config$out_dir, "tables", "final_filtering_summary_by_sample.csv"), row.names = FALSE)
plot_bar(
  final_counts,
  x = "sample",
  y = "n_cells",
  fill = "final_status",
  out_file = file.path(config$out_dir, "figures", "final_filtering_summary_by_sample.pdf")
)

if (config$annotation_col %in% colnames(obj@meta.data)) {
  comp <- obj@meta.data |>
    filter(qc_pass_final) |>
    mutate(
      sample = as.character(.data[[config$sample_col]]),
      annotation = as.character(.data[[config$annotation_col]])
    ) |>
    count(sample, annotation, name = "n_cells") |>
    group_by(sample) |>
    mutate(percent = n_cells / sum(n_cells) * 100) |>
    ungroup()
  write.csv(comp, file.path(config$out_dir, "tables", "cell_type_composition_after_QC.csv"), row.names = FALSE)

  p_comp <- ggplot(comp, aes(x = sample, y = percent, fill = annotation)) +
    geom_col(width = 0.78) +
    labs(x = NULL, y = "Cell proportion after QC (%)", fill = NULL) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  ggsave(file.path(config$out_dir, "figures", "cell_type_composition_after_QC_percent.pdf"), p_comp, width = 10, height = 5.8, units = "in", dpi = 300)
}

filtered_obj <- subset(obj, cells = colnames(obj)[obj$qc_pass_final])
qs2::qs_save(obj, file.path(config$out_dir, "outputs", "GBM.RNA.qc_doubletfinder_annotated.qs2"))
qs2::qs_save(filtered_obj, file.path(config$out_dir, "outputs", "GBM.RNA.qc_doubletfinder.filtered.qs2"))

msg("Done.")
