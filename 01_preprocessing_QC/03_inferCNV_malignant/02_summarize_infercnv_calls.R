suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(qs2)
  library(dplyr)
  library(ggplot2)
})

config <- list(
  input_qs2 = normalizePath(file.path("..", "02_scRNA_QC", "outputs", "GBM.RNA.qc_doubletfinder.filtered.qs2"), winslash = "\\", mustWork = FALSE),
  infercnv_dir = normalizePath(file.path(".", "outputs"), winslash = "\\", mustWork = TRUE),
  out_dir = normalizePath(".", winslash = "\\", mustWork = TRUE),
  sample_col = "Pt_number",
  annotation_col = "anno_ident",
  reference_groups = c(
    "T cells", "NK cells", "B cells",
    "Macrophages", "Microglial", "Monocytes", "cDCs", "pDCs",
    "Endothelial", "Mural cells"
  ),
  cnv_threshold_nmads = 3,
  sensitivity_nmads = c(2, 3, 5),
  min_putative_malignant_cells = 20,
  putative_top_fraction = 0.1,
  min_reference_cells = 50,
  marker_validation_genes = c(
    "EGFR", "SOX2", "OLIG2", "GFAP",
    "PDGFRA", "CD44", "CHI3L1", "VIM",
    "MKI67",
    "PTPRC", "CD3D", "NKG7", "CD68", "CX3CR1", "PECAM1",
    "PTEN", "CDKN2A",
    "HSPA1A", "HSPA1B", "JUN", "FOS"
  )
)

dir.create(file.path(config$out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$out_dir, "outputs"), showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

read_object <- function(path) {
  if (grepl("\\.qs2$", path, ignore.case = TRUE)) {
    return(qs2::qs_read(path))
  }
  readRDS(path)
}

join_rna_layers_if_needed <- function(obj, assay = "RNA") {
  obj <- UpdateSeuratObject(obj)
  DefaultAssay(obj) <- assay
  if (inherits(obj[[assay]], "Assay5")) {
    obj[[assay]] <- JoinLayers(obj[[assay]])
  }
  obj
}

robust_upper <- function(x, nmads = 3) {
  med <- median(x, na.rm = TRUE)
  mad_value <- mad(x, constant = 1, na.rm = TRUE)
  if (is.na(mad_value) || mad_value == 0) {
    mad_value <- stats::IQR(x, na.rm = TRUE) / 1.349
  }
  if (is.na(mad_value) || mad_value == 0) {
    return(stats::quantile(x, 0.95, na.rm = TRUE))
  }
  med + nmads * mad_value
}

robust_center_scale <- function(x) {
  med <- median(x, na.rm = TRUE)
  mad_value <- mad(x, constant = 1, na.rm = TRUE)
  if (is.na(mad_value) || mad_value == 0) {
    mad_value <- stats::IQR(x, na.rm = TRUE) / 1.349
  }
  if (is.na(mad_value) || mad_value == 0) {
    mad_value <- stats::sd(x, na.rm = TRUE)
  }
  if (is.na(mad_value) || mad_value == 0) {
    mad_value <- NA_real_
  }
  list(center = med, scale = mad_value)
}

find_infercnv_rds <- function(sample_dir) {
  candidates <- c(
    list.files(sample_dir, pattern = "_infercnv_obj\\.qs2$", full.names = TRUE),
    list.files(sample_dir, pattern = "_infercnv_obj\\.rds$", full.names = TRUE),
    list.files(sample_dir, pattern = "^run\\.final\\.infercnv_obj$", full.names = TRUE)
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) {
    return(NA_character_)
  }
  candidates[1]
}

infercnv_cell_scores <- function(infercnv_obj) {
  expr <- infercnv_obj@expr.data
  # inferCNV expression is centered around the reference baseline after smoothing.
  # This is a relative CNV burden score, not an absolute copy-number estimate.
  score <- colMeans((expr - 1)^2)
  data.frame(
    cell = names(score),
    cnv_burden = as.numeric(score),
    n_genes_infercnv = nrow(expr),
    stringsAsFactors = FALSE
  )
}

safe_correlation_to_profile <- function(expr, profile) {
  profile <- as.numeric(profile)
  profile_sd <- stats::sd(profile, na.rm = TRUE)
  if (is.na(profile_sd) || profile_sd == 0) {
    return(rep(NA_real_, ncol(expr)))
  }
  apply(expr, 2, function(x) {
    suppressWarnings(stats::cor(as.numeric(x), profile, method = "pearson", use = "pairwise.complete.obs"))
  })
}

input_path <- config$input_qs2
if (!file.exists(input_path)) {
  stop("QC-filtered qs2 object not found. Run 02_scRNA_QC first: ", input_path, call. = FALSE)
}

msg("Loading object:", input_path)
obj <- read_object(input_path)
obj <- join_rna_layers_if_needed(obj, "RNA")
md <- obj@meta.data
stopifnot(config$sample_col %in% colnames(md))
stopifnot(config$annotation_col %in% colnames(md))

samples <- sort(unique(as.character(md[[config$sample_col]])))
all_scores <- list()
threshold_rows <- list()

for (sample_id in samples) {
  sample_dir <- file.path(config$infercnv_dir, sample_id)
  infercnv_rds <- find_infercnv_rds(sample_dir)
  if (is.na(infercnv_rds)) {
    warning("No inferCNV RDS found for sample: ", sample_id)
    next
  }

  msg("Summarizing inferCNV:", sample_id)
  infercnv_obj <- read_object(infercnv_rds)
  expr <- infercnv_obj@expr.data
  scores <- infercnv_cell_scores(infercnv_obj)
  scores$sample <- sample_id

  scores$annotation <- as.character(md[scores$cell, config$annotation_col])
  scores$is_reference_group <- scores$annotation %in% config$reference_groups

  ref_scores <- scores$cnv_burden[scores$is_reference_group]
  if (length(ref_scores) < config$min_reference_cells) {
    threshold <- NA_real_
    ref_center <- NA_real_
    ref_scale <- NA_real_
    cor_threshold <- NA_real_
    scores$cnv_burden_z <- NA_real_
    scores$cnv_correlation <- NA_real_
    scores$cnv_correlation_ref_z <- NA_real_
    scores$infercnv_burden_call <- "undetermined_low_reference"
    scores$infercnv_call <- "undetermined_low_reference"
    scores$is_malignant_for_downstream <- FALSE
  } else {
    ref_stats <- robust_center_scale(ref_scores)
    ref_center <- ref_stats$center
    ref_scale <- ref_stats$scale
    scores$cnv_burden_z <- (scores$cnv_burden - ref_center) / ref_scale
    threshold <- ref_center + config$cnv_threshold_nmads * ref_scale

    scores$infercnv_burden_call <- ifelse(
      scores$is_reference_group,
      "non_malignant_reference",
      ifelse(
        scores$cnv_burden_z > config$cnv_threshold_nmads,
        "malignant_like_CNV_high",
        "non_malignant_like_CNV_low"
      )
    )

    non_ref <- which(!scores$is_reference_group)
    putative <- non_ref[which(scores$cnv_burden_z[non_ref] > config$cnv_threshold_nmads)]
    if (length(putative) < config$min_putative_malignant_cells && length(non_ref) > 0) {
      n_top <- max(config$min_putative_malignant_cells, ceiling(length(non_ref) * config$putative_top_fraction))
      n_top <- min(n_top, length(non_ref))
      putative <- non_ref[order(scores$cnv_burden_z[non_ref], decreasing = TRUE)[seq_len(n_top)]]
    }

    if (length(putative) > 0) {
      mean_cnv_profile <- rowMeans(expr[, scores$cell[putative], drop = FALSE])
      cnv_cor <- safe_correlation_to_profile(expr[, scores$cell, drop = FALSE], mean_cnv_profile)
      scores$cnv_correlation <- as.numeric(cnv_cor)

      ref_cor <- scores$cnv_correlation[scores$is_reference_group]
      ref_cor_stats <- robust_center_scale(ref_cor)
      cor_threshold <- ref_cor_stats$center + config$cnv_threshold_nmads * ref_cor_stats$scale
      scores$cnv_correlation_ref_z <- (scores$cnv_correlation - ref_cor_stats$center) / ref_cor_stats$scale
      two_axis_malignant <- !scores$is_reference_group &
        scores$cnv_burden_z > config$cnv_threshold_nmads &
        scores$cnv_correlation_ref_z > config$cnv_threshold_nmads
    } else {
      cor_threshold <- NA_real_
      scores$cnv_correlation <- NA_real_
      scores$cnv_correlation_ref_z <- NA_real_
      two_axis_malignant <- rep(FALSE, nrow(scores))
    }

    scores$infercnv_call <- ifelse(
      scores$is_reference_group,
      "non_malignant_reference",
      ifelse(
        two_axis_malignant,
        "malignant_like_CNV_high_confidence",
        ifelse(
          scores$cnv_burden_z > config$cnv_threshold_nmads,
          "malignant_like_CNV_burden_only",
          "non_malignant_like_CNV_low"
        )
      )
    )
    scores$is_malignant_for_downstream <- two_axis_malignant
  }

  threshold_rows[[sample_id]] <- data.frame(
    sample = sample_id,
    n_cells = nrow(scores),
    n_reference_cells = length(ref_scores),
    cnv_burden_threshold = threshold,
    cnv_burden_ref_median = ref_center,
    cnv_burden_ref_mad = ref_scale,
    cnv_burden_z_threshold = config$cnv_threshold_nmads,
    cnv_correlation_ref_z_threshold = config$cnv_threshold_nmads,
    cnv_correlation_threshold = cor_threshold,
    stringsAsFactors = FALSE
  )
  all_scores[[sample_id]] <- scores
}

score_df <- bind_rows(all_scores)
threshold_df <- bind_rows(threshold_rows)
if (nrow(score_df) == 0) {
  stop("No inferCNV cell scores were found. Check sample-wise inferCNV outputs first.", call. = FALSE)
}

write.csv(score_df, file.path(config$out_dir, "tables", "infercnv_cell_cnv_burden_and_calls.csv"), row.names = FALSE)
write.csv(threshold_df, file.path(config$out_dir, "tables", "infercnv_samplewise_thresholds.csv"), row.names = FALSE)

sensitivity_df <- bind_rows(lapply(config$sensitivity_nmads, function(nmads) {
  score_df |>
    filter(!is_reference_group, !is.na(cnv_burden_z)) |>
    group_by(sample) |>
    summarise(
      threshold_nmads = nmads,
      n_non_reference = n(),
      n_malignant_like = sum(cnv_burden_z > nmads, na.rm = TRUE),
      malignant_like_percent = n_malignant_like / n_non_reference * 100,
      .groups = "drop"
    )
}))
write.csv(sensitivity_df, file.path(config$out_dir, "tables", "infercnv_threshold_sensitivity_2_3_5MAD.csv"), row.names = FALSE)

call_summary <- score_df |>
  mutate(
    infercnv_call = factor(
      infercnv_call,
      levels = c(
        "non_malignant_reference",
        "non_malignant_like_CNV_low",
        "malignant_like_CNV_burden_only",
        "malignant_like_CNV_high_confidence",
        "undetermined_low_reference"
      )
    )
  ) |>
  count(sample, infercnv_call, name = "n_cells") |>
  group_by(sample) |>
  mutate(percent = n_cells / sum(n_cells) * 100) |>
  ungroup()
write.csv(call_summary, file.path(config$out_dir, "tables", "infercnv_call_summary_by_sample.csv"), row.names = FALSE)

annotation_call_cross_tab <- score_df |>
  filter(!is.na(infercnv_call), !is.na(annotation)) |>
  count(annotation, infercnv_call, name = "n_cells") |>
  group_by(annotation) |>
  mutate(percent_within_annotation = n_cells / sum(n_cells) * 100) |>
  ungroup()
write.csv(annotation_call_cross_tab, file.path(config$out_dir, "tables", "infercnv_call_by_annotation_cross_tab.csv"), row.names = FALSE)

sample_density_lines <- threshold_df |>
  select(sample, cnv_burden_z_threshold)

p_call <- ggplot(call_summary, aes(x = sample, y = percent, fill = infercnv_call)) +
  geom_col(width = 0.78) +
  scale_fill_manual(
    values = c(
      non_malignant_reference = "#587A8C",
      non_malignant_like_CNV_low = "#8E9A5B",
      malignant_like_CNV_burden_only = "#D8A03D",
      malignant_like_CNV_high_confidence = "#B24745",
      undetermined_low_reference = "#A7A7A7"
    ),
    drop = FALSE
  ) +
  labs(x = NULL, y = "Cell proportion (%)", fill = NULL) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
ggsave(file.path(config$out_dir, "figures", "infercnv_malignant_call_composition_by_sample.pdf"), p_call, width = 10, height = 5.8, units = "in", dpi = 300)

p_density <- ggplot(score_df |> filter(!is.na(cnv_burden_z)), aes(x = cnv_burden_z, fill = is_reference_group)) +
  geom_density(alpha = 0.45, linewidth = 0.25) +
  geom_vline(data = sample_density_lines, aes(xintercept = cnv_burden_z_threshold), linetype = "dashed", linewidth = 0.25) +
  facet_wrap(~ sample, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c(`TRUE` = "#587A8C", `FALSE` = "#B24745"), labels = c(`TRUE` = "Reference", `FALSE` = "Non-reference")) +
  labs(x = "CNV burden z-score relative to reference", y = "Density", fill = NULL) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(config$out_dir, "figures", "infercnv_cnv_burden_z_density_by_sample.pdf"), p_density, width = 11, height = 7.5, units = "in", dpi = 300)

p_two_axis <- ggplot(score_df |> filter(!is.na(cnv_burden_z), !is.na(cnv_correlation_ref_z)), aes(x = cnv_burden_z, y = cnv_correlation_ref_z, color = infercnv_call)) +
  geom_point(size = 0.25, alpha = 0.45) +
  geom_vline(xintercept = config$cnv_threshold_nmads, linetype = "dashed", linewidth = 0.25) +
  geom_hline(yintercept = config$cnv_threshold_nmads, linetype = "dashed", linewidth = 0.25) +
  facet_wrap(~ sample, scales = "free", ncol = 6) +
  scale_color_manual(
    values = c(
      non_malignant_reference = "#587A8C",
      non_malignant_like_CNV_low = "#8E9A5B",
      malignant_like_CNV_burden_only = "#D8A03D",
      malignant_like_CNV_high_confidence = "#B24745",
      undetermined_low_reference = "#A7A7A7"
    ),
    drop = FALSE
  ) +
  labs(x = "CNV burden z-score", y = "CNV correlation z-score", color = NULL) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
ggsave(file.path(config$out_dir, "figures", "infercnv_two_axis_cnv_burden_correlation_by_sample.pdf"), p_two_axis, width = 11, height = 7.5, units = "in", dpi = 300)

p_cross <- ggplot(annotation_call_cross_tab, aes(x = infercnv_call, y = annotation, fill = percent_within_annotation)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "#F4F1EA", high = "#B24745", name = "% within\nannotation") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
ggsave(file.path(config$out_dir, "figures", "infercnv_call_by_annotation_cross_tab_heatmap.pdf"), p_cross, width = 8.2, height = 5.8, units = "in", dpi = 300)

obj$infercnv_cnv_burden <- NA_real_
obj$infercnv_cnv_burden_z <- NA_real_
obj$infercnv_cnv_correlation <- NA_real_
obj$infercnv_cnv_correlation_ref_z <- NA_real_
obj$infercnv_burden_call <- NA_character_
obj$infercnv_call <- NA_character_
obj$is_malignant_for_downstream <- FALSE
rownames(score_df) <- score_df$cell
shared_cells <- intersect(colnames(obj), score_df$cell)
obj@meta.data[shared_cells, "infercnv_cnv_burden"] <- score_df[shared_cells, "cnv_burden"]
obj@meta.data[shared_cells, "infercnv_cnv_burden_z"] <- score_df[shared_cells, "cnv_burden_z"]
obj@meta.data[shared_cells, "infercnv_cnv_correlation"] <- score_df[shared_cells, "cnv_correlation"]
obj@meta.data[shared_cells, "infercnv_cnv_correlation_ref_z"] <- score_df[shared_cells, "cnv_correlation_ref_z"]
obj@meta.data[shared_cells, "infercnv_burden_call"] <- score_df[shared_cells, "infercnv_burden_call"]
obj@meta.data[shared_cells, "infercnv_call"] <- score_df[shared_cells, "infercnv_call"]
obj@meta.data[shared_cells, "is_malignant_for_downstream"] <- score_df[shared_cells, "is_malignant_for_downstream"]
qs2::qs_save(obj, file.path(config$out_dir, "outputs", "GBM.RNA.qc_doubletfinder.infercnv_calls.qs2"))

genes_present <- intersect(config$marker_validation_genes, rownames(obj))
if (length(genes_present) > 0) {
  DefaultAssay(obj) <- "RNA"
  if (!"data" %in% Layers(obj[["RNA"]])) {
    obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
  }
  marker_plot <- DotPlot(
    obj,
    features = genes_present,
    group.by = "infercnv_call",
    assay = "RNA"
  ) +
    RotatedAxis() +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank())
  ggsave(file.path(config$out_dir, "figures", "marker_validation_by_infercnv_call.pdf"), marker_plot, width = 8, height = 4.8, units = "in", dpi = 300)

  obj$infercnv_call_annotation <- ifelse(
    is.na(obj$infercnv_call) | is.na(obj@meta.data[[config$annotation_col]]),
    NA_character_,
    paste(obj$infercnv_call, obj@meta.data[[config$annotation_col]], sep = " | ")
  )
  group_sizes <- sort(table(obj$infercnv_call_annotation), decreasing = TRUE)
  keep_groups <- names(group_sizes)[group_sizes >= 20]
  if (length(keep_groups) > 0) {
    marker_split_plot <- DotPlot(
      subset(obj, cells = rownames(obj@meta.data)[obj$infercnv_call_annotation %in% keep_groups]),
      features = genes_present,
      group.by = "infercnv_call_annotation",
      assay = "RNA"
    ) +
      RotatedAxis() +
      labs(x = NULL, y = NULL) +
      theme_bw(base_size = 7) +
      theme(panel.grid.minor = element_blank())
    ggsave(file.path(config$out_dir, "figures", "marker_validation_by_infercnv_call_and_annotation.pdf"), marker_split_plot, width = 10, height = 8.5, units = "in", dpi = 300)
  }
}

msg("Done.")
