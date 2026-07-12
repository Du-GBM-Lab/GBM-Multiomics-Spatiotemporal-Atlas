suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

config <- list(
  input_qs2 = normalizePath(file.path("..", "02_scRNA_QC", "outputs", "GBM.RNA.qc_doubletfinder.filtered.qs2"), winslash = "\\", mustWork = FALSE),
  infercnv_dir = normalizePath(file.path(".", "outputs"), winslash = "\\", mustWork = TRUE),
  manifest_csv = normalizePath(file.path("tables", "samplewise_infercnv_immune_reference_manifest.csv"), winslash = "\\", mustWork = FALSE),
  broad_score_csv = normalizePath(file.path("..", "03_inferCNV_恶性识别", "tables", "infercnv_cell_cnv_burden_and_calls.csv"), winslash = "\\", mustWork = FALSE),
  out_dir = normalizePath(".", winslash = "\\", mustWork = TRUE),
  sample_col = "Pt_number",
  annotation_col = "anno_ident",
  cnv_threshold_nmads = 3,
  sensitivity_nmads = c(2, 3, 5),
  min_reference_cells = 50,
  min_putative_malignant_cells = 20,
  putative_top_fraction = 0.1
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
  if (length(candidates) == 0) NA_character_ else candidates[1]
}

infercnv_cell_scores <- function(infercnv_obj) {
  expr <- infercnv_obj@expr.data
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
  as.numeric(suppressWarnings(stats::cor(profile, expr, method = "pearson", use = "pairwise.complete.obs")))
}

classify_malignant_status <- function(infercnv_call) {
  dplyr::case_when(
    infercnv_call == "malignant_like_CNV_high_confidence" ~ "malignant_high_confidence",
    infercnv_call == "malignant_like_CNV_burden_only" ~ "ambiguous_burden_only",
    infercnv_call %in% c("non_malignant_reference", "non_malignant_like_CNV_low") ~ "non_malignant",
    infercnv_call == "undetermined_low_reference" ~ "undetermined",
    TRUE ~ "not_evaluated"
  )
}

call_levels <- c(
  "non_malignant_reference",
  "non_malignant_like_CNV_low",
  "malignant_like_CNV_burden_only",
  "malignant_like_CNV_high_confidence",
  "undetermined_low_reference"
)

call_palette <- c(
  non_malignant_reference = "#587A8C",
  non_malignant_like_CNV_low = "#8E9A5B",
  malignant_like_CNV_burden_only = "#D8A03D",
  malignant_like_CNV_high_confidence = "#B24745",
  undetermined_low_reference = "#A7A7A7"
)

theme_cnv_pub <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold")
    )
}

if (!file.exists(config$input_qs2)) {
  stop("QC-filtered qs2 object not found. Run 02_scRNA_QC first: ", config$input_qs2, call. = FALSE)
}
if (!file.exists(config$manifest_csv)) {
  stop("Sample-wise inferCNV manifest not found. Run 01_samplewise_infercnv_immune_reference.R first: ", config$manifest_csv, call. = FALSE)
}

msg("Loading object:", config$input_qs2)
obj <- read_object(config$input_qs2)
obj <- join_rna_layers_if_needed(obj, "RNA")
md <- obj@meta.data
stopifnot(config$sample_col %in% colnames(md))
stopifnot(config$annotation_col %in% colnames(md))

manifest_df <- read.csv(config$manifest_csv, stringsAsFactors = FALSE)
required_manifest_cols <- c("sample", "n_reference_cells", "reference_tier", "reference_groups", "status", "out_dir")
if (!all(required_manifest_cols %in% colnames(manifest_df))) {
  stop("Manifest missing required columns: ", paste(setdiff(required_manifest_cols, colnames(manifest_df)), collapse = ", "), call. = FALSE)
}
reference_lookup <- manifest_df |>
  mutate(
    reference_groups_list = strsplit(reference_groups, ";", fixed = TRUE)
  )

samples <- sort(unique(as.character(md[[config$sample_col]])))
all_scores <- list()
threshold_rows <- list()

for (sample_id in samples) {
  manifest_row <- reference_lookup[reference_lookup$sample == sample_id, , drop = FALSE]
  if (nrow(manifest_row) == 0) {
    warning("No manifest row found for sample: ", sample_id)
    next
  }
  ref_groups_sample <- unlist(manifest_row$reference_groups_list[[1]])
  ref_groups_sample <- ref_groups_sample[nzchar(ref_groups_sample)]
  reference_tier <- manifest_row$reference_tier[[1]]
  manifest_status <- manifest_row$status[[1]]

  if (identical(reference_tier, "insufficient") || grepl("^skipped", manifest_status)) {
    warning("Skipping sample with insufficient reference cells: ", sample_id)
    next
  }

  sample_dir <- file.path(config$infercnv_dir, sample_id)
  infercnv_rds <- find_infercnv_rds(sample_dir)
  if (is.na(infercnv_rds)) {
    warning("No inferCNV RDS found for sample: ", sample_id)
    next
  }

  msg("Summarizing immune-reference inferCNV:", sample_id)
  infercnv_obj <- read_object(infercnv_rds)
  expr <- infercnv_obj@expr.data
  scores <- infercnv_cell_scores(infercnv_obj)
  scores$sample <- sample_id
  scores$reference_tier <- reference_tier
  scores$reference_groups <- paste(ref_groups_sample, collapse = ";")
  scores$annotation <- as.character(md[scores$cell, config$annotation_col])
  scores$is_reference_group <- scores$annotation %in% ref_groups_sample

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
    scores$malignant_call_status <- "undetermined"
    scores$is_evaluable_for_malignancy <- FALSE
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
      ifelse(scores$cnv_burden_z > config$cnv_threshold_nmads, "malignant_like_CNV_high", "non_malignant_like_CNV_low")
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
      scores$cnv_correlation <- as.numeric(safe_correlation_to_profile(expr[, scores$cell, drop = FALSE], mean_cnv_profile))

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
        ifelse(scores$cnv_burden_z > config$cnv_threshold_nmads, "malignant_like_CNV_burden_only", "non_malignant_like_CNV_low")
      )
    )
    scores$malignant_call_status <- classify_malignant_status(scores$infercnv_call)
    scores$is_evaluable_for_malignancy <- scores$malignant_call_status != "undetermined"
    scores$is_malignant_for_downstream <- two_axis_malignant
  }

  threshold_rows[[sample_id]] <- data.frame(
    sample = sample_id,
    n_cells = nrow(scores),
    n_reference_cells = length(ref_scores),
    reference_tier = reference_tier,
    reference_groups = paste(ref_groups_sample, collapse = ";"),
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
  stop("No inferCNV cell scores were found. Run 01_samplewise_infercnv_immune_reference.R first.", call. = FALSE)
}

score_df$infercnv_call <- factor(score_df$infercnv_call, levels = call_levels)

write.csv(score_df, file.path(config$out_dir, "tables", "infercnv_immune_reference_cell_cnv_burden_and_calls.csv"), row.names = FALSE)
write.csv(threshold_df, file.path(config$out_dir, "tables", "infercnv_immune_reference_samplewise_thresholds.csv"), row.names = FALSE)

sensitivity_df <- bind_rows(lapply(config$sensitivity_nmads, function(nmads) {
  score_df |>
    filter(!is_reference_group, !is.na(cnv_burden_z)) |>
    group_by(sample) |>
    summarise(
      threshold_nmads = nmads,
      n_non_reference = n(),
      n_burden_only = sum(cnv_burden_z > nmads, na.rm = TRUE),
      n_two_axis = sum(cnv_burden_z > nmads & cnv_correlation_ref_z > nmads, na.rm = TRUE),
      burden_only_percent = n_burden_only / n_non_reference * 100,
      two_axis_percent = n_two_axis / n_non_reference * 100,
      .groups = "drop"
    )
}))
write.csv(sensitivity_df, file.path(config$out_dir, "tables", "infercnv_immune_reference_threshold_sensitivity_2_3_5MAD.csv"), row.names = FALSE)

call_summary <- score_df |>
  count(sample, infercnv_call, name = "n_cells") |>
  group_by(sample) |>
  mutate(percent = n_cells / sum(n_cells) * 100) |>
  ungroup()
write.csv(call_summary, file.path(config$out_dir, "tables", "infercnv_immune_reference_call_summary_by_sample.csv"), row.names = FALSE)

annotation_call_cross_tab <- score_df |>
  filter(!is.na(infercnv_call), !is.na(annotation)) |>
  count(annotation, infercnv_call, name = "n_cells") |>
  group_by(annotation) |>
  mutate(percent_within_annotation = n_cells / sum(n_cells) * 100) |>
  ungroup()
write.csv(annotation_call_cross_tab, file.path(config$out_dir, "tables", "infercnv_immune_reference_call_by_annotation_cross_tab.csv"), row.names = FALSE)

if (file.exists(config$broad_score_csv)) {
  broad_raw <- read.csv(config$broad_score_csv, stringsAsFactors = FALSE)
  required_cols <- c("cell", "infercnv_call")
  if (!all(required_cols %in% colnames(broad_raw))) {
    warning("Broad reference CSV missing required columns; skipping concordance: ", config$broad_score_csv)
  } else {
    broad_df <- broad_raw |>
      select(cell, broad_infercnv_call = infercnv_call) |>
      mutate(broad_high_confidence = broad_infercnv_call == "malignant_like_CNV_high_confidence")

    concordance_df <- score_df |>
      select(cell, immune_infercnv_call = infercnv_call) |>
      mutate(immune_high_confidence = immune_infercnv_call == "malignant_like_CNV_high_confidence") |>
      left_join(broad_df, by = "cell") |>
      mutate(
        concordance = case_when(
          is.na(broad_infercnv_call) ~ "missing_broad_call",
          immune_high_confidence & broad_high_confidence ~ "both_high_confidence",
          immune_high_confidence & !broad_high_confidence ~ "immune_only_high_confidence",
          !immune_high_confidence & broad_high_confidence ~ "broad_only_high_confidence",
          TRUE ~ "neither_high_confidence"
        )
      )

    concordance_summary <- concordance_df |>
      count(concordance, name = "n_cells") |>
      mutate(percent = n_cells / sum(n_cells) * 100)

    write.csv(concordance_df, file.path(config$out_dir, "tables", "infercnv_broad_vs_immune_reference_cell_concordance.csv"), row.names = FALSE)
    write.csv(concordance_summary, file.path(config$out_dir, "tables", "infercnv_broad_vs_immune_reference_concordance_summary.csv"), row.names = FALSE)
  }
}

p_call <- ggplot(call_summary, aes(x = sample, y = percent, fill = infercnv_call)) +
  geom_col(width = 0.78) +
  scale_fill_manual(values = call_palette, drop = FALSE) +
  labs(x = NULL, y = "Cell proportion (%)", fill = NULL) +
  theme_cnv_pub(base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "right")
ggsave(file.path(config$out_dir, "figures", "infercnv_immune_reference_call_composition_by_sample.pdf"), p_call, width = 10, height = 5.8, units = "in", dpi = 300, device = grDevices::cairo_pdf)

p_two_axis <- ggplot(score_df |> filter(!is.na(cnv_burden_z), !is.na(cnv_correlation_ref_z)), aes(x = cnv_burden_z, y = cnv_correlation_ref_z, color = infercnv_call)) +
  geom_point(size = 0.25, alpha = 0.45) +
  geom_vline(xintercept = config$cnv_threshold_nmads, linetype = "dashed", linewidth = 0.25) +
  geom_hline(yintercept = config$cnv_threshold_nmads, linetype = "dashed", linewidth = 0.25) +
  facet_wrap(~ sample, scales = "free", ncol = 6) +
  scale_color_manual(values = call_palette, drop = FALSE) +
  labs(x = "CNV burden z-score", y = "CNV correlation z-score", color = NULL) +
  theme_cnv_pub(base_size = 8) +
  theme(legend.position = "bottom")
ggsave(file.path(config$out_dir, "figures", "infercnv_immune_reference_two_axis_by_sample.pdf"), p_two_axis, width = 11, height = 7.5, units = "in", dpi = 300, device = grDevices::cairo_pdf)

p_cross <- ggplot(annotation_call_cross_tab, aes(x = infercnv_call, y = annotation, fill = percent_within_annotation)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(low = "#F4F1EA", high = "#B24745", name = "% within\nannotation") +
  labs(x = NULL, y = NULL) +
  theme_cnv_pub(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(config$out_dir, "figures", "infercnv_immune_reference_call_by_annotation_cross_tab_heatmap.pdf"), p_cross, width = 8.2, height = 5.8, units = "in", dpi = 300, device = grDevices::cairo_pdf)

obj$infercnv_immune_cnv_burden <- NA_real_
obj$infercnv_immune_cnv_burden_z <- NA_real_
obj$infercnv_immune_cnv_correlation <- NA_real_
obj$infercnv_immune_cnv_correlation_ref_z <- NA_real_
obj$infercnv_immune_burden_call <- NA_character_
obj$infercnv_immune_call <- NA_character_
obj$infercnv_immune_reference_tier <- NA_character_
obj$infercnv_immune_reference_groups <- NA_character_
obj$malignant_call_status_immune <- "not_evaluated"
obj$is_evaluable_for_malignancy_immune <- FALSE
obj$is_malignant_for_downstream_immune <- NA

rownames(score_df) <- score_df$cell
shared_cells <- intersect(colnames(obj), score_df$cell)
obj@meta.data[shared_cells, "infercnv_immune_cnv_burden"] <- score_df[shared_cells, "cnv_burden"]
obj@meta.data[shared_cells, "infercnv_immune_cnv_burden_z"] <- score_df[shared_cells, "cnv_burden_z"]
obj@meta.data[shared_cells, "infercnv_immune_cnv_correlation"] <- score_df[shared_cells, "cnv_correlation"]
obj@meta.data[shared_cells, "infercnv_immune_cnv_correlation_ref_z"] <- score_df[shared_cells, "cnv_correlation_ref_z"]
obj@meta.data[shared_cells, "infercnv_immune_burden_call"] <- as.character(score_df[shared_cells, "infercnv_burden_call"])
obj@meta.data[shared_cells, "infercnv_immune_call"] <- as.character(score_df[shared_cells, "infercnv_call"])
obj@meta.data[shared_cells, "infercnv_immune_reference_tier"] <- score_df[shared_cells, "reference_tier"]
obj@meta.data[shared_cells, "infercnv_immune_reference_groups"] <- score_df[shared_cells, "reference_groups"]
obj@meta.data[shared_cells, "malignant_call_status_immune"] <- score_df[shared_cells, "malignant_call_status"]
obj@meta.data[shared_cells, "is_evaluable_for_malignancy_immune"] <- score_df[shared_cells, "is_evaluable_for_malignancy"]
obj@meta.data[shared_cells, "is_malignant_for_downstream_immune"] <- score_df[shared_cells, "is_malignant_for_downstream"]

qs2::qs_save(obj, file.path(config$out_dir, "outputs", "GBM.RNA.qc_doubletfinder.infercnv_immune_reference_calls.qs2"))

msg("Done.")
