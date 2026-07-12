suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(qs2)
  library(dplyr)
})

config <- list(
  input_qs2 = normalizePath(file.path("..", "02_scRNA_QC", "outputs", "GBM.RNA.qc_doubletfinder.filtered.qs2"), winslash = "\\", mustWork = FALSE),
  out_dir = normalizePath(".", winslash = "\\", mustWork = TRUE),
  sample_col = "Pt_number",
  annotation_col = "anno_ident",
  assay = "RNA",
  gene_order_file = normalizePath(file.path("..", "03_inferCNV_恶性识别", "reference", "hg38_gencode_v27.txt"), winslash = "\\", mustWork = FALSE),
  reference_tier1 = c("T cells", "NK cells", "B cells", "pDCs"),
  reference_tier2 = c("Microglial"),
  reference_tier3 = c("Macrophages", "Monocytes", "cDCs"),
  min_reference_cells = 50,
  cutoff = 0.1,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = FALSE,
  num_threads = 4,
  keep_intermediate = FALSE,
  rerun_completed = FALSE
)

dir.create(file.path(config$out_dir, "outputs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(config$out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

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

cleanup_infercnv_intermediates <- function(sample_out, annotation_file) {
  keep_regex <- paste(c(
    "^run\\.final\\.infercnv_obj$",
    paste0("^", gsub("([.])", "\\\\\\1", basename(annotation_file)), "$"),
    "^infercnv\\.observations\\.txt$",
    "^infercnv\\.references\\.txt$",
    "^infercnv\\.png$",
    "^infercnv_subclusters\\.png$"
  ), collapse = "|")
  all_files <- list.files(sample_out, full.names = TRUE, recursive = FALSE, all.files = FALSE)
  to_remove <- all_files[!grepl(keep_regex, basename(all_files))]
  if (length(to_remove) > 0) {
    unlink(to_remove, recursive = TRUE, force = TRUE)
  }
}

get_sample_reference <- function(sample_md, annotation_col, config) {
  groups <- unique(as.character(sample_md[[annotation_col]]))
  annotations <- as.character(sample_md[[annotation_col]])

  ref <- intersect(config$reference_tier1, groups)
  n <- sum(annotations %in% ref)
  if (n >= config$min_reference_cells) {
    return(list(ref = ref, tier = "lymphoid_pdc_only", n = n))
  }

  ref <- intersect(c(config$reference_tier1, config$reference_tier2), groups)
  n <- sum(annotations %in% ref)
  if (n >= config$min_reference_cells) {
    return(list(ref = ref, tier = "lymphoid_pdc_plus_microglia", n = n))
  }

  ref <- intersect(c(config$reference_tier1, config$reference_tier2, config$reference_tier3), groups)
  n <- sum(annotations %in% ref)
  if (n >= config$min_reference_cells) {
    return(list(ref = ref, tier = "all_immune_fallback", n = n))
  }

  list(ref = ref, tier = "insufficient", n = n)
}

if (!requireNamespace("infercnv", quietly = TRUE)) {
  stop("Required package not installed: infercnv", call. = FALSE)
}
if (!file.exists(config$input_qs2)) {
  stop("QC-filtered qs2 object not found. Run 02_scRNA_QC first: ", config$input_qs2, call. = FALSE)
}
if (!file.exists(config$gene_order_file)) {
  stop("Gene order file not found: ", config$gene_order_file, call. = FALSE)
}

msg("Loading QC-filtered object:", config$input_qs2)
obj <- read_object(config$input_qs2)
obj <- join_rna_layers_if_needed(obj, config$assay)

md <- obj@meta.data
stopifnot(config$sample_col %in% colnames(md))
stopifnot(config$annotation_col %in% colnames(md))
if ("qc_pass_final" %in% colnames(md)) {
  keep_cells <- colnames(obj)[which(obj$qc_pass_final)]
  obj <- subset(obj, cells = keep_cells)
  md <- obj@meta.data
}

samples <- sort(unique(as.character(md[[config$sample_col]])))
manifest <- list()

for (sample_id in samples) {
  msg("Preparing immune-reference inferCNV:", sample_id)
  sample_cells <- colnames(obj)[as.character(obj@meta.data[[config$sample_col]]) == sample_id]
  sample_obj <- subset(obj, cells = sample_cells)
  sample_md <- sample_obj@meta.data

  sample_md$infercnv_group <- as.character(sample_md[[config$annotation_col]])
  ref_info <- get_sample_reference(sample_md, config$annotation_col, config)
  ref_groups_present <- ref_info$ref
  n_ref <- ref_info$n

  sample_out <- file.path(config$out_dir, "outputs", sample_id)
  dir.create(sample_out, showWarnings = FALSE, recursive = TRUE)

  final_obj <- file.path(sample_out, "run.final.infercnv_obj")
  if (file.exists(final_obj) && !config$rerun_completed) {
    manifest[[sample_id]] <- data.frame(
      sample = sample_id,
      n_cells = ncol(sample_obj),
      n_reference_cells = n_ref,
      reference_tier = ref_info$tier,
      reference_groups = paste(ref_groups_present, collapse = ";"),
      status = "existing_completed",
      out_dir = sample_out,
      stringsAsFactors = FALSE
    )
    next
  }

  if (identical(ref_info$tier, "insufficient")) {
    manifest[[sample_id]] <- data.frame(
      sample = sample_id,
      n_cells = ncol(sample_obj),
      n_reference_cells = n_ref,
      reference_tier = ref_info$tier,
      reference_groups = paste(ref_groups_present, collapse = ";"),
      status = "skipped_low_reference_cells",
      out_dir = sample_out,
      stringsAsFactors = FALSE
    )
    next
  }

  counts <- GetAssayData(sample_obj, assay = config$assay, layer = "counts")
  stopifnot(identical(rownames(sample_md), colnames(counts)))

  annotation_file <- file.path(sample_out, paste0(sample_id, "_cell_annotations.txt"))
  annotation_df <- data.frame(
    cell = colnames(counts),
    group = sample_md[colnames(counts), "infercnv_group"],
    stringsAsFactors = FALSE
  )
  stopifnot(identical(annotation_df$cell, colnames(counts)))
  write.table(annotation_df, annotation_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  infercnv_status <- tryCatch({
    infercnv_obj <- infercnv::CreateInfercnvObject(
      raw_counts_matrix = counts,
      annotations_file = annotation_file,
      delim = "\t",
      gene_order_file = config$gene_order_file,
      ref_group_names = ref_groups_present
    )

    infercnv::run(
      infercnv_obj,
      cutoff = config$cutoff,
      out_dir = sample_out,
      cluster_by_groups = config$cluster_by_groups,
      denoise = config$denoise,
      HMM = config$HMM,
      num_threads = config$num_threads
    )

    if (!config$keep_intermediate) {
      cleanup_infercnv_intermediates(sample_out, annotation_file)
    }
    "completed"
  }, error = function(e) {
    msg("inferCNV failed for", sample_id, ":", conditionMessage(e))
    paste0("failed: ", conditionMessage(e))
  })

  manifest[[sample_id]] <- data.frame(
    sample = sample_id,
    n_cells = ncol(sample_obj),
    n_reference_cells = n_ref,
    reference_tier = ref_info$tier,
    reference_groups = paste(ref_groups_present, collapse = ";"),
    status = infercnv_status,
    out_dir = sample_out,
    stringsAsFactors = FALSE
  )
}

manifest_df <- bind_rows(manifest)
write.csv(manifest_df, file.path(config$out_dir, "tables", "samplewise_infercnv_immune_reference_manifest.csv"), row.names = FALSE)

msg("Done.")
