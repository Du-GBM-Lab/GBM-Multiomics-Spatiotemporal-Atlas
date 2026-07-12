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
  gene_order_file = normalizePath(file.path("reference", "hg38_gencode_v27.txt"), winslash = "\\", mustWork = FALSE),
  reference_groups = c(
    "T cells", "NK cells", "B cells",
    "Macrophages", "Microglial", "Monocytes", "cDCs", "pDCs",
    "Endothelial", "Mural cells"
  ),
  min_reference_cells = 50,
  cutoff = 0.1,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = FALSE,
  num_threads = 4,
  keep_intermediate = FALSE
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

if (!requireNamespace("infercnv", quietly = TRUE)) {
  stop(
    "Required package not installed: infercnv\n",
    "Install infercnv before running this script. This script is written but not executed automatically.",
    call. = FALSE
  )
}

if (!file.exists(config$gene_order_file)) {
  stop(
    "Please set config$gene_order_file to a real human gene order file before running inferCNV.\n",
    "The file should contain gene, chromosome, start, and end columns compatible with infercnv.",
    call. = FALSE
  )
}

if (!file.exists(config$input_qs2)) {
  stop("QC-filtered qs2 object not found. Run 02_scRNA_QC first: ", config$input_qs2, call. = FALSE)
}
input_path <- config$input_qs2

msg("Loading object:", input_path)
obj <- read_object(input_path)
obj <- join_rna_layers_if_needed(obj, config$assay)

md <- obj@meta.data
stopifnot(config$sample_col %in% colnames(md))
stopifnot(config$annotation_col %in% colnames(md))

if ("qc_pass_final" %in% colnames(md)) {
  obj <- subset(obj, cells = colnames(obj)[obj$qc_pass_final])
  md <- obj@meta.data
}

samples <- sort(unique(as.character(md[[config$sample_col]])))
manifest <- list()

for (sample_id in samples) {
  msg("Preparing sample-wise inferCNV:", sample_id)
  sample_cells <- colnames(obj)[as.character(obj@meta.data[[config$sample_col]]) == sample_id]
  sample_obj <- subset(obj, cells = sample_cells)
  sample_md <- sample_obj@meta.data

  sample_md$infercnv_group <- as.character(sample_md[[config$annotation_col]])
  ref_groups_present <- intersect(config$reference_groups, unique(sample_md$infercnv_group))
  n_ref <- sum(sample_md$infercnv_group %in% ref_groups_present)

  sample_out <- file.path(config$out_dir, "outputs", sample_id)
  dir.create(sample_out, showWarnings = FALSE, recursive = TRUE)

  if (n_ref < config$min_reference_cells) {
    manifest[[sample_id]] <- data.frame(
      sample = sample_id,
      n_cells = ncol(sample_obj),
      n_reference_cells = n_ref,
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
      keep_files <- c(
        "run.final.infercnv_obj",
        basename(annotation_file),
        "infercnv.png",
        "infercnv_subclusters.png"
      )
      remove_patterns <- c(
        "^[0-9]+_.*\\.infercnv_obj$",
        "preliminary",
        "_HMM_",
        "\\.dat$",
        "\\.txt$",
        "_chr_",
        "infercnv\\.[0-9]+_",
        "BayesNetOutput"
      )
      to_remove <- list.files(
        sample_out,
        pattern = paste(remove_patterns, collapse = "|"),
        full.names = TRUE,
        recursive = FALSE
      )
      to_remove <- to_remove[!basename(to_remove) %in% keep_files]
      if (length(to_remove) > 0) {
        unlink(to_remove, recursive = TRUE, force = TRUE)
      }
    }
    "completed"
  }, error = function(e) {
    msg("inferCNV failed for", sample_id, ":", conditionMessage(e))
    paste0("failed: ", conditionMessage(e))
  })

  if (!identical(infercnv_status, "completed")) {
    manifest[[sample_id]] <- data.frame(
      sample = sample_id,
      n_cells = ncol(sample_obj),
      n_reference_cells = n_ref,
      reference_groups = paste(ref_groups_present, collapse = ";"),
      status = infercnv_status,
      out_dir = sample_out,
      stringsAsFactors = FALSE
    )
    next
  }

  manifest[[sample_id]] <- data.frame(
    sample = sample_id,
    n_cells = ncol(sample_obj),
    n_reference_cells = n_ref,
    reference_groups = paste(ref_groups_present, collapse = ";"),
    status = "completed",
    out_dir = sample_out,
    stringsAsFactors = FALSE
  )
}

manifest_df <- bind_rows(manifest)
write.csv(manifest_df, file.path(config$out_dir, "tables", "samplewise_infercnv_manifest.csv"), row.names = FALSE)

msg("Done.")
