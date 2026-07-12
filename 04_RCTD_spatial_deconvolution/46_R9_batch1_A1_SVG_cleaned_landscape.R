#!/usr/bin/env Rscript
# =============================================================================
# R9 batch 1 | A1 SVG cleaned landscape
#
# Input is the already completed SPARK-X SVG table. This script does not rerun
# SPARK-X. It filters only pre-specified technical genes, audits expression
# dependence, checks locked MES-V and vascular marker-set SVG hits, and reports
# set-level high-expression region overlap without smoothing or correlation.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
})

base_dir <- getwd()
out_dir <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A1_SVG_cleaned")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

svg_long_path <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A1_SVG/A1_SVG_SPARKX_all_slices_long.csv")
svg_consensus_path <- file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A1_SVG/A1_SVG_SPARKX_consensus_summary.csv")
spatial_obj_path <- "<DATA_ROOT>/项目/分型/分型代码/0.对象/5.ST_merge.rds"
mesv_marker_path <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/05_恶性细胞分亚群与Neftel对照/tables/10d_S3_markers_top20.csv"
vascular_marker_path <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/05_恶性细胞分亚群与Neftel对照/tables/non_malignant_signature_genes_used.csv"

required_files <- c(svg_long_path, svg_consensus_path, spatial_obj_path, mesv_marker_path, vascular_marker_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) stop("Missing required input(s):\n", paste(missing_files, collapse = "\n"))

params <- data.table(
  parameter = c(
    "svg_hit_threshold",
    "technical_filter_mito",
    "technical_filter_ribo",
    "technical_filter_fixed_housekeeping",
    "technical_filter_hsp_family",
    "expression_drive_top_n",
    "high_expression_region_rule",
    "permutation_n",
    "random_seed"
  ),
  value = c(
    "adjustedPval < 0.05",
    "^MT-",
    "^RP[SL]",
    "MALAT1, NEAT1, FTL, FTH1, TMSB4X, TMSB10, B2M, ACTB, EEF1A1, EEF2, GAPDH, UBC",
    "^HSP",
    "top 200 and top 500 by adjustedPval within slice",
    "per marker gene, detected spots with raw count >= 90th percentile among detected spots; gene-set region is union across SVG-hit marker genes",
    "10000",
    "20260612"
  )
)
fwrite(params, file.path(out_dir, "A1_SVG_cleaned_parameters.csv"))

set.seed(20260612)
svg <- fread(svg_long_path)
consensus <- fread(svg_consensus_path)

needed_cols <- c("gene", "combinedPval", "adjustedPval", "slice", "n_spots", "n_genes_tested", "detected_fraction")
if (!all(needed_cols %in% names(svg))) {
  stop("A1 SVG long table lacks expected columns: ", paste(setdiff(needed_cols, names(svg)), collapse = ", "))
}

svg[, adjustedPval := as.numeric(adjustedPval)]
svg[, combinedPval := as.numeric(combinedPval)]
svg[, detected_fraction := as.numeric(detected_fraction)]
svg[, svg_rank := frank(adjustedPval, ties.method = "average", na.last = "keep"), by = slice]
svg[, neg_log10_fdr := -log10(pmax(adjustedPval, .Machine$double.xmin))]
svg[, is_svg_fdr_005 := adjustedPval < 0.05]

fixed_housekeeping <- c("MALAT1", "NEAT1", "FTL", "FTH1", "TMSB4X", "TMSB10",
                        "B2M", "ACTB", "EEF1A1", "EEF2", "GAPDH", "UBC")
svg[, filter_mito := grepl("^MT-", gene)]
svg[, filter_ribo := grepl("^RP[SL]", gene)]
svg[, filter_fixed_housekeeping := gene %in% fixed_housekeeping]
svg[, filter_hsp_family := grepl("^HSP", gene)]
svg[, technical_filtered := filter_mito | filter_ribo | filter_fixed_housekeeping | filter_hsp_family]
svg_clean <- svg[technical_filtered == FALSE]

filter_summary <- rbindlist(list(
  svg[, .(scope = "all_rows", n_rows = .N, n_genes = uniqueN(gene), n_svg_hits = sum(is_svg_fdr_005, na.rm = TRUE))],
  svg_clean[, .(scope = "cleaned_rows", n_rows = .N, n_genes = uniqueN(gene), n_svg_hits = sum(is_svg_fdr_005, na.rm = TRUE))]
), fill = TRUE)
fwrite(filter_summary, file.path(out_dir, "A1_filter_before_after_summary.csv"))

filter_by_category <- svg[, .(
  n_rows = .N,
  n_genes = uniqueN(gene),
  n_svg_hits = sum(is_svg_fdr_005, na.rm = TRUE),
  mito_rows = sum(filter_mito),
  ribo_rows = sum(filter_ribo),
  fixed_housekeeping_rows = sum(filter_fixed_housekeeping),
  hsp_family_rows = sum(filter_hsp_family),
  any_technical_rows = sum(technical_filtered)
), by = slice]
fwrite(filter_by_category, file.path(out_dir, "A1_filter_by_slice_category_counts.csv"))

clean_gene_summary <- svg_clean[, .(
  n_slices_tested = uniqueN(slice),
  n_slices_svg = sum(is_svg_fdr_005, na.rm = TRUE),
  best_fdr = min(adjustedPval, na.rm = TRUE),
  median_fdr = median(adjustedPval, na.rm = TRUE),
  mean_detected_fraction = mean(detected_fraction, na.rm = TRUE)
), by = gene][order(-n_slices_svg, best_fdr)]
fwrite(clean_gene_summary, file.path(out_dir, "A1_cleaned_SVG_gene_slice_hits.csv"))

raw_consensus_rank_col <- intersect(c("n_slices_svg", "n_slices_fdr_005", "n_slices_significant", "hit_slices"), names(consensus))[1]
if (is.na(raw_consensus_rank_col)) {
  consensus[, consensus_order := .I]
  top_raw <- consensus[order(consensus_order)][1:min(.N, 50)]
} else {
  top_raw <- consensus[order(-get(raw_consensus_rank_col))][1:min(.N, 50)]
}
fwrite(top_raw, file.path(out_dir, "A1_raw_consensus_top50_for_audit.csv"))
fwrite(clean_gene_summary[1:min(.N, 50)], file.path(out_dir, "A1_cleaned_consensus_top50.csv"))

message("Reading Spatial object: ", spatial_obj_path)
st <- readRDS(spatial_obj_path)
if (!"Spatial" %in% names(st@assays)) stop("Spatial assay not found in Seurat object.")

slice_levels <- sort(unique(svg$slice))
image_names <- names(st@images)
image_slice_map <- data.table(
  image = image_names,
  slice = vapply(image_names, function(img) {
    cells <- Cells(st@images[[img]])
    slice_hit <- unique(st$orig.ident[cells])
    if (length(slice_hit) != 1) NA_character_ else as.character(slice_hit)
  }, character(1))
)
image_slice_map <- image_slice_map[!is.na(slice)]
if (!all(slice_levels %in% image_slice_map$slice)) {
  stop("Could not map all SVG slices to Spatial images: ", paste(setdiff(slice_levels, image_slice_map$slice), collapse = ", "))
}
fwrite(image_slice_map, file.path(out_dir, "A1_SVG_cleaned_image_slice_map.csv"))

get_counts_for_slice <- function(slice_id) {
  img <- image_slice_map[slice == slice_id, image][1]
  cells <- Cells(st@images[[img]])
  layer_candidates <- c(paste0("counts.", slice_id), paste0("counts.", gsub("^#", "", slice_id)), "counts")
  layer_candidates <- layer_candidates[layer_candidates %in% Layers(st[["Spatial"]])]
  if (length(layer_candidates) == 0) {
    stop("No Spatial counts layer found for slice: ", slice_id)
  }
  mat <- LayerData(st[["Spatial"]], layer = layer_candidates[1])
  present_cells <- intersect(cells, colnames(mat))
  if (length(present_cells) == 0) stop("No image cells found in counts layer for slice: ", slice_id)
  mat[, present_cells, drop = FALSE]
}

expr_stats_list <- vector("list", length(slice_levels))
names(expr_stats_list) <- slice_levels
overlap_region_rows <- list()

message("Computing expression audit and set-level high-expression regions.")
for (slice_id in slice_levels) {
  mat <- get_counts_for_slice(slice_id)
  genes_needed <- intersect(svg[slice == slice_id, gene], rownames(mat))
  mat_sub <- mat[genes_needed, , drop = FALSE]
  gene_mean <- Matrix::rowMeans(mat_sub)
  gene_detected <- Matrix::rowSums(mat_sub > 0) / ncol(mat_sub)
  expr_stats_list[[slice_id]] <- data.table(
    slice = slice_id,
    gene = genes_needed,
    gene_mean_raw_count = as.numeric(gene_mean),
    detected_fraction_from_counts = as.numeric(gene_detected),
    n_cells_in_counts = ncol(mat_sub)
  )
}
expr_stats <- rbindlist(expr_stats_list, use.names = TRUE)
fwrite(expr_stats, file.path(out_dir, "A1_gene_expression_audit_by_slice.csv"))

svg <- merge(svg, expr_stats, by = c("slice", "gene"), all.x = TRUE, sort = FALSE)
svg_clean <- merge(svg_clean, expr_stats, by = c("slice", "gene"), all.x = TRUE, sort = FALSE)
fwrite(svg_clean, file.path(out_dir, "A1_SVG_SPARKX_all_slices_cleaned_long.csv"))

expr_drive_cor <- svg[, .(
  n_genes = sum(!is.na(svg_rank) & !is.na(gene_mean_raw_count)),
  spearman_rank_vs_mean = suppressWarnings(cor(svg_rank, gene_mean_raw_count, method = "spearman", use = "complete.obs")),
  spearman_rank_vs_detection = suppressWarnings(cor(svg_rank, detected_fraction_from_counts, method = "spearman", use = "complete.obs")),
  spearman_signal_vs_mean = suppressWarnings(cor(neg_log10_fdr, gene_mean_raw_count, method = "spearman", use = "complete.obs")),
  spearman_signal_vs_detection = suppressWarnings(cor(neg_log10_fdr, detected_fraction_from_counts, method = "spearman", use = "complete.obs"))
), by = slice]

expr_drive_top <- rbindlist(lapply(c(200L, 500L), function(n_top) {
  rbindlist(lapply(slice_levels, function(slice_id) {
    slice_all <- svg[slice == slice_id]
    mean_cut <- quantile(slice_all$gene_mean_raw_count, 0.9, na.rm = TRUE, names = FALSE)
    det_cut <- quantile(slice_all$detected_fraction_from_counts, 0.9, na.rm = TRUE, names = FALSE)
    slice_top <- slice_all[order(adjustedPval)][1:min(.N, n_top)]
    slice_top[, .(
      slice = slice_id,
      n_top = n_top,
      n_genes = .N,
      n_top_expression_decile = sum(gene_mean_raw_count >= mean_cut, na.rm = TRUE),
      n_top_detection_decile = sum(detected_fraction_from_counts >= det_cut, na.rm = TRUE),
      n_technical_filtered = sum(technical_filtered),
      n_svg_fdr_005 = sum(is_svg_fdr_005, na.rm = TRUE)
    )]
  }), use.names = TRUE)
}), use.names = TRUE)

fwrite(expr_drive_cor, file.path(out_dir, "A1_expression_drive_spearman_by_slice.csv"))
fwrite(expr_drive_top, file.path(out_dir, "A1_expression_drive_topSVG_audit.csv"))

mesv_markers <- unique(fread(mesv_marker_path)$gene)
vascular_sig <- fread(vascular_marker_path)
vascular_markers <- unique(vascular_sig[signature %in% c("Endothelial", "Pericyte", "vSMC") & used_in_object == TRUE, gene])

gene_set_provenance <- rbindlist(list(
  data.table(gene_set = "MES-V_program", gene = mesv_markers, source_path = mesv_marker_path,
             source_rule = "all genes in locked 10d_S3_markers_top20.csv"),
  data.table(gene_set = "vascular", gene = vascular_markers, source_path = vascular_marker_path,
             source_rule = "union of Endothelial, Pericyte, and vSMC signatures with used_in_object == TRUE")
), use.names = TRUE)
gene_set_provenance[, in_a1_tested_universe := gene %in% unique(svg$gene)]
gene_set_provenance[, in_cleaned_universe := gene %in% unique(svg_clean$gene)]
fwrite(gene_set_provenance, file.path(out_dir, "A1_target_gene_set_provenance.csv"))

universe_clean <- unique(svg_clean$gene)
svg_gene_any <- clean_gene_summary[n_slices_svg > 0, gene]

target_hit_long <- merge(
  gene_set_provenance[in_cleaned_universe == TRUE, .(gene_set, gene)],
  svg_clean[, .(slice, gene, adjustedPval, is_svg_fdr_005, detected_fraction, gene_mean_raw_count, detected_fraction_from_counts)],
  by = "gene",
  allow.cartesian = TRUE,
  all.x = TRUE,
  sort = FALSE
)
target_hit_long[, marker_svg_hit := is_svg_fdr_005 == TRUE]
fwrite(target_hit_long, file.path(out_dir, "A1_target_gene_set_SVG_hits_long.csv"))

target_gene_summary <- target_hit_long[, .(
  n_slices_tested = uniqueN(slice[!is.na(slice)]),
  n_slices_svg = sum(marker_svg_hit, na.rm = TRUE),
  best_fdr = min(adjustedPval, na.rm = TRUE),
  median_fdr = median(adjustedPval, na.rm = TRUE)
), by = .(gene_set, gene)][order(gene_set, -n_slices_svg, best_fdr)]
fwrite(target_gene_summary, file.path(out_dir, "A1_target_gene_set_gene_hit_summary.csv"))

target_slice_summary <- target_hit_long[, .(
  n_markers_in_cleaned_universe = uniqueN(gene),
  n_svg_markers = uniqueN(gene[marker_svg_hit == TRUE]),
  svg_marker_genes = paste(sort(unique(gene[marker_svg_hit == TRUE])), collapse = ";")
), by = .(gene_set, slice)][order(gene_set, slice)]
fwrite(target_slice_summary, file.path(out_dir, "A1_target_gene_set_per_slice_hits.csv"))

hypergeom_one <- function(markers, svg_hits, universe) {
  markers <- intersect(unique(markers), universe)
  svg_hits <- intersect(unique(svg_hits), universe)
  k <- length(intersect(markers, svg_hits))
  m <- length(svg_hits)
  n <- length(universe) - m
  q <- k - 1
  data.table(
    marker_n_in_universe = length(markers),
    svg_hit_genes_in_universe = m,
    marker_svg_overlap_n = k,
    hypergeom_p_enriched = phyper(q = q, m = m, n = n, k = length(markers), lower.tail = FALSE)
  )
}

perm_one <- function(markers, hit_table, universe, n_perm = 10000L) {
  markers <- intersect(unique(markers), universe)
  hit_counts <- hit_table[match(universe, gene), n_slices_svg]
  hit_counts[is.na(hit_counts)] <- 0L
  names(hit_counts) <- universe
  observed_overlap <- sum(markers %in% names(hit_counts)[hit_counts > 0])
  observed_sum_slices <- sum(hit_counts[markers])
  if (length(markers) == 0) {
    return(data.table(perm_p_overlap = NA_real_, perm_p_sum_slices = NA_real_,
                      observed_overlap = 0L, observed_sum_slices = 0L,
                      perm_mean_overlap = NA_real_, perm_mean_sum_slices = NA_real_))
  }
  perm_overlap <- integer(n_perm)
  perm_sum <- integer(n_perm)
  for (i in seq_len(n_perm)) {
    sampled <- sample(universe, length(markers), replace = FALSE)
    perm_overlap[i] <- sum(hit_counts[sampled] > 0)
    perm_sum[i] <- sum(hit_counts[sampled])
  }
  data.table(
    perm_p_overlap = (sum(perm_overlap >= observed_overlap) + 1) / (n_perm + 1),
    perm_p_sum_slices = (sum(perm_sum >= observed_sum_slices) + 1) / (n_perm + 1),
    observed_overlap = observed_overlap,
    observed_sum_slices = observed_sum_slices,
    perm_mean_overlap = mean(perm_overlap),
    perm_mean_sum_slices = mean(perm_sum)
  )
}

global_enrich <- rbindlist(lapply(unique(gene_set_provenance$gene_set), function(gs) {
  markers <- gene_set_provenance[gene_set == gs, gene]
  cbind(
    data.table(gene_set = gs, level = "global_any_slice"),
    hypergeom_one(markers, svg_gene_any, universe_clean),
    perm_one(markers, clean_gene_summary, universe_clean, n_perm = 10000L)
  )
}), fill = TRUE)
global_enrich[, hypergeom_fdr_bh := p.adjust(hypergeom_p_enriched, method = "BH")]
global_enrich[, perm_overlap_fdr_bh := p.adjust(perm_p_overlap, method = "BH")]
global_enrich[, perm_sum_slices_fdr_bh := p.adjust(perm_p_sum_slices, method = "BH")]
fwrite(global_enrich, file.path(out_dir, "A1_target_gene_set_global_enrichment.csv"))

slice_enrich <- rbindlist(lapply(slice_levels, function(slice_id) {
  universe_slice <- unique(svg_clean[slice == slice_id, gene])
  svg_hits_slice <- unique(svg_clean[slice == slice_id & is_svg_fdr_005 == TRUE, gene])
  rbindlist(lapply(unique(gene_set_provenance$gene_set), function(gs) {
    markers <- gene_set_provenance[gene_set == gs, gene]
    cbind(data.table(gene_set = gs, slice = slice_id, level = "per_slice"),
          hypergeom_one(markers, svg_hits_slice, universe_slice))
  }), fill = TRUE)
}), fill = TRUE)
slice_enrich[, hypergeom_fdr_bh_by_gene_set := p.adjust(hypergeom_p_enriched, method = "BH"), by = gene_set]
fwrite(slice_enrich, file.path(out_dir, "A1_target_gene_set_per_slice_enrichment.csv"))

high_region_for_genes <- function(mat, genes) {
  genes <- intersect(unique(genes), rownames(mat))
  if (length(genes) == 0) return(character(0))
  hit_cells <- character(0)
  for (g in genes) {
    vals <- as.numeric(mat[g, , drop = TRUE])
    positive <- vals > 0
    if (sum(positive) == 0) next
    cutoff <- as.numeric(stats::quantile(vals[positive], probs = 0.90, na.rm = TRUE, names = FALSE, type = 7))
    hit_cells <- union(hit_cells, colnames(mat)[positive & vals >= cutoff])
  }
  hit_cells
}

overlap_rows <- list()
for (slice_id in slice_levels) {
  mat <- get_counts_for_slice(slice_id)
  mesv_svg_genes <- target_hit_long[gene_set == "MES-V_program" & slice == slice_id & marker_svg_hit == TRUE, unique(gene)]
  vascular_svg_genes <- target_hit_long[gene_set == "vascular" & slice == slice_id & marker_svg_hit == TRUE, unique(gene)]
  mesv_region <- high_region_for_genes(mat, mesv_svg_genes)
  vascular_region <- high_region_for_genes(mat, vascular_svg_genes)
  inter_n <- length(intersect(mesv_region, vascular_region))
  union_n <- length(union(mesv_region, vascular_region))
  overlap_rows[[slice_id]] <- data.table(
    slice = slice_id,
    n_spots = ncol(mat),
    mesv_svg_gene_n = length(mesv_svg_genes),
    vascular_svg_gene_n = length(vascular_svg_genes),
    mesv_region_spot_n = length(mesv_region),
    vascular_region_spot_n = length(vascular_region),
    overlap_spot_n = inter_n,
    union_spot_n = union_n,
    jaccard = ifelse(union_n > 0, inter_n / union_n, NA_real_),
    overlap_fraction_of_mesv_region = ifelse(length(mesv_region) > 0, inter_n / length(mesv_region), NA_real_),
    overlap_fraction_of_vascular_region = ifelse(length(vascular_region) > 0, inter_n / length(vascular_region), NA_real_),
    mesv_svg_genes = paste(sort(mesv_svg_genes), collapse = ";"),
    vascular_svg_genes = paste(sort(vascular_svg_genes), collapse = ";")
  )
}
overlap_per_slice <- rbindlist(overlap_rows, fill = TRUE)
fwrite(overlap_per_slice, file.path(out_dir, "A1_target_SVG_set_level_region_overlap_per_slice.csv"))

overlap_summary <- overlap_per_slice[, .(
  n_slices = .N,
  n_slices_with_both_svg_sets = sum(mesv_svg_gene_n > 0 & vascular_svg_gene_n > 0),
  median_jaccard = median(jaccard, na.rm = TRUE),
  q25_jaccard = quantile(jaccard, 0.25, na.rm = TRUE),
  q75_jaccard = quantile(jaccard, 0.75, na.rm = TRUE),
  median_overlap_fraction_of_mesv_region = median(overlap_fraction_of_mesv_region, na.rm = TRUE),
  median_overlap_fraction_of_vascular_region = median(overlap_fraction_of_vascular_region, na.rm = TRUE)
)]
fwrite(overlap_summary, file.path(out_dir, "A1_target_SVG_set_level_region_overlap_summary.csv"))

stop_lines <- c(
  "# R9 Batch1 A1 SVG cleaned landscape STOP",
  "",
  paste0("- Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "- Scope: SPARK-X output cleanup only; SPARK-X was not rerun.",
  "- Technical filter: fixed rules only (`^MT-`, `^RP[SL]`, fixed housekeeping list, `^HSP`). Expression-driven audit was reported but not used for deletion.",
  "- Marker source: MES-V program from locked R6 `10d_S3_markers_top20.csv`; vascular from locked non-malignant Endothelial/Pericyte/vSMC signatures.",
  "- Spatial overlap rule: set-level high-expression region overlap counts only; no neighborhood smoothing and no continuous MES-V/vascular correlation.",
  "",
  "## Key output files",
  paste0("- ", normalizePath(file.path(out_dir, "A1_SVG_SPARKX_all_slices_cleaned_long.csv"), winslash = "/")),
  paste0("- ", normalizePath(file.path(out_dir, "A1_cleaned_SVG_gene_slice_hits.csv"), winslash = "/")),
  paste0("- ", normalizePath(file.path(out_dir, "A1_expression_drive_spearman_by_slice.csv"), winslash = "/")),
  paste0("- ", normalizePath(file.path(out_dir, "A1_target_gene_set_provenance.csv"), winslash = "/")),
  paste0("- ", normalizePath(file.path(out_dir, "A1_target_gene_set_gene_hit_summary.csv"), winslash = "/")),
  paste0("- ", normalizePath(file.path(out_dir, "A1_target_gene_set_global_enrichment.csv"), winslash = "/")),
  paste0("- ", normalizePath(file.path(out_dir, "A1_target_gene_set_per_slice_enrichment.csv"), winslash = "/")),
  paste0("- ", normalizePath(file.path(out_dir, "A1_target_SVG_set_level_region_overlap_per_slice.csv"), winslash = "/")),
  "",
  "## Self-placement placeholder",
  "- This is an A1 landscape/audit result. It can support supplementary landscape wording only if MES-V and vascular locked markers are both reproducibly present among cleaned SVGs and set-level overlap is cross-slice consistent.",
  "- It does not establish spatial relation-level evidence; relation-level claims remain delegated to locked v3-B crossK/nearest-distance/random-labeling null.",
  "- Forbidden wording remains: interaction, recruitment, drives, causality, single-cell colocalization, TAM recruitment."
)
writeLines(stop_lines, file.path(base_dir, "docs/R9_batch1_A1_SVG_cleaned_landscape_STOP.md"), useBytes = TRUE)

message("Done. Outputs written to: ", out_dir)
