options(stringsAsFactors = FALSE)
set.seed(42)

if (basename(getwd()) == "06_恶性细胞拟时序") {
  project_root <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(file.path(project_root, "06_恶性细胞拟时序"))

dir.create("tables", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/source_data", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

log_file <- file.path("logs", "08_tradeseq_driver_genes.log")
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  txt <- paste(...)
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", txt)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_msg("Phase 4 tradeSeq driver gene analysis started")
log_msg("set.seed = 42")

suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(Matrix)
  library(SingleCellExperiment)
  library(tradeSeq)
  library(BiocParallel)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

source("scripts/_naming.R")

obj_path <- file.path(project_root, "05_恶性细胞分亚群与Neftel对照", "outputs", "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2")
pseudotime_path <- file.path("tables", "slingshot_pseudotime_per_cell.csv")
weight_path <- file.path("tables", "slingshot_lineage_weight_per_cell.csv")
blacklist_path <- file.path("data", "gene_blacklist.txt")
endpoint_source_path <- file.path("figures/source_data", "06_panel_D_endpoint_boxplot.csv")

log_msg("Input object:", obj_path)
log_msg("Pseudotime table:", pseudotime_path)
log_msg("Lineage weight table:", weight_path)
log_msg("Gene blacklist:", blacklist_path)

read_blacklist_patterns <- function(path) {
  if (!file.exists(path)) return(character())
  x <- readLines(path, warn = FALSE, encoding = "UTF-8")
  x <- trimws(x)
  x <- x[nzchar(x) & !grepl("^#", x)]
  x
}

is_blacklisted <- function(genes, patterns) {
  if (length(patterns) == 0) return(rep(FALSE, length(genes)))
  vapply(genes, function(g) any(grepl(paste(patterns, collapse = "|"), g, perl = TRUE)), logical(1))
}

get_assay_layer <- function(obj, layer_name) {
  tryCatch(
    Seurat::GetAssayData(obj, assay = "RNA", layer = layer_name),
    error = function(e) {
      slot_name <- ifelse(layer_name == "counts", "counts", "data")
      Seurat::GetAssayData(obj, assay = "RNA", slot = slot_name)
    }
  )
}

save_plot <- function(plot, stem, width = 6, height = 4) {
  pdf_path <- file.path("figures", paste0(stem, ".pdf"))
  png_path <- file.path("figures", paste0(stem, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf)
  ggsave(png_path, plot, width = width, height = height, dpi = 300)
  log_msg("Wrote", pdf_path)
  log_msg("Wrote", png_path)
}

extract_pair_test <- function(test_df, pair = c("lineage1", "lineage2")) {
  df <- as.data.frame(test_df) |>
    tibble::rownames_to_column("gene")
  p_cols <- grep("pvalue|pVal|p.value|p.value", colnames(df), ignore.case = TRUE, value = TRUE)
  stat_cols <- grep("wald|stat", colnames(df), ignore.case = TRUE, value = TRUE)
  pair_regex <- paste(pair, collapse = ".*|")
  pair_cols <- grep(paste(pair, collapse = "|"), colnames(df), ignore.case = TRUE, value = TRUE)
  if (length(pair_cols) > 0) {
    p_pair <- intersect(p_cols, pair_cols)
    stat_pair <- intersect(stat_cols, pair_cols)
  } else {
    p_pair <- p_cols
    stat_pair <- stat_cols
  }
  p_col <- if (length(p_pair) > 0) p_pair[1] else if (length(p_cols) > 0) p_cols[1] else NA_character_
  stat_col <- if (length(stat_pair) > 0) stat_pair[1] else if (length(stat_cols) > 0) stat_cols[1] else NA_character_
  out <- df |>
    dplyr::transmute(
      gene,
      waldStat = if (!is.na(stat_col)) .data[[stat_col]] else NA_real_,
      pvalue = if (!is.na(p_col)) .data[[p_col]] else NA_real_
    ) |>
    dplyr::mutate(p_adj_BH = p.adjust(.data$pvalue, method = "BH"))
  attr(out, "p_col") <- p_col
  attr(out, "stat_col") <- stat_col
  out
}

extract_lineage_test <- function(test_df, lineage_cols = paste0("lineage", 1:3)) {
  df <- as.data.frame(test_df) |>
    tibble::rownames_to_column("gene")
  out <- lapply(lineage_cols, function(lin) {
    lin_num <- sub("^lineage", "", lin)
    p_cols <- grep("pvalue|pVal|p.value|p.value", colnames(df), ignore.case = TRUE, value = TRUE)
    stat_cols <- grep("wald|stat", colnames(df), ignore.case = TRUE, value = TRUE)
    lin_cols <- grep(lin, colnames(df), ignore.case = TRUE, value = TRUE)
    if (length(lin_cols) == 0) {
      lin_cols <- grep(paste0("(_", lin_num, "$|_", lin_num, "_|", lin_num, "$)"), colnames(df), ignore.case = TRUE, value = TRUE)
    }
    p_col <- intersect(p_cols, lin_cols)
    stat_col <- intersect(stat_cols, lin_cols)
    if (length(p_col) == 0 && lin == "lineage1" && "pvalue" %in% colnames(df)) p_col <- "pvalue"
    if (length(stat_col) == 0 && lin == "lineage1" && "waldStat" %in% colnames(df)) stat_col <- "waldStat"
    tibble::tibble(
      gene = df$gene,
      lineage = lin,
      waldStat = if (length(stat_col) > 0) df[[stat_col[1]]] else NA_real_,
      pvalue = if (length(p_col) > 0) df[[p_col[1]]] else NA_real_
    )
  }) |>
    dplyr::bind_rows() |>
    dplyr::group_by(.data$lineage) |>
    dplyr::mutate(p_adj_BH = p.adjust(.data$pvalue, method = "BH")) |>
    dplyr::ungroup()
  out
}

load_start <- Sys.time()
obj <- qs2::qs_read(obj_path)
log_msg("Object load seconds:", round(as.numeric(difftime(Sys.time(), load_start, units = "secs")), 2))

meta <- obj@meta.data |>
  tibble::rownames_to_column("cellID") |>
  recode_subtype()

pst <- readr::read_csv(pseudotime_path, show_col_types = FALSE)
weights_long <- readr::read_csv(weight_path, show_col_types = FALSE)
lineages <- paste0("lineage", 1:3)

cells_use <- intersect(colnames(obj), pst$cellID)
obj <- obj[, cells_use]
pst <- pst |> dplyr::filter(.data$cellID %in% cells_use) |> dplyr::arrange(match(.data$cellID, cells_use))
stopifnot(identical(pst$cellID, colnames(obj)))

endpoint_sampling_cells <- readr::read_csv(endpoint_source_path, show_col_types = FALSE) |>
  dplyr::distinct(.data$cellID, .data$endpoint_group) |>
  dplyr::filter(.data$endpoint_group %in% c("MES-I terminal", "MES-V terminal")) |>
  dplyr::pull(.data$cellID) |>
  intersect(colnames(obj))
set.seed(42)
sampled_by_lineage <- pst |>
  dplyr::filter(!.data$cellID %in% endpoint_sampling_cells) |>
  dplyr::group_by(.data$lineage_assignment) |>
  dplyr::group_modify(~ dplyr::slice_sample(.x, n = min(1400, nrow(.x)))) |>
  dplyr::ungroup() |>
  dplyr::pull(.data$cellID)
cells_tradeSeq <- unique(c(endpoint_sampling_cells, sampled_by_lineage))
pst <- pst |> dplyr::filter(.data$cellID %in% cells_tradeSeq) |> dplyr::arrange(match(.data$cellID, cells_tradeSeq))
obj <- obj[, pst$cellID]
log_msg("tradeSeq fitting cells:", ncol(obj), "(all strict MES endpoints retained; remaining cells stratified by assigned lineage)")

pseudotime <- as.matrix(pst[, paste0(lineages, "_pst")])
colnames(pseudotime) <- lineages
cellWeights <- as.matrix(pst[, paste0(lineages, "_weight")])
colnames(cellWeights) <- lineages
pseudotime[is.na(pseudotime) & cellWeights == 0] <- 0
if (anyNA(pseudotime)) {
  pseudotime[is.na(pseudotime)] <- 0
  log_msg("WARNING: remaining NA pseudotime values were set to 0.")
}
cellWeights[is.na(cellWeights)] <- 0

phase4_extra_blacklist <- c("^TRAC$", "^TRBC", "^TRDC$", "^TRGC", "^CD3D$", "^CD3E$", "^CD3G$")
blacklist_patterns <- unique(c(read_blacklist_patterns(blacklist_path), phase4_extra_blacklist))
log_msg("Blacklist pattern count:", length(blacklist_patterns))

key_genes <- unique(c(
  "SOX2", "SOX4", "SOX9", "SOX11", "OLIG1", "OLIG2", "PDGFRA", "PTPRZ1",
  "MBP", "MOG", "PLP1", "MKI67", "TOP2A", "PCNA",
  "CD74", "HLA-DRA", "HLA-DPB1", "HLA-DPA1", "HLA-DRB1", "B2M",
  "RGS5", "ACTA2", "TAGLN", "HIGD1B", "CHI3L1", "CHI3L2", "SERPINE1", "VIM"
))
de_marker_files <- file.path(project_root, "05_恶性细胞分亚群与Neftel对照", "tables", c(
  "10e_Subtype1_vs_rest_wilcoxon_downsampled_ranked_gene_list.csv",
  "10e_Subtype2_vs_rest_wilcoxon_downsampled_ranked_gene_list.csv",
  "10e_Subtype3_vs_rest_wilcoxon_downsampled_ranked_gene_list.csv",
  "10e_Subtype4_vs_rest_wilcoxon_downsampled_ranked_gene_list.csv",
  "10d_S3_vs_S4_passing_markers.csv"
))
de_genes <- character()
for (f in de_marker_files[file.exists(de_marker_files)]) {
  tmp <- readr::read_csv(f, show_col_types = FALSE, n_max = 200)
  gene_col <- intersect(c("gene", "Gene", "features", "feature"), colnames(tmp))
  if (length(gene_col) > 0) de_genes <- c(de_genes, tmp[[gene_col[1]]])
}

hvg <- Seurat::VariableFeatures(obj)
gene_candidates <- unique(c(hvg, de_genes, key_genes))
gene_candidates <- intersect(gene_candidates, rownames(obj))
gene_candidates <- gene_candidates[!is_blacklisted(gene_candidates, blacklist_patterns)]

counts_all <- get_assay_layer(obj, "counts")
counts_sel_sparse <- counts_all[gene_candidates, , drop = FALSE]
detected <- Matrix::rowSums(counts_sel_sparse > 0)
gene_candidates <- rownames(counts_sel_sparse)[detected >= max(20, round(0.002 * ncol(counts_sel_sparse)))]
max_tradeSeq_genes <- 600
if (length(gene_candidates) > max_tradeSeq_genes) {
  keep_key <- intersect(key_genes, gene_candidates)
  keep_rest <- setdiff(gene_candidates, keep_key)
  gene_candidates <- unique(c(keep_key, keep_rest[seq_len(min(length(keep_rest), max_tradeSeq_genes - length(keep_key)))]))
}
log_msg("Selected genes for tradeSeq:", length(gene_candidates))

readr::write_csv(
  tibble::tibble(gene = gene_candidates, is_key_gene = gene_candidates %in% key_genes),
  file.path("tables", "tradeseq_gene_universe.csv")
)

counts <- as.matrix(counts_all[gene_candidates, colnames(obj), drop = FALSE])
mode(counts) <- "integer"

selected_nknots <- 6
readr::write_csv(
  tibble::tibble(
    decision = "evaluateK skipped after initial 270-gene run exceeded 20 min before fitGAM",
    selected_nknots = selected_nknots,
    rationale = "Use fixed nknots=6 to keep Phase 4 tractable; lineage interpretation is driven by cached Slingshot topology and endpoint comparison."
  ),
  file.path("tables", "tradeseq_evaluateK_decision.csv")
)
log_msg("evaluateK skipped; initial 270-gene evaluateK attempt exceeded 20 minutes in prior run.")
log_msg("Selected nknots:", selected_nknots)

fit_start <- Sys.time()
set.seed(42)
fit_path <- file.path("outputs", "tradeseq_fit_hvg_driver_genes.qs2")
if (file.exists(fit_path)) {
  gam <- qs2::qs_read(fit_path)
  log_msg("Loaded existing", fit_path, "instead of refitting.")
} else {
  gam <- tradeSeq::fitGAM(
    counts = counts,
    pseudotime = pseudotime,
    cellWeights = cellWeights,
    nknots = selected_nknots,
    verbose = TRUE,
    BPPARAM = BiocParallel::SerialParam()
  )
  log_msg("fitGAM seconds:", round(as.numeric(difftime(Sys.time(), fit_start, units = "secs")), 2))
  qs2::qs_save(gam, fit_path)
  log_msg("Saved outputs/tradeseq_fit_hvg_driver_genes.qs2")
}

assoc_start <- Sys.time()
assoc_raw <- tradeSeq::associationTest(gam, lineages = TRUE)
assoc_tbl <- extract_lineage_test(assoc_raw, lineages)
assoc_tbl <- assoc_tbl |> dplyr::filter(.data$gene %in% gene_candidates)
log_msg("associationTest seconds:", round(as.numeric(difftime(Sys.time(), assoc_start, units = "secs")), 2))

start_end_start <- Sys.time()
start_end_raw <- tradeSeq::startVsEndTest(gam, lineages = TRUE)
start_end_tbl <- extract_lineage_test(start_end_raw, lineages)
start_end_tbl <- start_end_tbl |> dplyr::filter(.data$gene %in% gene_candidates)
log_msg("startVsEndTest seconds:", round(as.numeric(difftime(Sys.time(), start_end_start, units = "secs")), 2))

diff_start <- Sys.time()
diff_end_raw <- tradeSeq::diffEndTest(gam)
diff_end_tbl <- extract_pair_test(diff_end_raw, c("lineage1", "lineage2"))
diff_end_tbl <- diff_end_tbl |> dplyr::filter(.data$gene %in% gene_candidates)
log_msg("diffEndTest seconds:", round(as.numeric(difftime(Sys.time(), diff_start, units = "secs")), 2))
log_msg("diffEnd extracted p column:", attr(diff_end_tbl, "p_col"))
log_msg("diffEnd extracted stat column:", attr(diff_end_tbl, "stat_col"))

norm_data <- get_assay_layer(obj, "data")
norm_sel <- norm_data[gene_candidates, colnames(obj), drop = FALSE]

assigned_long <- pst |>
  dplyr::select("cellID", "subtype_short", "lineage_assignment", dplyr::all_of(paste0(lineages, "_pst"))) |>
  tidyr::pivot_longer(cols = dplyr::all_of(paste0(lineages, "_pst")), names_to = "lineage_col", values_to = "pseudotime") |>
  dplyr::mutate(lineage = sub("_pst$", "", .data$lineage_col)) |>
  dplyr::left_join(
    weights_long |> dplyr::select("cellID", "lineage_id", "lineage_weight"),
    by = c("cellID" = "cellID", "lineage" = "lineage_id")
  ) |>
  dplyr::filter(.data$lineage_assignment == .data$lineage, .data$lineage_weight > 0.5, is.finite(.data$pseudotime))

make_binned_expr <- function(lineage_id, genes) {
  cells <- assigned_long |> dplyr::filter(.data$lineage == lineage_id)
  if (nrow(cells) == 0) return(tibble::tibble())
  cells <- cells |>
    dplyr::mutate(bin = dplyr::ntile(.data$pseudotime, 20))
  genes <- intersect(genes, rownames(norm_sel))
  mat <- norm_sel[genes, cells$cellID, drop = FALSE]
  out <- lapply(sort(unique(cells$bin)), function(b) {
    b_cells <- cells$cellID[cells$bin == b]
    tibble::tibble(
      gene = genes,
      lineage = lineage_id,
      bin = b,
      pseudotime_mean = mean(cells$pseudotime[cells$bin == b], na.rm = TRUE),
      expr_mean = Matrix::rowMeans(mat[, b_cells, drop = FALSE])
    )
  }) |>
    dplyr::bind_rows()
  out
}

classify_patterns <- function(binned_df) {
  binned_df |>
    dplyr::group_by(.data$lineage, .data$gene) |>
    dplyr::summarise(
      early = mean(.data$expr_mean[.data$bin <= 3], na.rm = TRUE),
      late = mean(.data$expr_mean[.data$bin >= 18], na.rm = TRUE),
      mid = mean(.data$expr_mean[.data$bin > 3 & .data$bin < 18], na.rm = TRUE),
      max_bin = .data$bin[which.max(.data$expr_mean)][1],
      .groups = "drop"
    ) |>
    dplyr::mutate(
      pattern = dplyr::case_when(
        .data$late - .data$early > 0.25 ~ "late_up",
        .data$early - .data$late > 0.25 ~ "early_high",
        .data$max_bin > 4 & .data$max_bin < 17 & .data$mid > pmax(.data$early, .data$late) + 0.15 ~ "transient",
        TRUE ~ "flat_or_subtle"
      )
    ) |>
    dplyr::select("gene", "lineage", "pattern", "early", "late", "mid", "max_bin")
}

top_assoc_for_patterns <- assoc_tbl |>
  dplyr::filter(!is.na(.data$pvalue)) |>
  dplyr::group_by(.data$lineage) |>
  dplyr::slice_min(.data$p_adj_BH, n = 120, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::pull(.data$gene) |>
  unique()
binned_for_patterns <- lapply(lineages, make_binned_expr, genes = top_assoc_for_patterns) |>
  dplyr::bind_rows()
pattern_tbl <- classify_patterns(binned_for_patterns)

assoc_tbl <- assoc_tbl |>
  dplyr::left_join(pattern_tbl |> dplyr::select("gene", "lineage", "pattern"), by = c("gene", "lineage"))
start_end_tbl <- start_end_tbl |>
  dplyr::left_join(pattern_tbl |> dplyr::select("gene", "lineage", "pattern"), by = c("gene", "lineage"))

endpoint_cells <- readr::read_csv(endpoint_source_path, show_col_types = FALSE) |>
  dplyr::distinct(.data$cellID, .data$endpoint_group) |>
  dplyr::filter(.data$endpoint_group %in% c("MES-I terminal", "MES-V terminal"))
endpoint_expr <- norm_data[intersect(gene_candidates, rownames(norm_data)), endpoint_cells$cellID, drop = FALSE]
mes_i_cells <- endpoint_cells$cellID[endpoint_cells$endpoint_group == "MES-I terminal"]
mes_v_cells <- endpoint_cells$cellID[endpoint_cells$endpoint_group == "MES-V terminal"]
endpoint_fc <- tibble::tibble(
  gene = rownames(endpoint_expr),
  mean_MESI = Matrix::rowMeans(endpoint_expr[, mes_i_cells, drop = FALSE]),
  mean_MESV = Matrix::rowMeans(endpoint_expr[, mes_v_cells, drop = FALSE])
) |>
  dplyr::mutate(log2FC_MESV_vs_MESI = log2((.data$mean_MESV + 0.05) / (.data$mean_MESI + 0.05)))
diff_end_tbl <- diff_end_tbl |>
  dplyr::left_join(endpoint_fc, by = "gene") |>
  dplyr::mutate(
    direction = dplyr::case_when(
      .data$log2FC_MESV_vs_MESI > 0 ~ "MES-V high",
      .data$log2FC_MESV_vs_MESI < 0 ~ "MES-I high",
      TRUE ~ "no_direction"
    ),
    contamination_sensitive_marker = .data$gene %in% c("CD74", "HLA-DRA", "HLA-DPB1", "HLA-DPA1", "HLA-DRB1", "B2M", "RGS5", "ACTA2", "TAGLN", "HIGD1B")
  ) |>
  dplyr::arrange(.data$p_adj_BH)

readr::write_csv(assoc_tbl, file.path("tables", "tradeseq_association_per_lineage.csv"))
readr::write_csv(start_end_tbl, file.path("tables", "tradeseq_startVsEnd_per_lineage.csv"))
readr::write_csv(diff_end_tbl, file.path("tables", "tradeseq_diffEnd_MESI_vs_MESV.csv"))
readr::write_csv(pattern_tbl, file.path("tables", "driver_gene_clusters.csv"))
readr::write_csv(binned_for_patterns, file.path("figures/source_data", "08_binned_expression_for_patterns.csv"))

enrich_input <- dplyr::bind_rows(
  assoc_tbl |>
    dplyr::filter(.data$lineage == "lineage3", .data$p_adj_BH < 0.05) |>
    dplyr::arrange(.data$p_adj_BH) |>
    dplyr::slice_head(n = 150) |>
    dplyr::mutate(gene_set = "L3_association"),
  diff_end_tbl |>
    dplyr::filter(.data$p_adj_BH < 0.05, .data$direction == "MES-I high") |>
    dplyr::arrange(.data$p_adj_BH) |>
    dplyr::slice_head(n = 150) |>
    dplyr::mutate(lineage = "diffEnd", gene_set = "MES-I_endpoint_high"),
  diff_end_tbl |>
    dplyr::filter(.data$p_adj_BH < 0.05, .data$direction == "MES-V high") |>
    dplyr::arrange(.data$p_adj_BH) |>
    dplyr::slice_head(n = 150) |>
    dplyr::mutate(lineage = "diffEnd", gene_set = "MES-V_endpoint_high")
) |>
  dplyr::select("gene", "gene_set") |>
  dplyr::distinct()

enrichment_tbl <- tibble::tibble()
if (nrow(enrich_input) > 0) {
  gene_map <- clusterProfiler::bitr(unique(enrich_input$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
  universe_map <- clusterProfiler::bitr(gene_candidates, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
  enrich_universe <- unique(universe_map$ENTREZID)
  readr::write_csv(
    tibble::tibble(
      background_type = "tradeSeq_gene_universe",
      n_symbol = length(unique(gene_candidates)),
      n_entrez = length(enrich_universe)
    ),
    file.path("tables", "driver_enrichment_background.csv")
  )
  log_msg("GO enrichment background: tradeSeq gene universe; SYMBOL n =", length(unique(gene_candidates)), "; ENTREZ n =", length(enrich_universe))
  for (gs in unique(enrich_input$gene_set)) {
    entrez <- gene_map$ENTREZID[match(enrich_input$gene[enrich_input$gene_set == gs], gene_map$SYMBOL)]
    entrez <- unique(entrez[!is.na(entrez)])
    if (length(entrez) >= 10) {
      ego <- tryCatch(
        clusterProfiler::enrichGO(gene = entrez, universe = enrich_universe, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.2, readable = TRUE),
        error = function(e) NULL
      )
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        enrichment_tbl <- dplyr::bind_rows(enrichment_tbl, as.data.frame(ego) |> dplyr::mutate(gene_set = gs))
      }
    }
  }
}
readr::write_csv(enrichment_tbl, file.path("tables", "driver_enrichment_per_cluster.csv"))

top_l3_genes <- assoc_tbl |>
  dplyr::filter(.data$lineage == "lineage3", !is.na(.data$p_adj_BH)) |>
  dplyr::arrange(.data$p_adj_BH) |>
  dplyr::slice_head(n = 35) |>
  dplyr::pull(.data$gene)
l3_heat <- make_binned_expr("lineage3", top_l3_genes) |>
  dplyr::group_by(.data$gene) |>
  dplyr::mutate(z = as.numeric(scale(.data$expr_mean))) |>
  dplyr::ungroup() |>
  dplyr::mutate(gene = factor(.data$gene, levels = rev(top_l3_genes)))
readr::write_csv(l3_heat, file.path("figures/source_data", "08_L3_driver_heatmap.csv"))
p_l3 <- ggplot(l3_heat, aes(bin, gene, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Z") +
  labs(title = "L3 driver genes along NPC-P to OPC-M pseudotime", x = "Pseudotime bin", y = NULL) +
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 5.2))
save_plot(p_l3, "08_L3_driver_heatmap", width = 5.2, height = 6.2)

volcano_data <- diff_end_tbl |>
  dplyr::mutate(
    neglog10_BH = -log10(pmax(.data$p_adj_BH, 1e-300)),
    neglog10_BH_capped = pmin(.data$neglog10_BH, 60),
    label_gene = .data$gene %in% c("CD74", "HLA-DRA", "HLA-DPB1", "HLA-DPA1", "HLA-DRB1", "RGS5", "ACTA2", "TAGLN", "HIGD1B", "SOX9", "CHI3L1")
  )
readr::write_csv(volcano_data, file.path("figures/source_data", "08_MES_diffEnd_volcano.csv"))
p_volcano <- ggplot(volcano_data, aes(log2FC_MESV_vs_MESI, neglog10_BH_capped)) +
  geom_point(aes(color = direction), alpha = 0.65, size = 0.8) +
  geom_text(data = dplyr::filter(volcano_data, .data$label_gene), aes(label = gene), size = 2.1, vjust = -0.6, check_overlap = TRUE) +
  scale_color_manual(values = c("MES-I high" = "#BC3C29", "MES-V high" = "#20854E", "no_direction" = "grey65"), name = NULL) +
  labs(
    title = "MES-I vs MES-V terminal diffEnd genes",
    subtitle = "HLA/CD74 and vascular markers are retained with Fig 07 validation caveat",
    x = "Endpoint log2FC (MES-V / MES-I)",
    y = "-log10(BH), capped at 60"
  ) +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
save_plot(p_volcano, "08_MES_diffEnd_volcano", width = 5.4, height = 4.0)

if (nrow(enrichment_tbl) > 0) {
  enrich_plot <- enrichment_tbl |>
    dplyr::group_by(.data$gene_set) |>
    dplyr::arrange(.data$p.adjust, .by_group = TRUE) |>
    dplyr::slice_head(n = 6) |>
    dplyr::ungroup() |>
    dplyr::mutate(Description = stringr::str_trunc(.data$Description, 55), Description = factor(.data$Description, levels = rev(unique(.data$Description))))
  readr::write_csv(enrich_plot, file.path("figures/source_data", "08_driver_enrichment_dotplot.csv"))
  p_enrich <- ggplot(enrich_plot, aes(-log10(p.adjust), Description, size = Count, color = gene_set)) +
    geom_point(alpha = 0.85) +
    facet_wrap(~gene_set, scales = "free_y") +
    labs(title = "Driver gene GO enrichment", x = "-log10(BH)", y = NULL, color = NULL) +
    theme_bw(base_size = 7) +
    theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 5.4), legend.position = "bottom")
  save_plot(p_enrich, "08_driver_enrichment_dotplot", width = 8.0, height = 4.5)
} else {
  log_msg("No enrichment terms returned.")
}

l1_l2_support <- assoc_tbl |>
  dplyr::filter(.data$lineage %in% c("lineage1", "lineage2"), !is.na(.data$p_adj_BH)) |>
  dplyr::group_by(.data$lineage) |>
  dplyr::arrange(.data$p_adj_BH, .by_group = TRUE) |>
  dplyr::slice_head(n = 20) |>
  dplyr::ungroup()
readr::write_csv(l1_l2_support, file.path("figures/source_data", "08_L1L2_association_caveat.csv"))
p_l1l2 <- ggplot(l1_l2_support, aes(-log10(pmax(p_adj_BH, 1e-300)), reorder(gene, -log10(pmax(p_adj_BH, 1e-300))), fill = lineage)) +
  geom_col(width = 0.65) +
  facet_wrap(~lineage, scales = "free_y") +
  scale_fill_manual(values = c("lineage1" = "#BC3C29", "lineage2" = "#20854E"), guide = "none") +
  labs(title = "L1/L2 association genes (supporting only)", subtitle = "Fine-ordering caveat applies to MES-like lineages", x = "-log10(BH)", y = NULL) +
  theme_bw(base_size = 7) +
  theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 5.5))
save_plot(p_l1l2, "08_L1L2_association_caveat", width = 6.0, height = 5.0)

if (requireNamespace("patchwork", quietly = TRUE)) {
  preview <- (p_l3 | p_volcano) / (if (exists("p_enrich")) p_enrich else p_l1l2)
  ggsave(file.path("figures", "08_driver_genes_preview.png"), preview, width = 10.5, height = 8.5, dpi = 200)
  log_msg("Wrote figures/08_driver_genes_preview.png")
}

log_msg("Top diffEnd genes:")
capture.output(print(diff_end_tbl |> dplyr::slice_head(n = 25))) |> paste(collapse = "\n") |> log_msg()
log_msg("Top L3 association genes:")
capture.output(print(assoc_tbl |> dplyr::filter(.data$lineage == "lineage3") |> dplyr::arrange(.data$p_adj_BH) |> dplyr::slice_head(n = 25))) |> paste(collapse = "\n") |> log_msg()
writeLines(capture.output(sessionInfo()), file.path("logs", "08_session_info.txt"))
log_msg("Session info: logs/08_session_info.txt")
log_msg("STOP: Phase 4 tradeSeq complete. Await audit.")
