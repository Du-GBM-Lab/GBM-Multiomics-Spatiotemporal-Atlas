suppressPackageStartupMessages({
  library(qs2)
  library(Matrix)
  library(Seurat)
  library(data.table)
})

set.seed(1)

root <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
big_path <- file.path(root, "outputs/R9_A2_RCTD_weights_allslices_long.qs2")
path_mal <- file.path(root, "data/reference/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2")
path_markers <- file.path(root, "data/reference/sc_subtype_markers.rds")
out_dir <- file.path(root, "outputs")
tab_dir <- file.path(root, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

mv <- "Subtype3"
mi <- "Subtype4"
thresholds <- c(0.05, 0.03)

dominance_stats <- function(dt, thr) {
  V <- dt[[mv]]
  I <- dt[[mi]]
  mes_pos <- (V > thr) | (I > thr)
  dom_v <- V > thr & V > 2 * I
  dom_i <- I > thr & I > 2 * V
  both <- V > thr & I > thr & pmax(V, I) <= 2 * pmin(V, I)
  other_mes <- mes_pos & !(dom_v | dom_i | both)

  data.table(
    threshold = thr,
    n_total = length(V),
    n_mes_pos = sum(mes_pos),
    pct_mes_pos_all = sum(mes_pos) / length(V),
    cor_all_spearman = suppressWarnings(cor(V, I, method = "spearman")),
    cor_mes_pos_spearman = if (sum(mes_pos) > 50) {
      suppressWarnings(cor(V[mes_pos], I[mes_pos], method = "spearman"))
    } else NA_real_,
    n_dom_v = sum(dom_v),
    n_dom_i = sum(dom_i),
    n_both = sum(both),
    n_other_mes = sum(other_mes),
    pct_dom_v_mespos = if (sum(mes_pos)) sum(dom_v) / sum(mes_pos) else NA_real_,
    pct_dom_i_mespos = if (sum(mes_pos)) sum(dom_i) / sum(mes_pos) else NA_real_,
    pct_both_mespos = if (sum(mes_pos)) sum(both) / sum(mes_pos) else NA_real_,
    pct_other_mespos = if (sum(mes_pos)) sum(other_mes) / sum(mes_pos) else NA_real_
  )
}

per_slice_stats <- function(dt, thr) {
  ans <- dt[, {
    V <- get(mv)
    I <- get(mi)
    mes_pos <- (V > thr) | (I > thr)
    dom_v <- V > thr & V > 2 * I
    dom_i <- I > thr & I > 2 * V
    both <- V > thr & I > thr & pmax(V, I) <= 2 * pmin(V, I)
    .(
      n_spots = .N,
      mean_mesv = mean(V),
      mean_mesi = mean(I),
      max_mesv = max(V),
      max_mesi = max(I),
      n_mes_pos = sum(mes_pos),
      n_dom_v = sum(dom_v),
      n_dom_i = sum(dom_i),
      n_both = sum(both),
      cor_vi_spearman = if (.N > 2) suppressWarnings(cor(V, I, method = "spearman")) else NA_real_,
      cor_vi_mespos_spearman = if (sum(mes_pos) > 10) {
        suppressWarnings(cor(V[mes_pos], I[mes_pos], method = "spearman"))
      } else NA_real_
    )
  }, by = slice]
  ans[, threshold := thr]
  setcolorder(ans, c("threshold", "slice"))
  ans[]
}

big <- as.data.table(qs2::qs_read(big_path))
stopifnot(all(c("spot_id", "slice", "x", "y", mv, mi) %in% colnames(big)))

cat("== A2 weights:", nrow(big), "spots x", ncol(big), "columns\n")
cat("== MES columns:", mv, mi, "\n")

V <- big[[mv]]
I <- big[[mi]]
cat("\n== MES-V / MES-I distribution:\n")
dist_tab <- data.table(
  subtype = c("MES_V_Subtype3", "MES_I_Subtype4"),
  mean = c(mean(V), mean(I)),
  q50 = c(quantile(V, 0.5), quantile(I, 0.5)),
  q75 = c(quantile(V, 0.75), quantile(I, 0.75)),
  q90 = c(quantile(V, 0.9), quantile(I, 0.9)),
  q95 = c(quantile(V, 0.95), quantile(I, 0.95)),
  q99 = c(quantile(V, 0.99), quantile(I, 0.99)),
  max = c(max(V), max(I))
)
print(dist_tab)

dom_tab <- rbindlist(lapply(thresholds, function(thr) dominance_stats(big, thr)))
cat("\n== global separability stats:\n")
print(dom_tab)

slice_tab <- rbindlist(lapply(thresholds, function(thr) per_slice_stats(big, thr)))
cat("\n== per-slice stats:\n")
print(slice_tab)

markers <- readRDS(path_markers)
stopifnot(all(c("MES_V", "MES_I") %in% names(markers)))
mesv_markers <- unique(as.character(markers$MES_V))
mesi_markers <- unique(as.character(markers$MES_I))
sig_overlap <- intersect(mesv_markers, mesi_markers)

cat("\n== marker overlap:\n")
cat("MES_V markers:", length(mesv_markers), "| MES_I markers:", length(mesi_markers),
    "| overlap:", length(sig_overlap), "\n")
if (length(sig_overlap)) print(sig_overlap)

mal <- qs2::qs_read(path_mal)
DefaultAssay(mal) <- "RNA"
mal <- NormalizeData(mal, verbose = FALSE)
mesv_use <- intersect(mesv_markers, rownames(mal))
mesi_use <- intersect(mesi_markers, rownames(mal))
cat("== markers present in malignant object: MES_V", length(mesv_use),
    "| MES_I", length(mesi_use), "\n")
mal <- AddModuleScore(mal, features = list(mesv_use), name = "mvS", seed = 1)
mal <- AddModuleScore(mal, features = list(mesi_use), name = "miS", seed = 1)
sig_cor_all <- suppressWarnings(cor(mal$mvS1, mal$miS1, method = "spearman"))

subtype_col <- if ("subtype_k4" %in% colnames(mal@meta.data)) "subtype_k4" else NA_character_
if (!is.na(subtype_col)) {
  sig_cor_by_subtype <- data.table(
    subtype = names(tapply(mal$mvS1, mal[[subtype_col]][, 1], length)),
    n = as.integer(tapply(mal$mvS1, mal[[subtype_col]][, 1], length)),
    cor_mv_mi = as.numeric(tapply(seq_along(mal$mvS1), mal[[subtype_col]][, 1], function(idx) {
      suppressWarnings(cor(mal$mvS1[idx], mal$miS1[idx], method = "spearman"))
    })),
    mean_mv_score = as.numeric(tapply(mal$mvS1, mal[[subtype_col]][, 1], mean)),
    mean_mi_score = as.numeric(tapply(mal$miS1, mal[[subtype_col]][, 1], mean))
  )
} else {
  sig_cor_by_subtype <- data.table()
}

cat("== reference malignant module-score Spearman(MES_V vs MES_I):",
    round(sig_cor_all, 3), "\n")
if (nrow(sig_cor_by_subtype)) {
  cat("== by subtype module-score summary:\n")
  print(sig_cor_by_subtype)
}

data.table::fwrite(dist_tab, file.path(tab_dir, "R9_B_MESV_MESI_distribution.csv"))
data.table::fwrite(dom_tab, file.path(tab_dir, "R9_B_MESV_MESI_global_separability.csv"))
data.table::fwrite(slice_tab, file.path(tab_dir, "R9_B_MESV_MESI_perslice_separability.csv"))
if (nrow(sig_cor_by_subtype)) {
  data.table::fwrite(sig_cor_by_subtype, file.path(tab_dir, "R9_B_reference_signature_correlation_by_subtype.csv"))
}
data.table::fwrite(
  data.table(overlap_gene = sig_overlap),
  file.path(tab_dir, "R9_B_MESV_MESI_marker_overlap.csv")
)

qs2::qs_save(
  list(
    distribution = dist_tab,
    global = dom_tab,
    perslice = slice_tab,
    marker_overlap = sig_overlap,
    mesv_markers_present = mesv_use,
    mesi_markers_present = mesi_use,
    signature_cor_all = sig_cor_all,
    signature_cor_by_subtype = sig_cor_by_subtype
  ),
  file.path(out_dir, "R9_B_MESV_MESI_separability.qs2")
)

cat("\n[STOP B] 回报:MES-V/MES-I distribution; global separability at 0.05 and 0.03; per-slice stats; marker overlap; reference module-score correlation.\n")
