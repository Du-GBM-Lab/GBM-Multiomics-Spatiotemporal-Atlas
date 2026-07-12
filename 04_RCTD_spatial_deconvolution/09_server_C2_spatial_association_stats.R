#!/usr/bin/env Rscript
# =============================================================================
# R9 | C2: RCTD weights spot-level co-enrichment statistics
# Input : A2 spot x celltype RCTD weights, already row-normalized.
# Output: global and per-slice association statistics for MES states vs
#         myeloid / vascular weights.
# STOP  : quantify same-spot deconvolution co-enrichment only. No PAN/niche,
#         no PLAU/PLAUR/FOSL1, no causal or physical-contact inference.
# =============================================================================

suppressPackageStartupMessages({
  library(qs2)
  library(data.table)
})

base_dir <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
in_path  <- file.path(base_dir, "outputs/R9_A2_RCTD_weights_allslices_long.qs2")
out_dir  <- file.path(base_dir, "outputs/R9_C2_spatial_association")
tab_dir  <- file.path(base_dir, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

message("== reading: ", in_path)
big <- as.data.table(qs2::qs_read(in_path))

required <- c("spot_id", "slice", "x", "y",
              "Subtype3", "Subtype4",
              "Macrophages", "Microglial", "Monocytes",
              "Endothelial", "Mural cells")
miss <- setdiff(required, names(big))
if (length(miss)) stop("missing required columns: ", paste(miss, collapse = ", "))

big[, `MES-V`       := Subtype3]
big[, `MES-I`       := Subtype4]
big[, `MES-lineage` := Subtype3 + Subtype4]
big[, Myeloid       := Macrophages + Microglial + Monocytes]
big[, Vascular      := Endothelial + `Mural cells`]

pairs <- data.table(
  mes_feature = c("MES-V", "MES-V", "MES-I", "MES-I",
                  "MES-lineage", "MES-lineage",
                  "MES-I", "MES-I", "MES-I",
                  "MES-V", "MES-V"),
  context_feature = c("Vascular", "Myeloid", "Vascular", "Myeloid",
                      "Vascular", "Myeloid",
                      "Macrophages", "Microglial", "Monocytes",
                      "Endothelial", "Mural cells")
)

safe_cor <- function(x, y) {
  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman"))
}

safe_fisher <- function(high_x, high_y) {
  tab <- table(factor(high_x, levels = c(FALSE, TRUE)),
               factor(high_y, levels = c(FALSE, TRUE)))
  ft <- tryCatch(fisher.test(tab), error = function(e) NULL)
  if (is.null(ft)) return(list(or = NA_real_, p = NA_real_))
  list(or = unname(ft$estimate), p = ft$p.value)
}

one_pair_stats <- function(dt, mes_feature, context_feature, group_label) {
  x <- dt[[mes_feature]]
  y <- dt[[context_feature]]
  n <- length(x)

  qx90 <- as.numeric(quantile(x, 0.90, na.rm = TRUE, names = FALSE))
  qy90 <- as.numeric(quantile(y, 0.90, na.rm = TRUE, names = FALSE))
  high_x <- x >= qx90
  high_y <- y >= qy90
  ft <- safe_fisher(high_x, high_y)

  xpos05 <- x > 0.05
  xpos03 <- x > 0.03

  data.table(
    group = group_label,
    n_spots = n,
    mes_feature = mes_feature,
    context_feature = context_feature,
    mean_mes = mean(x, na.rm = TRUE),
    mean_context = mean(y, na.rm = TRUE),
    cor_all_spearman = safe_cor(x, y),
    n_mes_pos_005 = sum(xpos05, na.rm = TRUE),
    cor_mes_pos_005 = safe_cor(x[xpos05], y[xpos05]),
    n_mes_pos_003 = sum(xpos03, na.rm = TRUE),
    cor_mes_pos_003 = safe_cor(x[xpos03], y[xpos03]),
    q90_mes = qx90,
    q90_context = qy90,
    n_high_mes = sum(high_x, na.rm = TRUE),
    n_high_context = sum(high_y, na.rm = TRUE),
    n_top_decile_overlap = sum(high_x & high_y, na.rm = TRUE),
    expected_overlap = sum(high_x, na.rm = TRUE) * sum(high_y, na.rm = TRUE) / n,
    overlap_over_expected = sum(high_x & high_y, na.rm = TRUE) /
      (sum(high_x, na.rm = TRUE) * sum(high_y, na.rm = TRUE) / n),
    fisher_or = ft$or,
    fisher_p = ft$p
  )
}

global_stats <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
  one_pair_stats(big, pairs$mes_feature[i], pairs$context_feature[i], "all_slices")
}))

per_slice_stats <- rbindlist(lapply(unique(big$slice), function(sl) {
  dt <- big[slice == sl]
  rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    one_pair_stats(dt, pairs$mes_feature[i], pairs$context_feature[i], sl)
  }))
}))

## Per-slice stability summaries for primary six pairs.
primary <- pairs[mes_feature %in% c("MES-V", "MES-I", "MES-lineage") &
                   context_feature %in% c("Vascular", "Myeloid")]
primary_key <- paste(primary$mes_feature, primary$context_feature, sep = "__")
per_slice_stats[, pair := paste(mes_feature, context_feature, sep = "__")]

stability <- per_slice_stats[pair %in% primary_key, .(
  n_slices = .N,
  n_cor_positive = sum(cor_all_spearman > 0, na.rm = TRUE),
  n_cor_negative = sum(cor_all_spearman < 0, na.rm = TRUE),
  median_cor_all = median(cor_all_spearman, na.rm = TRUE),
  q25_cor_all = quantile(cor_all_spearman, 0.25, na.rm = TRUE),
  q75_cor_all = quantile(cor_all_spearman, 0.75, na.rm = TRUE),
  n_or_gt1 = sum(fisher_or > 1, na.rm = TRUE),
  median_fisher_or = median(fisher_or, na.rm = TRUE),
  q25_fisher_or = quantile(fisher_or, 0.25, na.rm = TRUE),
  q75_fisher_or = quantile(fisher_or, 0.75, na.rm = TRUE)
), by = .(mes_feature, context_feature)]

qs2::qs_save(
  list(
    source = "R9_A2_RCTD_weights_allslices_long.qs2",
    derived_columns = c("MES-V=Subtype3", "MES-I=Subtype4",
                        "MES-lineage=Subtype3+Subtype4",
                        "Myeloid=Macrophages+Microglial+Monocytes",
                        "Vascular=Endothelial+Mural cells"),
    pairs = pairs,
    global_stats = global_stats,
    per_slice_stats = per_slice_stats,
    stability = stability
  ),
  file.path(out_dir, "R9_C2_spatial_association_stats.qs2")
)

fwrite(global_stats, file.path(tab_dir, "R9_C2_global_spot_association.csv"))
fwrite(per_slice_stats, file.path(tab_dir, "R9_C2_perslice_spot_association.csv"))
fwrite(stability, file.path(tab_dir, "R9_C2_perslice_stability_summary.csv"))

cat("\n== C2 global primary pairs:\n")
print(global_stats[mes_feature %in% c("MES-V", "MES-I", "MES-lineage") &
                     context_feature %in% c("Vascular", "Myeloid"),
                   .(mes_feature, context_feature, n_spots,
                     cor_all_spearman, cor_mes_pos_005,
                     n_top_decile_overlap, expected_overlap,
                     overlap_over_expected, fisher_or, fisher_p)])

cat("\n== C2 per-slice stability primary pairs:\n")
print(stability)

cat("\n[STOP C2] Report global associations, per-slice stability, and whether any\n",
    "MES state shows stable same-spot co-enrichment with myeloid/vascular weights.\n",
    "Do not infer physical contact or causality.\n", sep = "")
