#!/usr/bin/env Rscript
# =============================================================================
# R9 | Stage 3: IVY/hypoxia scoring + redundancy test
# Purpose:
#   Score locked IVY-like and hypoxia signatures with one unified method, then
#   run per-slice redundancy correlations before any proximity analysis.
#
# Method:
#   Seurat AddModuleScore on SCT data slot, nbin=24, ctrl=100, seed=1.
#
# Outputs:
#   - per-spot signature scores
#   - score parameter/source data
#   - per-slice redundancy correlations
#   - cross-slice redundancy summary
#
# STOP:
#   No MES-V proximity on these maps is run in this script.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(qs2)
})

set.seed(1)

base_dir <- getwd()
st_path <- "<DATA_ROOT>/项目/分型/分型代码/0.对象/5.ST_merge.rds"
weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
sig_path <- file.path(base_dir, "tables/R9_stage3_signature_sets/R9_stage3_locked_signature_genes_IVY_hypoxia.csv")
out_dir <- file.path(base_dir, "tables/R9_stage3_signature_scores")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

signature_order <- c("CT", "CTmvp", "CTpan", "LE",
                     "BUFFA_HYPOXIA_METAGENE", "HALLMARK_HYPOXIA")

## ---- Load inputs -----------------------------------------------------------
st <- readRDS(st_path)
DefaultAssay(st) <- "SCT"
weights <- as.data.table(qs2::qs_read(weights_path))
sig <- fread(sig_path)
stopifnot(all(signature_order %in% unique(sig$signature)))

sct_genes <- rownames(st[["SCT"]])
features <- lapply(signature_order, function(s) {
  unique(sig[signature == s, gene])
})
names(features) <- signature_order
features_use <- lapply(features, intersect, y = sct_genes)

coverage <- rbindlist(lapply(signature_order, function(s) {
  data.table(
    signature = s,
    n_input = length(features[[s]]),
    n_used = length(features_use[[s]]),
    frac_used = length(features_use[[s]]) / length(features[[s]]),
    missing = paste(setdiff(features[[s]], features_use[[s]]), collapse = ";")
  )
}))
if (any(coverage$frac_used < 0.80)) {
  stop("One or more signatures have <80% SCT gene coverage; inspect coverage table.")
}

## ---- Score signatures using one method ------------------------------------
score_prefix <- "R9sig"
st <- AddModuleScore(
  object = st,
  features = features_use,
  assay = "SCT",
  name = score_prefix,
  seed = 1,
  nbin = 24,
  ctrl = 100,
  slot = "data"
)
score_cols <- paste0(score_prefix, seq_along(signature_order))
stopifnot(all(score_cols %in% colnames(st@meta.data)))

score_col_map <- data.table(
  signature = signature_order,
  score_col = score_cols,
  output_col = c("IVY_CT", "IVY_CTmvp", "IVY_CTpan", "IVY_LE",
                 "Hypoxia_Buffa", "Hypoxia_Hallmark")
)

meta <- as.data.table(st@meta.data, keep.rownames = "spot_id")
score_dt <- meta[, c("spot_id", score_cols), with = FALSE]
setnames(score_dt, score_cols, score_col_map$output_col)

## Add slice/x/y and RCTD niche context.
context <- weights[, .(
  spot_id, slice, image, x, y,
  vascular = Endothelial + `Mural cells`,
  myeloid = Macrophages + Microglial + Monocytes,
  microglia = Microglial,
  macrophage_monocyte = Macrophages + Monocytes,
  neuron_control = Neurons,
  `MES-lineage` = Subtype3 + Subtype4,
  `MES-V` = Subtype3,
  `MES-I` = Subtype4
)]
score_dt <- merge(context, score_dt, by = "spot_id", all.x = TRUE)
stopifnot(nrow(score_dt) == nrow(weights))

## ---- Redundancy correlations ----------------------------------------------
score_vars <- c("IVY_CT", "IVY_CTmvp", "IVY_CTpan", "IVY_LE",
                "Hypoxia_Buffa", "Hypoxia_Hallmark")
context_vars <- c("vascular", "myeloid", "microglia",
                  "macrophage_monocyte", "neuron_control")
all_vars <- c(score_vars, context_vars)

pair_grid <- CJ(var1 = all_vars, var2 = all_vars)
pair_grid <- pair_grid[var1 < var2]

cor_one <- function(x, y) {
  suppressWarnings(cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
}

per_slice <- score_dt[, {
  lapply(seq_len(nrow(pair_grid)), function(i) {
    v1 <- pair_grid$var1[i]
    v2 <- pair_grid$var2[i]
    data.table(
      var1 = v1,
      var2 = v2,
      rho = cor_one(get(v1), get(v2)),
      n_spots = .N
    )
  }) |> rbindlist()
}, by = slice]

wilcox_p <- function(x, alternative = "two.sided") {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, alternative = alternative, exact = FALSE)$p.value)
}

summary <- per_slice[, .(
  n_slices = .N,
  median_rho = median(rho, na.rm = TRUE),
  q25_rho = as.numeric(quantile(rho, 0.25, na.rm = TRUE)),
  q75_rho = as.numeric(quantile(rho, 0.75, na.rm = TRUE)),
  min_rho = min(rho, na.rm = TRUE),
  max_rho = max(rho, na.rm = TRUE),
  n_positive = sum(rho > 0, na.rm = TRUE),
  n_negative = sum(rho < 0, na.rm = TRUE),
  wilcox_p_two_sided = wilcox_p(rho, "two.sided"),
  wilcox_p_greater = wilcox_p(rho, "greater"),
  wilcox_p_less = wilcox_p(rho, "less")
), by = .(var1, var2)]

## Mark focus pairs requested in the spec.
focus_pairs <- data.table(
  pair_label = c("CTmvp_vs_vascular", "Buffa_vs_vascular",
                 "CTpan_vs_Buffa", "Buffa_vs_Hallmark"),
  var_a = c("IVY_CTmvp", "Hypoxia_Buffa", "IVY_CTpan", "Hypoxia_Buffa"),
  var_b = c("vascular", "vascular", "Hypoxia_Buffa", "Hypoxia_Hallmark")
)
focus_pairs[, `:=`(
  var1 = pmin(var_a, var_b),
  var2 = pmax(var_a, var_b)
)]
focus_summary <- merge(focus_pairs[, .(pair_label, var1, var2)],
                       summary, by = c("var1", "var2"), all.x = TRUE)

params <- data.table(
  parameter = c("scoring_method", "assay", "slot", "nbin", "ctrl", "seed",
                "score_unit", "redundancy_stat", "statistical_unit",
                "proximity_run"),
  value = c("Seurat::AddModuleScore", "SCT", "data", "24", "100", "1",
            "per-spot continuous score", "per-slice Spearman rho",
            "slice (n=18)", "not run in this script")
)

fwrite(score_dt, file.path(out_dir, "R9_stage3_IVY_hypoxia_scores_per_spot.csv"))
fwrite(coverage, file.path(out_dir, "R9_stage3_IVY_hypoxia_score_gene_coverage.csv"))
fwrite(score_col_map, file.path(out_dir, "R9_stage3_IVY_hypoxia_score_column_map.csv"))
fwrite(params, file.path(out_dir, "R9_stage3_IVY_hypoxia_scoring_parameters.csv"))
fwrite(per_slice, file.path(out_dir, "R9_stage3_redundancy_correlations_per_slice.csv"))
fwrite(summary, file.path(out_dir, "R9_stage3_redundancy_correlations_summary.csv"))
fwrite(focus_summary, file.path(out_dir, "R9_stage3_redundancy_focus_pairs_summary.csv"))

cat("== Scoring parameters:\n")
print(params)

cat("\n== Gene coverage:\n")
print(coverage[, .(signature, n_input, n_used, frac_used)])

cat("\n== Focus redundancy pairs:\n")
print(focus_summary[order(pair_label)])

cat("\n== IVY/hypoxia/vascular redundancy matrix summary written. STOP before proximity.\n")
