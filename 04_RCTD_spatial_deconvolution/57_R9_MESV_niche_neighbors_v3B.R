#!/usr/bin/env Rscript
# =============================================================================
# R9 | MES-V niche neighbor profile under strict v3-B
#
# Nature:
#   Expand locked v3-B from "MES-V near vascular" to a conservative MES-V
#   neighborhood composition profile. This is spatial association only. It does
#   not infer recruitment, interaction, physical contact, single-cell
#   colocalization, or causality.
#
# Method:
#   - Focal response: MES-V top10 spots.
#   - Predictors: non-malignant class top10 spots.
#   - Neuron response control is run for every predictor.
#   - Random-labeling null fixes observed spot geometry and point counts.
#   - CSR and label permutation over expression/category labels are not used.
#   - Focal self-overlap is excluded from predictor-neighbor calculations.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(qs2)
  library(Rcpp)
})

base_dir <- getwd()
out_dir <- file.path(base_dir, "tables/R9_MESV_niche_neighbors")
doc_dir <- file.path(base_dir, "docs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(doc_dir, showWarnings = FALSE, recursive = TRUE)

weights_path <- file.path(base_dir, "tables/C3_4_local_niche/R9_A2_RCTD_weights_allslices_long.qs2")
threshold_path <- file.path(base_dir, "tables/C3_niche_preflight/R9_C3_1_neighborhood_niche_scores.csv")
set.seed(20260612 + 57)

n_perm <- 999L
k_values <- c(6L, 12L)
top_fraction <- 0.10

predictor_meta <- data.table(
  predictor = c("vascular", "myeloid_TAM", "Astrocyte", "Oligo_lineage", "T_cell", "Other_nonmalignant"),
  tier = c("primary_positive_anchor", "primary_TAM_null_control", rep("exploratory", 4)),
  interpretation_guard = c(
    "positive anchor expected to reproduce v3-B MES-V vascular proximity",
    "expected null; if significant, internal-audit until reviewed, never TAM recruitment",
    rep("exploratory; report only if neuron clean and 15-17/18 stable", 4)
  )
)

Rcpp::cppFunction('
  List run_obs_null_cpp(const NumericMatrix& D, const LogicalMatrix& Adj,
                        const IntegerVector& A_obs, const IntegerVector& B_obs,
                        int n, int nA, int nB, int n_perm) {
    RNGScope scope;
    auto sample_k_internal = [](int n, int k) {
      std::vector<int> vals(n);
      for (int i = 0; i < n; ++i) vals[i] = i;
      for (int i = 0; i < k; ++i) {
        int j = i + (int) floor(R::runif(0.0, 1.0) * (n - i));
        std::swap(vals[i], vals[j]);
      }
      IntegerVector out(k);
      for (int i = 0; i < k; ++i) out[i] = vals[i];
      return out;
    };
    auto calc_one_internal = [](const NumericMatrix& D, const LogicalMatrix& Adj,
                                const IntegerVector& A, const IntegerVector& B) {
      int nA = A.size();
      int nB = B.size();
      double cross_sum = 0.0;
      double nn_sum = 0.0;
      int nn_count = 0;
      for (int ai = 0; ai < nA; ++ai) {
        int a = A[ai];
        int row_count = 0;
        double min_d = R_PosInf;
        for (int bi = 0; bi < nB; ++bi) {
          int b = B[bi];
          if (a == b) continue;
          if (Adj(a, b)) row_count++;
          double d = D(a, b);
          if (d < min_d) min_d = d;
        }
        cross_sum += row_count;
        if (R_finite(min_d)) {
          nn_sum += min_d;
          nn_count++;
        }
      }
      NumericVector out(2);
      out[0] = nA > 0 ? cross_sum / nA : NA_REAL;
      out[1] = nn_count > 0 ? nn_sum / nn_count : NA_REAL;
      return out;
    };
    NumericVector obs = calc_one_internal(D, Adj, A_obs, B_obs);
    NumericMatrix out(n_perm, 2);
    for (int p = 0; p < n_perm; ++p) {
      IntegerVector A = sample_k_internal(n, nA);
      IntegerVector B = sample_k_internal(n, nB);
      NumericVector st = calc_one_internal(D, Adj, A, B);
      out(p, 0) = st[0];
      out(p, 1) = st[1];
    }
    return List::create(Named("obs") = obs, Named("null") = out);
  }
')

make_hot <- function(x, frac) {
  x >= as.numeric(quantile(x, 1 - frac, na.rm = TRUE))
}

emp_p <- function(obs, nul, side = c("greater", "less")) {
  side <- match.arg(side)
  if (side == "greater") return((1 + sum(nul >= obs, na.rm = TRUE)) / (sum(is.finite(nul)) + 1))
  (1 + sum(nul <= obs, na.rm = TRUE)) / (sum(is.finite(nul)) + 1)
}

calc_stats <- function(D, Adj, A_idx, B_idx) {
  if (!length(A_idx) || !length(B_idx)) {
    return(c(crossK_count = NA_real_, mean_nearest_distance = NA_real_))
  }
  sub_adj <- Adj[A_idx, B_idx, drop = FALSE]
  sub_d <- D[A_idx, B_idx, drop = FALSE]
  overlap <- match(A_idx, B_idx, nomatch = 0L)
  if (any(overlap > 0L)) {
    rr <- which(overlap > 0L)
    cc <- overlap[overlap > 0L]
    sub_adj[cbind(rr, cc)] <- FALSE
    sub_d[cbind(rr, cc)] <- Inf
  }
  nn <- apply(sub_d, 1, min)
  c(
    crossK_count = mean(rowSums(sub_adj)),
    mean_nearest_distance = mean(nn[is.finite(nn)])
  )
}

run_one <- function(dt, D, radius, response, predictor) {
  Adj <- D <= radius & D > 0
  A <- make_hot(dt[[response]], top_fraction)
  B <- make_hot(dt[[predictor]], top_fraction)
  n <- nrow(dt)
  nA <- sum(A)
  nB <- sum(B)
  sim <- run_obs_null_cpp(D, Adj, which(A) - 1L, which(B) - 1L, n, nA, nB, n_perm)
  obs <- sim$obs
  null <- sim$null
  nullK <- null[, 1]
  nullNN <- null[, 2]

  data.table(
    response = response,
    predictor = predictor,
    top_rule = "top10",
    top_fraction = top_fraction,
    n_spots = n,
    n_response_hot = nA,
    n_predictor_hot = nB,
    radius = radius,
    count_preserved_all_perm = TRUE,
    observed_crossK_count = obs[[1]],
    null_crossK_median = median(nullK, na.rm = TRUE),
    delta_crossK = obs[[1]] - median(nullK, na.rm = TRUE),
    p_emp_crossK_greater = emp_p(obs[[1]], nullK, "greater"),
    observed_mean_nn_dist = obs[[2]],
    null_mean_nn_dist_median = median(nullNN, na.rm = TRUE),
    delta_nn_closer = median(nullNN, na.rm = TRUE) - obs[[2]],
    p_emp_nn_less = emp_p(obs[[2]], nullNN, "less"),
    n_perm = n_perm
  )
}

sign_p <- function(x, alternative = "greater") {
  x <- x[is.finite(x) & x != 0]
  if (!length(x)) return(NA_real_)
  binom.test(sum(x > 0), length(x), p = 0.5, alternative = alternative)$p.value
}

wilcox_p <- function(x, alternative = "greater") {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  suppressWarnings(wilcox.test(x, mu = 0, alternative = alternative, exact = FALSE)$p.value)
}

weights <- as.data.table(qs2::qs_read(weights_path))
thresholds <- fread(threshold_path)[, .(distance_threshold = median(distance_threshold, na.rm = TRUE)), by = .(slice, k)]

required <- c(
  "spot_id", "slice", "x", "y", "Subtype3",
  "Endothelial", "Mural cells", "Macrophages", "Microglial", "Monocytes",
  "Astrocytes", "Oligodendrocytes", "OPCs", "T cells", "B cells", "NK cells",
  "cDCs", "pDCs", "Ependymal cells", "Radial glial", "Ambiguous", "Neurons"
)
missing <- setdiff(required, names(weights))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

feat <- data.table(
  spot_id = weights$spot_id,
  slice = weights$slice,
  x = weights$x,
  y = weights$y,
  `MES-V` = weights$Subtype3,
  neuron_control = weights$Neurons,
  vascular = weights$Endothelial + weights[["Mural cells"]],
  myeloid_TAM = weights$Macrophages + weights$Microglial + weights$Monocytes,
  Astrocyte = weights$Astrocytes,
  Oligo_lineage = weights$Oligodendrocytes + weights$OPCs,
  T_cell = weights[["T cells"]],
  Other_nonmalignant = weights[["B cells"]] + weights[["NK cells"]] + weights$cDCs + weights$pDCs +
    weights[["Ependymal cells"]] + weights[["Radial glial"]] + weights$Ambiguous
)

responses <- c("MES-V", "neuron_control")
predictors <- predictor_meta$predictor

params <- data.table(
  parameter = c("input_weights", "distance_threshold_source", "responses", "predictors", "top_fraction", "k_values", "n_perm", "null", "focal_overlap_handling", "TAM_guard"),
  value = c(
    weights_path,
    threshold_path,
    paste(responses, collapse = ", "),
    paste(predictors, collapse = ", "),
    as.character(top_fraction),
    paste(k_values, collapse = ", "),
    as.character(n_perm),
    "v3-B-style random labeling on fixed observed spot geometry with response and predictor point counts preserved; no CSR",
    "same focal spot is excluded when response and predictor hot sets overlap",
    "myeloid_TAM is predeclared primary null-control; significant result is internal-audit until reviewed, never recruitment"
  )
)
fwrite(params, file.path(out_dir, "MESV_niche_neighbors_parameters_source.csv"))
fwrite(predictor_meta, file.path(out_dir, "MESV_niche_neighbors_predictor_tiers.csv"))

res <- list()
i <- 0L
for (sl in sort(unique(feat$slice))) {
  dt <- feat[slice == sl]
  D <- as.matrix(dist(as.matrix(dt[, .(x, y)])))
  diag(D) <- 0
  for (kk in k_values) {
    radius <- thresholds[slice == sl & k == kk, distance_threshold][1]
    if (!length(radius) || is.na(radius)) stop("Missing radius for slice=", sl, " k=", kk)
    for (resp in responses) {
      for (pred in predictors) {
        i <- i + 1L
        ans <- run_one(dt, D, radius, resp, pred)
        ans[, `:=`(slice = sl, k = kk)]
        res[[i]] <- ans
      }
    }
  }
  message("done ", sl)
}

perslice <- rbindlist(res, use.names = TRUE)
perslice <- merge(perslice, predictor_meta, by = "predictor", all.x = TRUE)
setcolorder(perslice, c("slice", "k", "response", "predictor", "tier"))
perslice[, FDR_crossK_by_slice_response_k := p.adjust(p_emp_crossK_greater, method = "BH"), by = .(slice, response, k)]
perslice[, FDR_nn_by_slice_response_k := p.adjust(p_emp_nn_less, method = "BH"), by = .(slice, response, k)]
perslice[, crossK_pass_FDR := delta_crossK > 0 & FDR_crossK_by_slice_response_k < 0.05]
perslice[, nn_pass_FDR := delta_nn_closer > 0 & FDR_nn_by_slice_response_k < 0.05]

summary <- perslice[, .(
  n_slices = .N,
  median_delta_crossK = median(delta_crossK, na.rm = TRUE),
  q25_delta_crossK = as.numeric(quantile(delta_crossK, 0.25, na.rm = TRUE)),
  q75_delta_crossK = as.numeric(quantile(delta_crossK, 0.75, na.rm = TRUE)),
  n_crossK_positive = sum(delta_crossK > 0, na.rm = TRUE),
  n_slice_crossK_emp_p_lt_0p05 = sum(delta_crossK > 0 & p_emp_crossK_greater < 0.05, na.rm = TRUE),
  n_slice_crossK_FDR_lt_0p05 = sum(crossK_pass_FDR, na.rm = TRUE),
  sign_p_crossK_greater = sign_p(delta_crossK),
  wilcox_p_crossK_greater = wilcox_p(delta_crossK),
  median_delta_nn_closer = median(delta_nn_closer, na.rm = TRUE),
  q25_delta_nn_closer = as.numeric(quantile(delta_nn_closer, 0.25, na.rm = TRUE)),
  q75_delta_nn_closer = as.numeric(quantile(delta_nn_closer, 0.75, na.rm = TRUE)),
  n_nn_closer_positive = sum(delta_nn_closer > 0, na.rm = TRUE),
  n_slice_nn_emp_p_lt_0p05 = sum(delta_nn_closer > 0 & p_emp_nn_less < 0.05, na.rm = TRUE),
  n_slice_nn_FDR_lt_0p05 = sum(nn_pass_FDR, na.rm = TRUE),
  sign_p_nn_greater = sign_p(delta_nn_closer),
  wilcox_p_nn_greater = wilcox_p(delta_nn_closer),
  count_preserved_all = all(count_preserved_all_perm),
  median_response_hot = median(n_response_hot, na.rm = TRUE),
  median_predictor_hot = median(n_predictor_hot, na.rm = TRUE)
), by = .(k, response, predictor, tier, top_rule, top_fraction)]

mesv <- summary[response == "MES-V"]
neuron <- summary[response == "neuron_control", .(
  k,
  predictor,
  neuron_crossK_FDR_pass_slices = n_slice_crossK_FDR_lt_0p05,
  neuron_nn_FDR_pass_slices = n_slice_nn_FDR_lt_0p05,
  neuron_crossK_positive_slices = n_crossK_positive,
  neuron_nn_positive_slices = n_nn_closer_positive
)]
gate <- merge(mesv, neuron, by = c("k", "predictor"), all.x = TRUE)
gate[, neuron_clean := neuron_crossK_FDR_pass_slices == 0 & neuron_nn_FDR_pass_slices == 0]
gate[, mesv_crossK_gate_15 := n_slice_crossK_FDR_lt_0p05 >= 15 & n_crossK_positive >= 15]
gate[, mesv_nn_gate_15 := n_slice_nn_FDR_lt_0p05 >= 15 & n_nn_closer_positive >= 15]
gate[, overall_v3B_neighbor_gate := neuron_clean & mesv_crossK_gate_15 & mesv_nn_gate_15 & count_preserved_all]
gate[, TAM_hard_gate_status := fifelse(
  predictor == "myeloid_TAM" & overall_v3B_neighbor_gate,
  "unexpected_significant_internal_audit_do_not_revive_TAM",
  fifelse(predictor == "myeloid_TAM", "null_or_not_gate_consistent_with_TAM_lock", "not_TAM")
)]

fwrite(perslice, file.path(out_dir, "MESV_niche_neighbors_v3B_perslice.csv"))
fwrite(summary, file.path(out_dir, "MESV_niche_neighbors_v3B_summary.csv"))
fwrite(gate, file.path(out_dir, "MESV_niche_neighbors_v3B_gate_summary.csv"))
qs2::qs_save(
  list(perslice = perslice, summary = summary, gate = gate, parameters = params, predictor_meta = predictor_meta),
  file.path(out_dir, "MESV_niche_neighbors_v3B_results.qs2")
)

fmt_gate <- function(pred, kk = 6L) {
  g <- gate[k == kk & predictor == pred]
  if (!nrow(g)) return("NA")
  paste0(
    pred, ": MES-V crossK ", g$n_slice_crossK_FDR_lt_0p05, "/18 FDR, ",
    "NN ", g$n_slice_nn_FDR_lt_0p05, "/18 FDR; ",
    "neuron crossK ", g$neuron_crossK_FDR_pass_slices, "/18, NN ", g$neuron_nn_FDR_pass_slices, "/18; ",
    "overall=", g$overall_v3B_neighbor_gate
  )
}

vascular_anchor <- gate[k == 6 & predictor == "vascular", overall_v3B_neighbor_gate][1]
myeloid_gate <- gate[k == 6 & predictor == "myeloid_TAM", overall_v3B_neighbor_gate][1]
nonvascular_pass <- gate[k == 6 & predictor != "vascular" & predictor != "myeloid_TAM" & overall_v3B_neighbor_gate == TRUE, predictor]
placement <- if (isTRUE(vascular_anchor) && length(nonvascular_pass)) {
  "supplement_with_possible_main_text_sentence_for_clean_nonvascular_neighbors"
} else if (isTRUE(vascular_anchor)) {
  "supplement_profile_vascular_anchor_reproduced_no_additional_clean_neighbors"
} else {
  "internal_audit_anchor_not_reproduced"
}

stop_lines <- c(
  "# R9 MES-V niche neighbor profile v3-B STOP",
  "",
  "性质声明：本分析以 MES-V 为主语，按 v3-B-style random-labeling null 描述其非恶性空间邻居画像；不证明 physical contact、single-cell colocalization、interaction、recruitment 或 causality。",
  "",
  "## Method",
  paste0("- responses：", paste(responses, collapse = ", "), "；predictors：", paste(predictors, collapse = ", "), "。"),
  paste0("- top rule：top10；k=", paste(k_values, collapse = "/"), "；n_perm=", n_perm, "；focal overlap excluded。"),
  "- FDR：within each slice/response/k across all non-malignant predictor classes, separately for cross-K and nearest-distance.",
  paste0("- count preservation all rows：", all(perslice$count_preserved_all_perm), " (", sum(perslice$count_preserved_all_perm), "/", nrow(perslice), ")。"),
  "",
  "## k=6 Gate Summary",
  paste0("- ", fmt_gate("vascular", 6L), "。"),
  paste0("- ", fmt_gate("myeloid_TAM", 6L), "。"),
  paste0("- ", fmt_gate("Astrocyte", 6L), "。"),
  paste0("- ", fmt_gate("Oligo_lineage", 6L), "。"),
  paste0("- ", fmt_gate("T_cell", 6L), "。"),
  paste0("- ", fmt_gate("Other_nonmalignant", 6L), "。"),
  "",
  "## k=12 Sensitivity Summary",
  paste0("- ", fmt_gate("vascular", 12L), "。"),
  paste0("- ", fmt_gate("myeloid_TAM", 12L), "。"),
  paste0("- ", fmt_gate("Astrocyte", 12L), "。"),
  paste0("- ", fmt_gate("Oligo_lineage", 12L), "。"),
  paste0("- ", fmt_gate("T_cell", 12L), "。"),
  paste0("- ", fmt_gate("Other_nonmalignant", 12L), "。"),
  "",
  "## TAM hard gate",
  if (isTRUE(myeloid_gate)) {
    "- myeloid_TAM unexpectedly meets the gate at k=6; per predeclaration this is internal-audit pending review and must not be written as TAM recruitment."
  } else {
    "- myeloid_TAM does not meet the strict MES-V neighbor gate at k=6; this remains consistent with the locked TAM spatial null."
  },
  "",
  "## Placement",
  paste0("- self-assessed placement：", placement, "。"),
  "- Vascular is the positive anchor. Any non-vascular class requires neuron clean + cross-K and NN 15-17/18 equivalent to be mentioned beyond supplement audit.",
  "",
  "## Forbidden wording",
  "- 禁写 recruitment、interaction、physical contact、single-cell colocalization、drives、causality、TAM recruitment。",
  "",
  "## Source Data",
  paste0("- input weights：`", weights_path, "`"),
  paste0("- output dir：`", out_dir, "`")
)
writeLines(stop_lines, file.path(doc_dir, "R9_MESV_niche_neighbors_v3B_STOP.md"), useBytes = TRUE)

message("Wrote MES-V niche neighbor outputs to: ", out_dir)
message("Wrote STOP to: ", file.path(doc_dir, "R9_MESV_niche_neighbors_v3B_STOP.md"))
