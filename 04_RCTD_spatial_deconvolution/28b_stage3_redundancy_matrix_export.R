#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

base <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/R9_空间转录组/tables/R9_stage3_signature_scores"
sm <- fread(file.path(base, "R9_stage3_redundancy_correlations_summary.csv"))
vars <- c(
  "IVY_CT", "IVY_CTmvp", "IVY_CTpan", "IVY_LE",
  "Hypoxia_Buffa", "Hypoxia_Hallmark",
  "vascular", "myeloid", "microglia", "macrophage_monocyte", "neuron_control"
)
mat <- matrix(NA_real_, length(vars), length(vars), dimnames = list(vars, vars))
for (i in seq_len(nrow(sm))) {
  mat[sm$var1[i], sm$var2[i]] <- sm$median_rho[i]
  mat[sm$var2[i], sm$var1[i]] <- sm$median_rho[i]
}
diag(mat) <- 1
fwrite(as.data.table(mat, keep.rownames = "variable"),
       file.path(base, "R9_stage3_redundancy_median_rho_matrix.csv"))
print(round(mat, 3))
