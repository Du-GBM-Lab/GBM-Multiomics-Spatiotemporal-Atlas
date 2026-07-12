# =============================================================================
# R2 · scWGCNA (hdWGCNA) —— STAGE 7 focus lookup + tau specificity
# -----------------------------------------------------------------------------
# Input: STAGE 5-6 tables only.
# Output: MES-V module focus-gene hits, focus-gene assignment across all modules,
# MES-V hub summary, and module specificity tau.
#
# Boundary: this script does not save or modify Seurat/hdWGCNA objects.
# =============================================================================

options(stringsAsFactors = FALSE)

## ===== CONFIG ================================================================
OUT_DIR <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/R2_scWGCNA"
tab_dir <- file.path(OUT_DIR, "tables")

FOCUS_GENES <- c(
  "PLAUR", "FOSL1", "FOSL2", "FOS", "FOSB",
  "JUN", "JUNB", "JUND", "BATF", "ATF3"
)
MESV_MODULES <- c("greenyellow", "green", "salmon")
## ===========================================================================

cat("\n########## STAGE 7: MES-V focus-gene lookup ##########\n")

modules <- read.csv(
  file.path(tab_dir, "stage5_modules_gene_assignment_kME.csv"),
  check.names = FALSE
)
hub_df <- read.csv(
  file.path(tab_dir, "stage5_hub_genes_top25.csv"),
  check.names = FALSE
)
assoc <- read.csv(
  file.path(tab_dir, "module_subtype_association.csv"),
  check.names = FALSE
)
me <- read.csv(
  file.path(tab_dir, "ME_mean_by_subtype.csv"),
  row.names = 1,
  check.names = FALSE
)

stopifnot(all(c("gene_name", "module", "color") %in% colnames(modules)))
stopifnot(all(c("gene_name", "module", "kME") %in% colnames(hub_df)))
stopifnot(all(MESV_MODULES %in% modules$module))

all_hits <- list()

for (mm in MESV_MODULES) {
  k_col <- paste0("kME_", mm)
  if (!k_col %in% colnames(modules)) {
    stop("Missing kME column for module ", mm, ": ", k_col)
  }

  in_mod <- modules[modules$module == mm, , drop = FALSE]
  in_mod <- in_mod[order(-in_mod[[k_col]]), , drop = FALSE]
  in_mod$kME <- in_mod[[k_col]]
  in_mod$kME_rank <- seq_len(nrow(in_mod))

  top_hubs <- head(in_mod[, c("gene_name", "module", "kME", "kME_rank")], 25)
  write.csv(
    top_hubs,
    file.path(tab_dir, paste0("stage7_", mm, "_top_hubs_by_kME.csv")),
    row.names = FALSE
  )

  hit <- in_mod[in_mod$gene_name %in% FOCUS_GENES, c("gene_name", "module", "kME", "kME_rank"), drop = FALSE]
  hit$queried_module <- rep(mm, nrow(hit))
  all_hits[[mm]] <- hit

  cat(sprintf("\n-- MES-V module %s: %d genes\n", mm, nrow(in_mod)))
  cat("Top hubs:\n")
  print(top_hubs)
  cat("Focus-gene hits in this module:\n")
  print(hit)
}

mesv_focus_hits <- do.call(rbind, all_hits)
if (is.null(mesv_focus_hits)) {
  mesv_focus_hits <- data.frame(
    gene_name = character(), module = character(), kME = numeric(),
    kME_rank = integer(), queried_module = character()
  )
}
write.csv(
  mesv_focus_hits,
  file.path(tab_dir, "stage7_MESV_modules_focus_gene_hits.csv"),
  row.names = FALSE
)

focus_assign <- modules[modules$gene_name %in% FOCUS_GENES, , drop = FALSE]
focus_assign$own_module_kME <- NA_real_
focus_assign$own_module_kME_rank <- NA_integer_
for (i in seq_len(nrow(focus_assign))) {
  mm <- focus_assign$module[i]
  k_col <- paste0("kME_", mm)
  if (k_col %in% colnames(modules)) {
    in_mod <- modules[modules$module == mm, , drop = FALSE]
    in_mod <- in_mod[order(-in_mod[[k_col]]), , drop = FALSE]
    ranks <- setNames(seq_len(nrow(in_mod)), in_mod$gene_name)
    focus_assign$own_module_kME[i] <- focus_assign[[k_col]][i]
    focus_assign$own_module_kME_rank[i] <- unname(ranks[focus_assign$gene_name[i]])
  }
}
focus_assign <- focus_assign[match(intersect(FOCUS_GENES, focus_assign$gene_name), focus_assign$gene_name), ]
write.csv(
  focus_assign,
  file.path(tab_dir, "focus_genes_module_assignment.csv"),
  row.names = FALSE
)

cat("\nFocus genes across all modules:\n")
print(focus_assign[, c("gene_name", "module", "own_module_kME", "own_module_kME_rank")])

mesv_hubs <- hub_df[hub_df$module %in% MESV_MODULES, , drop = FALSE]
write.csv(
  mesv_hubs,
  file.path(tab_dir, "stage7_MESV_modules_hub_genes_top25.csv"),
  row.names = FALSE
)

cat("\nHub genes top25 for MES-V-associated modules:\n")
print(mesv_hubs)

cat("\n########## tau specificity from ME_mean_by_subtype ##########\n")
tau <- apply(me, 1, function(x) {
  x <- x - min(x)
  if (max(x) == 0) return(0)
  sum(1 - x / max(x)) / (length(x) - 1)
})
spec <- data.frame(
  module = rownames(me),
  top_subtype = colnames(me)[max.col(me)],
  tau = round(tau, 3),
  stringsAsFactors = FALSE
)
spec <- spec[order(-spec$tau), ]
write.csv(
  spec,
  file.path(tab_dir, "module_specificity_tau.csv"),
  row.names = FALSE
)
print(spec)

cat("\n########## STOP after STAGE 7 tables; no object saved ##########\n")
