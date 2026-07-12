# R4 05: STOP2 audit for subtype gene-set coverage and ssGSEA scores.
# Read hgg_data.rds + sc_subtype_markers.rds; no clustering, KM, or Cox.

suppressPackageStartupMessages({
  library(dplyr)
})

root <- getwd()
raw_dir <- file.path(root, "data", "raw", "CGGA_old_manuscript")
proc_dir <- file.path(root, "data", "processed")
tab_dir <- file.path(root, "tables")
fig_dir <- file.path(root, "figures")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

FORBIDDEN <- c("PLAUR", "FOSL1", "FOSL2")
SUBTYPE_LEVELS <- c("NPC_P", "OPC_M", "MES_V", "MES_I")

hgg_path <- file.path(raw_dir, "hgg_data.rds")
gs_path <- file.path(proc_dir, "sc_subtype_markers.rds")
cat("hgg_path:", hgg_path, "\n")
cat("gene_set_path:", gs_path, "\n")
stopifnot(file.exists(hgg_path), file.exists(gs_path))

hgg <- readRDS(hgg_path)
clin <- hgg$clinical
expr_df <- hgg$expression
gene_sets <- readRDS(gs_path)
stopifnot(all(SUBTYPE_LEVELS %in% names(gene_sets)))

forbidden_hits <- intersect(toupper(unlist(gene_sets)), FORBIDDEN)
stopifnot(length(forbidden_hits) == 0)

stopifnot("Gene_Name" %in% colnames(expr_df))
genes <- as.character(expr_df$Gene_Name)
expr_mat <- as.matrix(expr_df[, setdiff(colnames(expr_df), "Gene_Name")])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- genes

na_n <- sum(is.na(expr_mat))
if (na_n > 0) {
  cat("NA expression values detected:", na_n, "; imputing row medians for audit.\n")
  row_med <- apply(expr_mat, 1, median, na.rm = TRUE)
  idx <- which(is.na(expr_mat), arr.ind = TRUE)
  expr_mat[idx] <- row_med[idx[, 1]]
}

dup_genes <- sum(duplicated(rownames(expr_mat)))
if (dup_genes > 0) {
  cat("duplicated gene symbols:", dup_genes, "; aggregating by row mean.\n")
  sums <- rowsum(expr_mat, group = rownames(expr_mat), reorder = FALSE)
  counts <- as.vector(table(factor(rownames(expr_mat), levels = rownames(sums))))
  expr_mat <- sweep(sums, 1, counts, "/")
}

expr_samples <- colnames(expr_mat)
clin_ids <- as.character(clin$CGGA_ID)
common <- intersect(expr_samples, clin_ids)
cat("expression samples:", length(expr_samples), "\n")
cat("clinical samples:", nrow(clin), "\n")
cat("overlap:", length(common), "\n")
stopifnot(length(common) == length(expr_samples), length(common) == nrow(clin))
clin <- clin[match(expr_samples, clin_ids), , drop = FALSE]

vals_peek <- as.numeric(expr_mat[seq_len(min(100, nrow(expr_mat))), seq_len(min(100, ncol(expr_mat)))])
integer_frac <- mean(abs(vals_peek - round(vals_peek)) < 1e-8, na.rm = TRUE)
scale_guess <- if (max(vals_peek, na.rm = TRUE) > 1000 && integer_frac > 0.9) {
  "raw-count-like"
} else if (max(vals_peek, na.rm = TRUE) > 100) {
  "unlogged abundance (FPKM/TPM/RSEM-like)"
} else {
  "log-like or low-range abundance"
}
scale_tbl <- data.frame(
  metric = c("min", "q25", "median", "mean", "q75", "max", "integer_fraction", "scale_guess"),
  value = c(
    min(vals_peek, na.rm = TRUE),
    quantile(vals_peek, 0.25, na.rm = TRUE),
    median(vals_peek, na.rm = TRUE),
    mean(vals_peek, na.rm = TRUE),
    quantile(vals_peek, 0.75, na.rm = TRUE),
    max(vals_peek, na.rm = TRUE),
    round(integer_frac, 3),
    scale_guess
  )
)
write.csv(scale_tbl, file.path(tab_dir, "STOP2_表达量纲审计.csv"), row.names = FALSE)

available <- rownames(expr_mat)
coverage <- do.call(rbind, lapply(names(gene_sets), function(nm) {
  gs <- unique(as.character(gene_sets[[nm]]))
  hit <- intersect(gs, available)
  data.frame(
    subtype = nm,
    n_input = length(gs),
    n_available = length(hit),
    coverage = round(length(hit) / length(gs), 3),
    missing_genes = paste(setdiff(gs, available), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
write.csv(coverage, file.path(tab_dir, "STOP2_gene_set_CGGA覆盖率.csv"), row.names = FALSE)

overlap_pairs <- data.frame()
for (i in seq_along(gene_sets)) {
  for (j in seq_along(gene_sets)) {
    a <- names(gene_sets)[i]
    b <- names(gene_sets)[j]
    inter <- intersect(gene_sets[[a]], gene_sets[[b]])
    union <- union(gene_sets[[a]], gene_sets[[b]])
    overlap_pairs <- rbind(overlap_pairs, data.frame(
      set1 = a,
      set2 = b,
      n_intersect = length(inter),
      jaccard = ifelse(length(union) == 0, NA, round(length(inter) / length(union), 3)),
      stringsAsFactors = FALSE
    ))
  }
}
write.csv(overlap_pairs, file.path(tab_dir, "STOP2_gene_set_overlap_jaccard.csv"), row.names = FALSE)

gene_sets_avail <- lapply(gene_sets, function(gs) intersect(gs, available))
cat("\n== coverage ==\n")
print(coverage[, c("subtype", "n_input", "n_available", "coverage")])
cat("\n== gene set overlap ==\n")
print(overlap_pairs)

cat("\n== ssGSEA package audit ==\n")
cat("GSVA installed:", requireNamespace("GSVA", quietly = TRUE), "\n")
if (!requireNamespace("GSVA", quietly = TRUE)) {
  stop("GSVA is not installed. Install Bioconductor GSVA before STOP2 ssGSEA.")
}
cat("GSVA version:", as.character(utils::packageVersion("GSVA")), "\n")

ssgsea_scores <- NULL
if (exists("ssgseaParam", where = asNamespace("GSVA"), mode = "function")) {
  cat("Using GSVA parameter-object API: gsva(ssgseaParam(...)).\n")
  param <- GSVA::ssgseaParam(expr_mat, gene_sets_avail, normalize = TRUE)
  ssgsea_scores <- GSVA::gsva(param, verbose = FALSE)
} else {
  cat("Using legacy GSVA API fallback.\n")
  ssgsea_scores <- GSVA::gsva(expr_mat, gene_sets_avail, method = "ssgsea", verbose = FALSE)
}
ssgsea_scores <- as.matrix(ssgsea_scores)
ssgsea_scores <- ssgsea_scores[SUBTYPE_LEVELS, , drop = FALSE]
saveRDS(ssgsea_scores, file.path(proc_dir, "ssgsea_scores.rds"))
write.csv(data.frame(subtype = rownames(ssgsea_scores), ssgsea_scores, check.names = FALSE),
          file.path(tab_dir, "STOP2_ssGSEA_scores.csv"), row.names = FALSE)

score_long <- data.frame(
  sample = rep(colnames(ssgsea_scores), each = nrow(ssgsea_scores)),
  subtype = rep(rownames(ssgsea_scores), times = ncol(ssgsea_scores)),
  score = as.vector(ssgsea_scores),
  stringsAsFactors = FALSE
)
score_summary <- score_long %>%
  group_by(subtype) %>%
  summarise(
    min = min(score, na.rm = TRUE),
    q25 = quantile(score, 0.25, na.rm = TRUE),
    median = median(score, na.rm = TRUE),
    mean = mean(score, na.rm = TRUE),
    q75 = quantile(score, 0.75, na.rm = TRUE),
    max = max(score, na.rm = TRUE),
    sd = sd(score, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(score_summary, file.path(tab_dir, "STOP2_ssGSEA_score分布.csv"), row.names = FALSE)

score_cor <- cor(t(ssgsea_scores), method = "spearman", use = "pairwise.complete.obs")
write.csv(score_cor, file.path(tab_dir, "STOP2_ssGSEA_score_Spearman相关.csv"))

score_df <- data.frame(CGGA_ID = colnames(ssgsea_scores), t(ssgsea_scores), check.names = FALSE)
ann <- c("PRS_type", "Grade", "IDH_mutation_status", "X1p19q_codeletion_status",
         "MGMTp_methylation_status", "Histology")
clin_score <- cbind(clin[, intersect(colnames(clin), ann), drop = FALSE], score_df)
sanity <- data.frame()
for (grp in intersect(colnames(clin), ann)) {
  for (st in SUBTYPE_LEVELS) {
    tmp <- data.frame(group = as.character(clin[[grp]]), score = as.numeric(ssgsea_scores[st, ]))
    tmp <- tmp[!is.na(tmp$group) & tmp$group != "" & !is.na(tmp$score), ]
    if (!nrow(tmp)) next
    gsum <- tmp %>%
      group_by(group) %>%
      summarise(n = n(), mean_score = mean(score), median_score = median(score), .groups = "drop")
    gsum$variable <- grp
    gsum$subtype <- st
    sanity <- rbind(sanity, gsum[, c("variable", "subtype", "group", "n", "mean_score", "median_score")])
  }
}
write.csv(sanity, file.path(tab_dir, "STOP2_ssGSEA_score_by_clinical_annotation.csv"), row.names = FALSE)

pdf(file.path(fig_dir, "STOP2_ssGSEA_score_distribution.pdf"), width = 7, height = 4)
boxplot(score ~ subtype, data = score_long, las = 1, col = "#D9E6F2",
        ylab = "ssGSEA score", xlab = "", main = "Subtype ssGSEA scores in CGGA HGG")
dev.off()

pdf(file.path(fig_dir, "STOP2_ssGSEA_score_correlation.pdf"), width = 5, height = 5)
if (requireNamespace("pheatmap", quietly = TRUE)) {
  pheatmap::pheatmap(score_cor, cluster_rows = FALSE, cluster_cols = FALSE,
                     main = "Spearman correlation of subtype scores")
} else {
  heatmap(score_cor, symm = TRUE, main = "Spearman correlation of subtype scores")
}
dev.off()

cat("\n== score distribution ==\n")
print(score_summary)
cat("\n== score Spearman correlation ==\n")
print(round(score_cor, 3))
cat("\n== selected sanity means ==\n")
print(sanity %>% filter(variable %in% c("IDH_mutation_status", "Grade", "PRS_type")) %>%
        arrange(variable, subtype, desc(mean_score)))

cat("\n########## STOP2 REPORT ##########\n")
cat("1. Gene-set coverage in CGGA: see STOP2_gene_set_CGGA覆盖率.csv and printed table.\n")
cat("2. MES_V/MES_I overlap/correlation: overlap table above; Spearman MES_V-MES_I = ",
    round(score_cor["MES_V", "MES_I"], 3), ".\n", sep = "")
cat("3. Score distribution: see STOP2_ssGSEA_score分布.csv and printed table.\n")
cat("4. Biological sanity: see STOP2_ssGSEA_score_by_clinical_annotation.csv; inspect IDH/grade/PRS direction before clustering.\n")
cat("5. Expression scale guess: ", scale_guess, "; see STOP2_表达量纲审计.csv.\n", sep = "")
cat("=> STOP2 complete. Do not proceed to clustering until reviewed.\n")
