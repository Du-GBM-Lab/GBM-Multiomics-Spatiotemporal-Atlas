# R4 06 / STOP2b: HGNC symbol harmonization before clustering.
# Goal: distinguish true low coverage from old/new gene-symbol mismatch.
# No gene-set curation, no clustering, no KM, no Cox.

suppressPackageStartupMessages({
  library(GSVA)
  library(HGNChelper)
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
stopifnot(file.exists(hgg_path), file.exists(gs_path))

hgg <- readRDS(hgg_path)
stopifnot("expression" %in% names(hgg), "Gene_Name" %in% colnames(hgg$expression))
expr_df <- hgg$expression
expr <- as.matrix(expr_df[, setdiff(colnames(expr_df), "Gene_Name")])
mode(expr) <- "numeric"
rownames(expr) <- as.character(expr_df$Gene_Name)

gene_sets <- readRDS(gs_path)
gene_sets <- gene_sets[SUBTYPE_LEVELS]
stopifnot(all(SUBTYPE_LEVELS %in% names(gene_sets)))
stopifnot(length(intersect(toupper(unlist(gene_sets)), FORBIDDEN)) == 0)

# Collapse duplicated expression row names before ssGSEA.
dup_genes <- sum(duplicated(rownames(expr)))
if (dup_genes > 0) {
  cat("duplicated expression symbols:", dup_genes, "; aggregating by mean.\n")
  sums <- rowsum(expr, group = rownames(expr), reorder = FALSE)
  counts <- as.vector(table(factor(rownames(expr), levels = rownames(sums))))
  expr <- sweep(sums, 1, counts, "/")
}

split_suggested <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  sub(" ///.*$", "", x)
}

to_current <- function(x) {
  x <- as.character(x)
  res <- suppressWarnings(HGNChelper::checkGeneSymbols(x, species = "human"))
  suggested <- split_suggested(res$Suggested.Symbol)
  # HGNChelper returns the original approved symbol in Suggested.Symbol for approved entries.
  out <- suggested
  out[is.na(out) | out == ""] <- NA_character_
  out
}

expr_cur <- to_current(rownames(expr))
expr_map <- data.frame(
  expr_symbol = rownames(expr),
  current_symbol = expr_cur,
  stringsAsFactors = FALSE
)
write.csv(expr_map, file.path(tab_dir, "STOP2b_CGGA表达基因_HGNC映射.csv"), row.names = FALSE)

matched_sets <- list()
harmonized_gene_sets <- list()
coverage_rows <- data.frame()
rescued_rows <- data.frame()

for (nm in names(gene_sets)) {
  gs <- unique(as.character(gene_sets[[nm]]))
  gs_cur <- to_current(gs)

  raw_hit <- intersect(gs, rownames(expr))
  cur_hit_symbols <- unique(gs_cur[!is.na(gs_cur) & gs_cur %in% expr_cur])
  expr_rows_for_set <- unique(rownames(expr)[!is.na(expr_cur) & expr_cur %in% cur_hit_symbols])

  matched_sets[[nm]] <- expr_rows_for_set
  harmonized_gene_sets[[nm]] <- cur_hit_symbols

  gs_map <- data.frame(
    subtype = nm,
    input_symbol = gs,
    current_symbol = gs_cur,
    raw_hit = gs %in% rownames(expr),
    harmonized_hit = !is.na(gs_cur) & gs_cur %in% expr_cur,
    stringsAsFactors = FALSE
  )
  rescued <- gs_map[!gs_map$raw_hit & gs_map$harmonized_hit, , drop = FALSE]
  rescued_rows <- rbind(rescued_rows, rescued)

  coverage_rows <- rbind(coverage_rows, data.frame(
    subtype = nm,
    n_input = length(gs),
    raw_n = length(raw_hit),
    raw_coverage = round(length(raw_hit) / length(gs), 3),
    harmonized_n = sum(gs_map$harmonized_hit),
    harmonized_coverage = round(sum(gs_map$harmonized_hit) / length(gs), 3),
    n_rescued = nrow(rescued),
    stringsAsFactors = FALSE
  ))
  write.csv(gs_map, file.path(tab_dir, paste0("STOP2b_", nm, "_gene_symbol_mapping.csv")), row.names = FALSE)
}

forbidden_after <- intersect(toupper(unlist(matched_sets)), FORBIDDEN)
stopifnot(length(forbidden_after) == 0)

write.csv(coverage_rows, file.path(tab_dir, "STOP2b_gene_set_覆盖率_符号归一前后.csv"), row.names = FALSE)
write.csv(rescued_rows, file.path(tab_dir, "STOP2b_HGNC符号归一捞回基因.csv"), row.names = FALSE)
write.csv(data.frame(subtype = rep(names(matched_sets), lengths(matched_sets)),
                     expr_symbol = unlist(matched_sets, use.names = FALSE)),
          file.path(tab_dir, "STOP2b_归一后实际入算基因.csv"), row.names = FALSE)

cat("== Coverage: raw intersect vs HGNC-harmonized ==\n")
print(coverage_rows)

cat("\n== Rescued genes by subtype ==\n")
for (nm in names(matched_sets)) {
  r <- rescued_rows[rescued_rows$subtype == nm, , drop = FALSE]
  cat("\n", nm, " rescued n=", nrow(r), "\n", sep = "")
  if (nrow(r)) print(r[, c("input_symbol", "current_symbol")])
}

cat("\n== Harmonized MES_V genes actually used in ssGSEA ==\n")
print(sort(matched_sets$MES_V))

cat("\n== ssGSEA package audit ==\n")
cat("GSVA:", as.character(utils::packageVersion("GSVA")), "\n")
cat("HGNChelper:", as.character(utils::packageVersion("HGNChelper")), "\n")

if (exists("ssgseaParam", where = asNamespace("GSVA"), mode = "function")) {
  param <- GSVA::ssgseaParam(expr, matched_sets, normalize = TRUE)
  es <- GSVA::gsva(param, verbose = FALSE)
} else {
  es <- GSVA::gsva(expr, matched_sets, method = "ssgsea", verbose = FALSE)
}
es <- as.matrix(es)
es <- es[SUBTYPE_LEVELS, , drop = FALSE]
saveRDS(es, file.path(proc_dir, "ssgsea_scores_harmonized.rds"))
write.csv(data.frame(subtype = rownames(es), es, check.names = FALSE),
          file.path(tab_dir, "STOP2b_ssGSEA_scores_harmonized.csv"), row.names = FALSE)

score_cor <- cor(t(es), method = "spearman", use = "pairwise.complete.obs")
write.csv(score_cor, file.path(tab_dir, "STOP2b_ssGSEA_score_harmonized_Spearman相关.csv"))

score_long <- data.frame(
  sample = rep(colnames(es), each = nrow(es)),
  subtype = rep(rownames(es), times = ncol(es)),
  score = as.vector(es),
  stringsAsFactors = FALSE
)
score_summary <- do.call(rbind, lapply(split(score_long$score, score_long$subtype), function(x) {
  data.frame(
    min = min(x, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE)
  )
}))
score_summary <- data.frame(subtype = rownames(score_summary), score_summary, row.names = NULL)
write.csv(score_summary, file.path(tab_dir, "STOP2b_ssGSEA_score_harmonized分布.csv"), row.names = FALSE)

pdf(file.path(fig_dir, "STOP2b_ssGSEA_harmonized_score_correlation.pdf"), width = 5, height = 5)
if (requireNamespace("pheatmap", quietly = TRUE)) {
  pheatmap::pheatmap(score_cor, cluster_rows = FALSE, cluster_cols = FALSE,
                     main = "Spearman correlation of harmonized subtype scores")
} else {
  heatmap(score_cor, symm = TRUE, main = "Spearman correlation of harmonized subtype scores")
}
dev.off()

pdf(file.path(fig_dir, "STOP2b_ssGSEA_harmonized_score_distribution.pdf"), width = 7, height = 4)
boxplot(score ~ subtype, data = score_long, las = 1, col = "#E6DCCB",
        ylab = "ssGSEA score", xlab = "", main = "HGNC-harmonized subtype ssGSEA scores")
dev.off()

cat("\n== Harmonized ssGSEA distribution ==\n")
print(score_summary)
cat("\n== Harmonized score Spearman ==\n")
print(round(score_cor, 3))

cat("\n\n########## STOP2b REPORT ##########\n")
cat("1. Harmonized coverage: see coverage table above; MES_V = ",
    coverage_rows$harmonized_n[coverage_rows$subtype == "MES_V"], "/",
    coverage_rows$n_input[coverage_rows$subtype == "MES_V"], " (",
    round(100 * coverage_rows$harmonized_coverage[coverage_rows$subtype == "MES_V"]), "%).\n", sep = "")
cat("2. Harmonized MES_V genes used are printed above and saved to STOP2b_归一后实际入算基因.csv.\n")
cat("3. Harmonized MES_V vs MES_I Spearman = ",
    round(score_cor["MES_V", "MES_I"], 3), ".\n", sep = "")
cat("=> STOP2b complete. Do not proceed to clustering until reviewed.\n")
