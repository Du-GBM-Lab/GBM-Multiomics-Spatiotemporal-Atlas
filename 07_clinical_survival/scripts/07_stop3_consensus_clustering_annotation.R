# R4 07 / STOP3: consensus clustering + k evidence + cluster annotation.
# Input: HGNC-harmonized ssGSEA scores and old manuscript CGGA HGG clinical data.
# No KM, Cox, or survival modeling in this script.

suppressPackageStartupMessages({
  library(ConsensusClusterPlus)
})

root <- getwd()
raw_dir <- file.path(root, "data", "raw", "CGGA_old_manuscript")
proc_dir <- file.path(root, "data", "processed")
tab_dir <- file.path(root, "tables")
fig_dir <- file.path(root, "figures")
cons_dir <- file.path(root, "consensus")
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cons_dir, recursive = TRUE, showWarnings = FALSE)

SUBTYPE_LEVELS <- c("NPC_P", "OPC_M", "MES_V", "MES_I")
K_REVIEW <- 4

score_path <- file.path(proc_dir, "ssgsea_scores_harmonized.rds")
hgg_path <- file.path(raw_dir, "hgg_data.rds")
stopifnot(file.exists(score_path), file.exists(hgg_path))

es <- readRDS(score_path)
es <- as.matrix(es)

# Expected current shape is 4 subtypes x 504 samples. Support either orientation.
if (all(SUBTYPE_LEVELS %in% rownames(es))) {
  es_subtype_sample <- es[SUBTYPE_LEVELS, , drop = FALSE]
  es_sample_subtype <- t(es_subtype_sample)
} else if (all(SUBTYPE_LEVELS %in% colnames(es))) {
  es_sample_subtype <- es[, SUBTYPE_LEVELS, drop = FALSE]
  es_subtype_sample <- t(es_sample_subtype)
} else {
  stop("Cannot identify subtype score orientation.")
}

hgg <- readRDS(hgg_path)
clin <- hgg$clinical
stopifnot("CGGA_ID" %in% colnames(clin))
clin <- clin[match(rownames(es_sample_subtype), as.character(clin$CGGA_ID)), , drop = FALSE]
stopifnot(all(rownames(es_sample_subtype) == as.character(clin$CGGA_ID)))

cat("samples:", nrow(es_sample_subtype), "\n")
cat("features:", paste(colnames(es_sample_subtype), collapse = ", "), "\n")
cat("ConsensusClusterPlus:", as.character(utils::packageVersion("ConsensusClusterPlus")), "\n")

# Cluster samples using z-scored subtype scores.
z_sample_subtype <- scale(es_sample_subtype)
cc_input <- t(z_sample_subtype) # features x samples, as required by ConsensusClusterPlus.

set.seed(1)
old_wd <- getwd()
setwd(cons_dir)
cc <- ConsensusClusterPlus(
  as.matrix(cc_input),
  maxK = 6,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "hc",
  distance = "euclidean",
  seed = 1,
  title = "STOP3_R4_consensus",
  plot = "pdf",
  writeTable = TRUE,
  verbose = FALSE
)
setwd(old_wd)
saveRDS(cc, file.path(proc_dir, "STOP3_consensus_cluster_results_k2to6.rds"))

pac <- data.frame(
  k = 2:6,
  PAC = sapply(2:6, function(k) {
    M <- cc[[k]]$consensusMatrix
    v <- M[lower.tri(M)]
    mean(v > 0.1 & v < 0.9)
  })
)
write.csv(pac, file.path(tab_dir, "STOP3_consensus_PAC_k2to6.csv"), row.names = FALSE)
cat("\n== PAC (lower = cleaner; use with CDF/delta plots) ==\n")
print(pac)

cluster_tables <- list()
for (k in 2:6) {
  cl <- factor(paste0("C", cc[[k]]$consensusClass))
  cluster_tables[[paste0("k", k)]] <- data.frame(
    CGGA_ID = rownames(es_sample_subtype),
    k = k,
    cluster = cl,
    stringsAsFactors = FALSE
  )
}
all_clusters <- do.call(rbind, cluster_tables)
write.csv(all_clusters, file.path(tab_dir, "STOP3_cluster_labels_k2to6.csv"), row.names = FALSE)

cl <- factor(paste0("C", cc[[K_REVIEW]]$consensusClass))
ann <- data.frame(
  CGGA_ID = rownames(es_sample_subtype),
  cluster = cl,
  es_sample_subtype,
  clin,
  check.names = FALSE
)
saveRDS(ann, file.path(proc_dir, sprintf("cluster_labels_k%d_review.rds", K_REVIEW)))
write.csv(ann, file.path(tab_dir, sprintf("STOP3_cluster_labels_k%d_review.csv", K_REVIEW)), row.names = FALSE)

fp_mean <- aggregate(es_sample_subtype, list(cluster = cl), mean)
fp_median <- aggregate(es_sample_subtype, list(cluster = cl), median)
write.csv(fp_mean, file.path(tab_dir, sprintf("STOP3_k%d_cluster_subtype_score_mean_fingerprint.csv", K_REVIEW)), row.names = FALSE)
write.csv(fp_median, file.path(tab_dir, sprintf("STOP3_k%d_cluster_subtype_score_median_fingerprint.csv", K_REVIEW)), row.names = FALSE)

cat("\n== k=", K_REVIEW, " cluster subtype score mean fingerprint ==\n", sep = "")
fp_mean_print <- fp_mean
fp_mean_print[, SUBTYPE_LEVELS] <- round(fp_mean_print[, SUBTYPE_LEVELS], 3)
print(fp_mean_print)
cat("\n== k=", K_REVIEW, " cluster subtype score median fingerprint ==\n", sep = "")
fp_median_print <- fp_median
fp_median_print[, SUBTYPE_LEVELS] <- round(fp_median_print[, SUBTYPE_LEVELS], 3)
print(fp_median_print)

clinical_fields <- c(
  "IDH_mutation_status",
  "Grade",
  "X1p19q_codeletion_status",
  "MGMTp_methylation_status",
  "PRS_type",
  "Histology",
  "Gender",
  "Radio_status..treated.1.un.treated.0.",
  "Chemo_status..TMZ.treated.1.un.treated.0."
)
clinical_fields <- intersect(clinical_fields, colnames(clin))

xtab_rows <- data.frame()
for (field in clinical_fields) {
  tab <- table(cluster = cl, value = clin[[field]], useNA = "ifany")
  cat("\n== k=", K_REVIEW, " cluster x ", field, " ==\n", sep = "")
  print(tab)
  df <- as.data.frame(tab, stringsAsFactors = FALSE)
  colnames(df) <- c("cluster", "value", "n")
  df$field <- field
  xtab_rows <- rbind(xtab_rows, df[, c("field", "cluster", "value", "n")])
}
write.csv(xtab_rows, file.path(tab_dir, sprintf("STOP3_k%d_cluster_clinical_xtabs.csv", K_REVIEW)), row.names = FALSE)

cluster_sizes <- data.frame(cluster = names(table(cl)), n = as.integer(table(cl)))
write.csv(cluster_sizes, file.path(tab_dir, sprintf("STOP3_k%d_cluster_sizes.csv", K_REVIEW)), row.names = FALSE)

# Simple review heatmap: four scores x samples ordered by review clusters.
pdf(file.path(fig_dir, sprintf("STOP3_k%d_subtype_score_heatmap_review.pdf", K_REVIEW)), width = 9, height = 4)
ord <- order(cl)
heatmap(
  t(z_sample_subtype[ord, SUBTYPE_LEVELS, drop = FALSE]),
  Rowv = NA,
  Colv = NA,
  scale = "none",
  col = colorRampPalette(c("#0072B5", "white", "#BC3C29"))(50),
  labCol = FALSE,
  margins = c(2, 6),
  main = paste0("Subtype ssGSEA score heatmap, k=", K_REVIEW)
)
dev.off()

cat("\n\n########## STOP3 REPORT ##########\n")
cat("1. PAC table printed above and saved to STOP3_consensus_PAC_k2to6.csv; consensus CDF/delta PDFs are in consensus/.\n")
cat("2. k=", K_REVIEW, " subtype score fingerprints are printed above and saved in tables/.\n", sep = "")
cat("3. k=", K_REVIEW, " clinical cross-tabs include IDH/grade/1p19q/MGMT and are printed above.\n", sep = "")
cat("4. k=", K_REVIEW, " PRS_type cross-tab is printed above; inspect recurrent enrichment before survival.\n", sep = "")
cat("=> STOP3 complete. Do not run survival until k and cluster interpretation are reviewed.\n")
