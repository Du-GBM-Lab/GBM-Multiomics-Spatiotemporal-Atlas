# STOP6b: main Panel A ordered heatmap.
# Rows = MES_V ECMcore genes; columns = CGGA HGG samples ordered by MES_V score.
# Column clustering is explicitly disabled because R4 claims a continuous landscape, not discrete clusters.

root <- getwd()
source(file.path(root, "scripts", "R4_helpers.R"))

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
})

tab_dir <- file.path(root, "tables")
fig_dir <- file.path(root, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

es <- as.data.frame(readRDS(file.path(root, "data", "processed", "ssgsea_scores_harmonized.rds")))
if (all(PROGS %in% rownames(es))) es <- as.data.frame(t(as.matrix(es)))
hgg <- readRDS(file.path(root, "data", "raw", "CGGA_old_manuscript", "hgg_data.rds"))
common_ids <- intersect(rownames(es), hgg$clinical$CGGA_ID)
es <- es[common_ids, , drop = FALSE]
clin <- hgg$clinical[match(common_ids, hgg$clinical$CGGA_ID), , drop = FALSE]
stopifnot(all(rownames(es) == clin$CGGA_ID))

expr_df <- hgg$expression
gene_col <- "Gene_Name"
stopifnot(gene_col %in% colnames(expr_df))
expr <- as.matrix(expr_df[, setdiff(colnames(expr_df), gene_col), drop = FALSE])
mode(expr) <- "numeric"
rownames(expr) <- expr_df[[gene_col]]
expr <- expr[, rownames(es), drop = FALSE]

ECM21 <- c("ACTA2", "ANXA2", "CCN2", "COL1A2", "COL4A1", "COL4A2", "COL8A1",
           "FN1", "FSTL1", "IGFBP7", "LUM", "MMP14", "MYL9", "PCOLCE",
           "SPARCL1", "TAGLN", "TFPI", "TGFB1I1", "TGM2", "TPM1", "TPM2")

expr_cur <- to_current(rownames(expr))
ecm_cur <- to_current(ECM21)
rows <- !is.na(expr_cur) & expr_cur %in% ecm_cur
m <- expr[rows, , drop = FALSE]
g <- expr_cur[rows]
sm <- rowsum(m, group = g, reorder = FALSE)
counts <- as.vector(table(factor(g, levels = rownames(sm))))
m <- sweep(sm, 1, counts, "/")
m <- m[intersect(ecm_cur, rownames(m)), , drop = FALSE]

# Log transform before row z-score to reduce leverage of extreme RSEM-like abundance values.
m_log <- log2(m + 1)
mz <- t(scale(t(m_log)))
mz[is.na(mz)] <- 0

ord <- order(es$MES_V)
mz <- mz[, ord, drop = FALSE]
cl <- clin[ord, , drop = FALSE]
sc <- es[ord, , drop = FALSE]

vital <- factor(ifelse(cl$Censor..alive.0..dead.1. == 1, "dead", "alive"), levels = c("alive", "dead"))
grade <- factor(cl$Grade, levels = c("WHO III", "WHO IV"))
idh <- factor(cl$IDH_mutation_status, levels = c("Mutant", "Wildtype"))
codel <- factor(cl$X1p19q_codeletion_status, levels = c("Codel", "Non-codel"))
prs <- factor(cl$PRS_type, levels = c("Primary", "Secondary", "Recurrent"))

cont <- function(x) {
  circlize::colorRamp2(as.numeric(stats::quantile(x, c(0.05, 0.5, 0.95), na.rm = TRUE)),
                       c("#2c7bb6", "white", "#d7191c"))
}

ann <- ComplexHeatmap::HeatmapAnnotation(
  MES_V = sc$MES_V,
  IDH = idh,
  Grade = grade,
  codel = codel,
  PRS = prs,
  Vital = vital,
  col = list(
    MES_V = cont(sc$MES_V),
    IDH = c(Mutant = "#377eb8", Wildtype = "#e41a1c"),
    Grade = c("WHO III" = "#fdae61", "WHO IV" = "#a50026"),
    codel = c(Codel = "#4daf4a", "Non-codel" = "grey80"),
    PRS = c(Primary = "#984ea3", Secondary = "grey50", Recurrent = "#ff7f00"),
    Vital = c(alive = "grey85", dead = "black")
  ),
  na_col = "white",
  annotation_name_side = "left",
  simple_anno_size = grid::unit(3.5, "mm")
)

ht <- ComplexHeatmap::Heatmap(
  mz,
  name = "z(expr)",
  top_annotation = ann,
  col = circlize::colorRamp2(c(-2, 0, 2), c("#2166ac", "white", "#b2182b")),
  cluster_columns = FALSE,
  show_column_names = FALSE,
  cluster_rows = TRUE,
  row_names_gp = grid::gpar(fontsize = 8),
  column_title = "CGGA HGG (n=504), samples ordered by MES_V score"
)

pdf(file.path(fig_dir, "Fig_R4_PanelA_ordered_heatmap.pdf"), width = 10, height = 6)
ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

write.csv(mz, file.path(tab_dir, "STOP6b_ECMcore21_ordered_heatmap_zmatrix.csv"))
write.csv(data.frame(
  CGGA_ID = rownames(sc),
  order_index = seq_len(nrow(sc)),
  sc,
  IDH = as.character(idh),
  Grade = as.character(grade),
  codel = as.character(codel),
  PRS = as.character(prs),
  Vital = as.character(vital),
  stringsAsFactors = FALSE
), file.path(tab_dir, "STOP6b_ordered_heatmap_column_annotation.csv"), row.names = FALSE)
write.csv(data.frame(requested_gene = ECM21, current_symbol = ecm_cur,
                     in_heatmap = ecm_cur %in% rownames(mz), stringsAsFactors = FALSE),
          file.path(tab_dir, "STOP6b_ECMcore21_gene_mapping.csv"), row.names = FALSE)

cat("\n########## STOP6b REPORT ##########\n")
cat("1. Heatmap written to figures/Fig_R4_PanelA_ordered_heatmap.pdf.\n")
cat("2. Columns are ordered by MES_V ascending; cluster_columns = FALSE.\n")
cat("3. Annotation levels used:\n")
cat("   Grade:", paste(levels(grade), collapse = ", "), "\n")
cat("   codel:", paste(levels(codel), collapse = ", "), "\n")
cat("   PRS:", paste(levels(prs), collapse = ", "), "\n")
cat("4. ECMcore genes in heatmap:", nrow(mz), "/", length(unique(ecm_cur)), "\n")
cat("=> Review whether the right side shows MES_V-high ECMcore expression and adverse clinical annotations.\n")
