# =============================================================================
# R2 · scWGCNA —— STAGE 9: green module deep dive
#   S3 = green module kME ranking, labeling top hubs and PLAUR.
#   S4 = green module GO BP enrichment.
# -----------------------------------------------------------------------------
# Read-only plotting/table step: reads module csv only, does not recalculate the
# network, modify objects, or save objects.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggrepel)
})

OUT_DIR <- "<DATA_ROOT>/项目/分型/修稿杠生信/重新分析/R2_scWGCNA"
MODULES_CSV <- file.path(OUT_DIR, "tables/stage5_modules_gene_assignment_kME.csv")
fig_dir <- file.path(OUT_DIR, "figures")
tab_dir <- file.path(OUT_DIR, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

MOD <- "green"
KME_COL <- paste0("kME_", MOD)
HIGHLIGHT <- c("PLAUR")
TOP_HUB_N <- 15

modules <- read.csv(MODULES_CSV, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(KME_COL %in% colnames(modules))

g <- modules %>%
  dplyr::filter(module == MOD) %>%
  dplyr::arrange(desc(.data[[KME_COL]])) %>%
  dplyr::mutate(rank = dplyr::row_number(), kME = .data[[KME_COL]])
n_g <- nrow(g)
cat(sprintf("green module: %d genes\n", n_g))

for (hg in HIGHLIGHT) {
  r <- g$rank[g$gene_name == hg]
  k <- g$kME[g$gene_name == hg]
  if (length(r)) {
    cat(sprintf("  %s: kME=%.3f rank=%d/%d\n", hg, k, r, n_g))
  } else {
    cat(sprintf("  %s: not found in green module\n", hg))
  }
}

write.csv(
  g[, c("gene_name", "kME", "rank")],
  file.path(tab_dir, "stage9_green_kME_ranked.csv"),
  row.names = FALSE
)

top <- g %>% dplyr::slice_head(n = TOP_HUB_N)
foc <- g %>% dplyr::filter(gene_name %in% HIGHLIGHT)
lab <- dplyr::bind_rows(top, foc) %>% dplyr::distinct(gene_name, .keep_all = TRUE)

p_s3 <- ggplot(g, aes(rank, kME)) +
  geom_line(color = "grey75") +
  geom_point(color = "grey80", size = 0.6) +
  geom_point(data = top, color = "#20854E", size = 2) +
  geom_point(data = foc, color = "#BC3C29", size = 3) +
  ggrepel::geom_text_repel(
    data = lab,
    aes(label = gene_name),
    size = 3,
    max.overlaps = 30,
    segment.size = 0.3
  ) +
  labs(
    x = "Intramodular rank (by kME)",
    y = "kME (green module)",
    title = "Green (MES-V mesenchymal) module connectivity",
    subtitle = "green = top hubs; red = PLAUR (member, not hub)"
  ) +
  theme_classic(base_size = 12)
ggsave(file.path(fig_dir, "S3_green_kME_ranking.pdf"), p_s3, width = 7, height = 5)

gene_map <- clusterProfiler::bitr(
  g$gene_name,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)
genes_entrez <- unique(gene_map$ENTREZID)

ego <- clusterProfiler::enrichGO(
  genes_entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE
)
ego_df <- as.data.frame(ego)
write.csv(ego_df, file.path(tab_dir, "stage9_green_GO_BP.csv"), row.names = FALSE)
cat(sprintf("\ngreen GO BP enriched terms (q<0.1): %d\n", nrow(ego_df)))

if (nrow(ego_df) > 0) {
  top_go <- ego_df %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::mutate(
      Description = factor(Description, levels = rev(Description)),
      GeneRatio_num = vapply(strsplit(GeneRatio, "/"), function(x) {
        as.numeric(x[1]) / as.numeric(x[2])
      }, numeric(1))
    )

  p_s4 <- ggplot(top_go, aes(GeneRatio_num, Description, size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#BC3C29", high = "#0072B5", name = "p.adjust") +
    labs(
      x = "Gene ratio",
      y = NULL,
      size = "Count",
      title = "Green module GO (Biological Process)"
    ) +
    theme_bw(base_size = 11)
  ggsave(file.path(fig_dir, "S4_green_GO_BP.pdf"), p_s4, width = 7.5, height = 5)
} else {
  cat("No significant GO BP term at q < 0.1.\n")
}

cat("\nSTAGE 9 done. S3/S4 and source data written; no object saved.\n")
