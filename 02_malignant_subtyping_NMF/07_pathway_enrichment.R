# 05_恶性细胞分亚群与Neftel对照/07_pathway_enrichment.R
# Marker-list over-representation analysis for subtype naming validation.
# Input: 05b subtype markers. Output: GO BP, KEGG, Hallmark enrichment tables and Hallmark dotplot.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(dplyr)
  library(ggplot2)
  library(qs2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(enrichplot)
})

set.seed(42)

proj <- "05_恶性细胞分亚群与Neftel对照"
mk_csv <- file.path(proj, "tables", "05b_subtype_markers_all.csv")
audit_csv <- file.path(proj, "tables", "06_subtype_naming_audit.csv")
tab_dir <- file.path(proj, "tables")
fig_dir <- file.path(proj, "figures")
out_dir <- file.path(proj, "outputs")

dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")
safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

get_hallmark_term2gene <- function() {
  hm <- tryCatch(
    msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
    error = function(e) {
      msg("msigdbr collection='H' failed; retrying deprecated category='H':", conditionMessage(e))
      msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    }
  )

  gene_col <- intersect(c("entrez_gene", "ncbi_gene", "ncbi_gene_id", "db_ncbi_gene"), colnames(hm))[1]
  if (is.na(gene_col) || !"gs_name" %in% colnames(hm)) {
    stop("Cannot identify Hallmark TERM2GENE columns in msigdbr output.", call. = FALSE)
  }

  hm |>
    transmute(
      gs_name = .data$gs_name,
      entrez_gene = as.character(.data[[gene_col]])
    ) |>
    filter(!is.na(entrez_gene), entrez_gene != "") |>
    distinct()
}

as_result_df <- function(x, subtype, original_cluster, collection) {
  if (is.null(x) || nrow(x@result) == 0) {
    return(NULL)
  }
  x@result |>
    mutate(
      subtype = subtype,
      original_cluster = original_cluster,
      collection = collection,
      .before = 1
    )
}

msg("Reading marker table:", mk_csv)
mk <- read.csv(mk_csv, stringsAsFactors = FALSE)
audit <- read.csv(audit_csv, stringsAsFactors = FALSE)

fc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(mk))[1]
if (is.na(fc_col)) {
  stop("No avg_log2FC / avg_logFC column found in marker table.", call. = FALSE)
}
required_cols <- c("cluster", "gene", "p_val_adj", fc_col)
missing_cols <- setdiff(required_cols, colnames(mk))
if (length(missing_cols) > 0) {
  stop("Missing marker columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

name_map <- setNames(audit$final_label, audit$original_label)
mk$subtype <- unname(name_map[as.character(mk$cluster)])
if (any(is.na(mk$subtype))) {
  stop("Some marker clusters cannot be mapped to final subtype names.", call. = FALSE)
}

mk_filt <- mk |>
  filter(p_val_adj < 0.05, .data[[fc_col]] > 0.5, !is.na(gene), gene != "")

cat("== Genes per subtype after marker filter ==\n")
print(table(mk_filt$subtype))

msg("Mapping SYMBOL to ENTREZID.")
gene_map <- bitr(
  unique(mk_filt$gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
) |>
  distinct(SYMBOL, ENTREZID)

mk_filt <- mk_filt |>
  inner_join(gene_map, by = c("gene" = "SYMBOL")) |>
  mutate(ENTREZID = as.character(ENTREZID))

cat("\n== Genes per subtype after ENTREZ mapping ==\n")
print(table(mk_filt$subtype))

msg("Loading MSigDB Hallmark TERM2GENE.")
hallmark <- get_hallmark_term2gene()

subtype_levels <- c("NPC-Cycling", "OPC-like", "MES-Perivascular", "MES-Inflammatory")
subtypes <- subtype_levels[subtype_levels %in% unique(mk_filt$subtype)]
all_results <- list()
all_result_tables <- list()

for (st in subtypes) {
  original_cluster <- audit$original_label[match(st, audit$final_label)]
  st_safe <- safe_name(st)
  cat(sprintf("\n>>> Enrichment for %s (%s)\n", st, original_cluster))

  gene_st <- mk_filt |>
    filter(subtype == st) |>
    pull(ENTREZID) |>
    unique()

  cat(sprintf("  Input ENTREZ genes: %d\n", length(gene_st)))
  if (length(gene_st) < 10) {
    cat(sprintf("  [skip] only %d mapped genes\n", length(gene_st)))
    next
  }

  go_bp <- tryCatch(
    enrichGO(
      gene = gene_st,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    ),
    error = function(e) {
      msg("GO BP failed for", st, ":", conditionMessage(e))
      NULL
    }
  )

  kegg <- tryCatch(
    enrichKEGG(
      gene = gene_st,
      organism = "hsa",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    ),
    error = function(e) {
      msg("KEGG failed for", st, ":", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }

  hm <- tryCatch(
    enricher(
      gene = gene_st,
      TERM2GENE = hallmark,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    ),
    error = function(e) {
      msg("Hallmark failed for", st, ":", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(hm) && nrow(hm@result) > 0) {
    hm <- setReadable(hm, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }

  all_results[[st]] <- list(GO_BP = go_bp, KEGG = kegg, Hallmark = hm)

  go_df <- as_result_df(go_bp, st, original_cluster, "GO_BP")
  kegg_df <- as_result_df(kegg, st, original_cluster, "KEGG")
  hm_df <- as_result_df(hm, st, original_cluster, "Hallmark")
  all_result_tables <- c(all_result_tables, list(go_df, kegg_df, hm_df))

  if (!is.null(go_df)) {
    write.csv(go_df, file.path(tab_dir, sprintf("07_%s_GO_BP.csv", st_safe)), row.names = FALSE)
  }
  if (!is.null(kegg_df)) {
    write.csv(kegg_df, file.path(tab_dir, sprintf("07_%s_KEGG.csv", st_safe)), row.names = FALSE)
  }
  if (!is.null(hm_df)) {
    write.csv(hm_df, file.path(tab_dir, sprintf("07_%s_Hallmark.csv", st_safe)), row.names = FALSE)
  }

  if (!is.null(hm_df)) {
    cat("  Top 5 Hallmark:\n")
    print(head(hm_df[, c("ID", "Description", "Count", "p.adjust")], 5))
  } else {
    cat("  Top 5 Hallmark: none passed cutoff.\n")
  }

  if (!is.null(go_df)) {
    cat("  Top 5 GO BP:\n")
    print(head(go_df[, c("ID", "Description", "Count", "p.adjust")], 5))
  } else {
    cat("  Top 5 GO BP: none passed cutoff.\n")
  }
}

msg("Saving enrichment objects.")
qs2::qs_save(all_results, file.path(out_dir, "07_pathway_enrichment_results.qs2"))

combined_results <- bind_rows(all_result_tables)
if (nrow(combined_results) > 0) {
  write.csv(combined_results, file.path(tab_dir, "07_all_enrichment_results.csv"), row.names = FALSE)
}

hm_top <- bind_rows(lapply(names(all_results), function(st) {
  r <- all_results[[st]]$Hallmark
  if (is.null(r) || nrow(r@result) == 0) {
    return(NULL)
  }
  r@result |>
    slice_min(p.adjust, n = 5, with_ties = FALSE) |>
    mutate(subtype = st)
}))

if (nrow(hm_top) > 0) {
  hm_top <- hm_top |>
    mutate(
      ID_clean = gsub("^HALLMARK_", "", ID),
      ID_clean = gsub("_", " ", ID_clean),
      subtype = factor(subtype, levels = subtype_levels),
      term_order = reorder(ID_clean, -log10(p.adjust))
    )

  p_dot <- ggplot(
    hm_top,
    aes(x = subtype, y = term_order, size = Count, color = -log10(p.adjust))
  ) +
    geom_point(alpha = 0.90) +
    scale_color_gradient(low = "grey78", high = "firebrick3") +
    labs(
      title = "Top 5 Hallmark terms per subtype",
      x = NULL,
      y = NULL,
      size = "Gene count",
      color = "-log10(p.adj)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "right"
    )

  ggsave(
    file.path(fig_dir, "07_hallmark_top5_dotplot.pdf"),
    p_dot,
    width = 8,
    height = 6,
    device = cairo_pdf
  )
  ggsave(
    file.path(fig_dir, "07_hallmark_top5_dotplot.png"),
    p_dot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

go_top <- bind_rows(lapply(names(all_results), function(st) {
  r <- all_results[[st]]$GO_BP
  if (is.null(r) || nrow(r@result) == 0) {
    return(NULL)
  }
  r@result |>
    slice_min(p.adjust, n = 5, with_ties = FALSE) |>
    mutate(subtype = st)
}))

if (nrow(go_top) > 0) {
  go_top <- go_top |>
    mutate(
      subtype = factor(subtype, levels = subtype_levels),
      term_order = reorder(Description, -log10(p.adjust))
    )

  p_go <- ggplot(
    go_top,
    aes(x = subtype, y = term_order, size = Count, color = -log10(p.adjust))
  ) +
    geom_point(alpha = 0.90) +
    scale_color_gradient(low = "grey78", high = "steelblue4") +
    labs(
      title = "Top 5 GO BP terms per subtype",
      x = NULL,
      y = NULL,
      size = "Gene count",
      color = "-log10(p.adj)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "right"
    )

  ggsave(
    file.path(fig_dir, "07_GO_BP_top5_dotplot.pdf"),
    p_go,
    width = 9,
    height = 7,
    device = cairo_pdf
  )
  ggsave(
    file.path(fig_dir, "07_GO_BP_top5_dotplot.png"),
    p_go,
    width = 9,
    height = 7,
    dpi = 300
  )
}

cat("\nDone. Outputs:\n")
cat("- tables/07_*_GO_BP.csv\n")
cat("- tables/07_*_KEGG.csv\n")
cat("- tables/07_*_Hallmark.csv\n")
cat("- tables/07_all_enrichment_results.csv\n")
cat("- figures/07_hallmark_top5_dotplot.pdf\n")
cat("- figures/07_GO_BP_top5_dotplot.pdf\n")
