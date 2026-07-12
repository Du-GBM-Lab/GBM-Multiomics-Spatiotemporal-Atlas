suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
})

step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir) {
  setwd(file.path(getwd(), step_dir))
}

de_path <- "tables/10e_one_vs_rest_FindAllMarkers_wilcoxon_downsampled.csv"
out_source <- "figures/source_data/F6_ORA_source.csv"
out_gene_input <- "figures/source_data/F6_ORA_gene_input.csv"
out_go_pdf <- "figures/F6_GO_BP_per_subtype.pdf"
out_kegg_pdf <- "figures/F6_KEGG_per_subtype.pdf"
out_sanity <- "figures/source_data/F6_ORA_sanity_checks.csv"
out_session <- "figures/source_data/F6_ORA_session_info.txt"

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/source_data", showWarnings = FALSE, recursive = TRUE)

subtype_levels <- paste0("Subtype", 1:4)
subtype_colors <- c(
  "Subtype1" = "#0072B5",
  "Subtype2" = "#E18727",
  "Subtype3" = "#20854E",
  "Subtype4" = "#BC3C29"
)

de <- readr::read_csv(de_path, show_col_types = FALSE)
stopifnot(all(c("subtype_k4", "gene", "avg_log2FC", "BH_q") %in% colnames(de)))

de <- de %>%
  mutate(subtype_k4 = factor(subtype_k4, levels = subtype_levels))

marker_tbl <- de %>%
  filter(BH_q < 0.05, avg_log2FC > 0.5) %>%
  distinct(subtype_k4, gene, .keep_all = TRUE)

universe_symbols <- unique(de$gene)
universe_map <- suppressMessages(clusterProfiler::bitr(
  universe_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
))
universe_entrez <- unique(universe_map$ENTREZID)

gene_input <- marker_tbl %>%
  left_join(universe_map, by = c("gene" = "SYMBOL")) %>%
  filter(!is.na(ENTREZID)) %>%
  distinct(subtype_k4, gene, ENTREZID, .keep_all = TRUE)

readr::write_csv(gene_input, out_gene_input)

run_ora_one <- function(st, db_name) {
  genes <- gene_input %>%
    filter(subtype_k4 == st) %>%
    pull(ENTREZID) %>%
    unique()
  if (length(genes) < 10) {
    return(tibble())
  }

  res <- NULL
  if (db_name == "GO_BP") {
    res <- tryCatch(
      clusterProfiler::enrichGO(
        gene = genes,
        universe = universe_entrez,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = TRUE
      ),
      error = function(e) e
    )
  } else if (db_name == "KEGG") {
    res <- tryCatch(
      clusterProfiler::enrichKEGG(
        gene = genes,
        universe = universe_entrez,
        organism = "hsa",
        keyType = "kegg",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1
      ),
      error = function(e) e
    )
  }

  if (inherits(res, "error") || is.null(res)) {
    return(tibble(
      subtype_k4 = as.character(st),
      database = db_name,
      ora_error = ifelse(inherits(res, "error"), conditionMessage(res), "null_result")
    ))
  }

  out <- as.data.frame(res)
  if (nrow(out) == 0) {
    return(tibble())
  }
  out %>%
    as_tibble() %>%
    mutate(
      subtype_k4 = as.character(st),
      database = db_name,
      ora_error = NA_character_
    )
}

ora_all <- bind_rows(lapply(subtype_levels, function(st) {
  bind_rows(run_ora_one(st, "GO_BP"), run_ora_one(st, "KEGG"))
}))

if (!"p.adjust" %in% colnames(ora_all)) {
  stop("ORA failed for all databases; no p.adjust column returned.")
}

ora_all <- ora_all %>%
  mutate(
    pvalue = suppressWarnings(as.numeric(pvalue)),
    p.adjust = suppressWarnings(as.numeric(p.adjust)),
    qvalue = suppressWarnings(as.numeric(qvalue)),
    subtype_k4 = factor(subtype_k4, levels = subtype_levels),
    neg_log10_BH_q = -log10(p.adjust + 1e-300),
    term_short = stringr::str_trunc(Description, width = 42)
  )

readr::write_csv(ora_all, out_source)

plot_ora <- function(db_name, out_pdf, width = 9.5, height = 7.0) {
  plot_tbl <- ora_all %>%
    filter(database == db_name, !is.na(p.adjust), p.adjust < 0.05) %>%
    group_by(subtype_k4) %>%
    arrange(p.adjust, desc(neg_log10_BH_q), .by_group = TRUE) %>%
    slice_head(n = 8) %>%
    ungroup() %>%
    group_by(subtype_k4) %>%
    mutate(term_plot = reorder(term_short, neg_log10_BH_q)) %>%
    ungroup()

  if (nrow(plot_tbl) == 0) {
    warning("No significant ", db_name, " terms for plotting.")
    return(invisible(NULL))
  }

  p <- ggplot(plot_tbl, aes(x = neg_log10_BH_q, y = term_plot, fill = subtype_k4)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.2) +
    facet_wrap(~ subtype_k4, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = subtype_colors, guide = "none") +
    labs(x = expression(-log[10]("BH q")), y = NULL) +
    theme_classic(base_size = 8) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 8, face = "bold"),
      axis.text.y = element_text(size = 6.5),
      axis.text.x = element_text(size = 7),
      axis.title.x = element_text(size = 8),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      panel.spacing = unit(7, "mm"),
      plot.margin = margin(4, 4, 4, 4)
    )

  ggsave(out_pdf, p, width = width, height = height, useDingbats = FALSE)
}

plot_ora("GO_BP", out_go_pdf, width = 9.8, height = 7.2)
plot_ora("KEGG", out_kegg_pdf, width = 9.2, height = 6.8)

sanity_tbl <- tibble::tibble(
  metric = c(
    "n_marker_rows",
    "n_unique_marker_genes",
    "n_universe_symbols",
    "n_universe_entrez",
    "n_gene_input_mapped_rows",
    "GO_BP_significant_terms",
    "KEGG_significant_terms"
  ),
  value = c(
    nrow(marker_tbl),
    length(unique(marker_tbl$gene)),
    length(universe_symbols),
    length(universe_entrez),
    nrow(gene_input),
    sum(ora_all$database == "GO_BP" & ora_all$p.adjust < 0.05, na.rm = TRUE),
    sum(ora_all$database == "KEGG" & ora_all$p.adjust < 0.05, na.rm = TRUE)
  )
)
readr::write_csv(sanity_tbl, out_sanity)
writeLines(capture.output(sessionInfo()), out_session)

cat("F6 ORA completed.\n")
cat("Marker genes by subtype after Entrez mapping:\n")
print(table(gene_input$subtype_k4))
cat("Significant terms:\n")
print(ora_all %>% filter(!is.na(p.adjust), p.adjust < 0.05) %>% count(database, subtype_k4))
cat("Output GO PDF:", out_go_pdf, "\n")
cat("Output KEGG PDF:", out_kegg_pdf, "\n")
cat("Source:", out_source, "\n")
