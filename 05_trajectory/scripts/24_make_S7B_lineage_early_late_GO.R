options(stringsAsFactors = FALSE)
set.seed(42)

if (basename(getwd()) == "06_恶性细胞拟时序") {
  project_root <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(file.path(project_root, "06_恶性细胞拟时序"))

dir.create("figures/final_panels/supplement", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/final_panels/source_data/supplement", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  .libPaths(c("<R_LIBS>", "<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

log_file <- "logs/24_make_S7B_lineage_early_late_GO.log"
if (file.exists(log_file)) file.remove(log_file)
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(...))
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

lineage_labels <- c(
  lineage1 = "L1: NPC-P to MES-I",
  lineage2 = "L2: NPC-P to MES-V",
  lineage3 = "L3: NPC-P to OPC-M"
)
timing_labels <- c(
  early_high = "Early-high",
  late_high = "Late-high"
)

assoc <- read_csv("tables/tradeseq_association_per_lineage.csv", show_col_types = FALSE) |>
  filter(.data$p_adj_BH < 0.05)

binned <- read_csv("figures/source_data/08_binned_expression_for_patterns.csv", show_col_types = FALSE) |>
  inner_join(assoc, by = c("gene", "lineage"))

gene_metrics <- binned |>
  group_by(.data$gene, .data$lineage) |>
  summarise(
    early_mean = mean(.data$expr_mean[.data$bin <= 4], na.rm = TRUE),
    late_mean = mean(.data$expr_mean[.data$bin >= 17], na.rm = TRUE),
    min_expr = min(.data$expr_mean, na.rm = TRUE),
    max_expr = max(.data$expr_mean, na.rm = TRUE),
    dynamic_range = .data$max_expr - .data$min_expr,
    early_late_delta = .data$late_mean - .data$early_mean,
    abs_delta = abs(.data$early_late_delta),
    waldStat = dplyr::first(.data$waldStat),
    p_adj_BH = dplyr::first(.data$p_adj_BH),
    pattern = dplyr::first(.data$pattern),
    .groups = "drop"
  ) |>
  mutate(
    lineage_label = recode(.data$lineage, !!!lineage_labels),
    timing = if_else(.data$early_late_delta < 0, "early_high", "late_high"),
    timing_label = recode(.data$timing, !!!timing_labels),
    visual_score = .data$abs_delta * 0.55 + .data$dynamic_range * 0.35 + log1p(.data$waldStat) * 0.10,
    gene_set = paste(.data$lineage, .data$timing, sep = "_")
  )

top_n <- 60
selected_genes <- gene_metrics |>
  group_by(.data$lineage, .data$lineage_label, .data$timing, .data$timing_label, .data$gene_set) |>
  arrange(desc(.data$visual_score), desc(.data$dynamic_range), .by_group = TRUE) |>
  slice_head(n = top_n) |>
  ungroup()

gene_set_summary <- selected_genes |>
  count(.data$lineage, .data$lineage_label, .data$timing, .data$timing_label, .data$gene_set, name = "n_symbol_selected")

universe <- read_csv("tables/tradeseq_gene_universe.csv", show_col_types = FALSE) |>
  pull(.data$gene) |>
  unique()

gene_map <- clusterProfiler::bitr(unique(selected_genes$gene), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
universe_map <- clusterProfiler::bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = TRUE)
enrich_universe <- unique(universe_map$ENTREZID)

enrichment <- tibble()
for (gs in unique(selected_genes$gene_set)) {
  genes <- selected_genes$gene[selected_genes$gene_set == gs]
  entrez <- unique(gene_map$ENTREZID[match(genes, gene_map$SYMBOL)])
  entrez <- entrez[!is.na(entrez)]
  if (length(entrez) < 8) {
    log_msg("Skip", gs, "because mapped ENTREZ genes <", 8)
    next
  }
  ego <- tryCatch(
    clusterProfiler::enrichGO(
      gene = entrez,
      universe = enrich_universe,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.2,
      qvalueCutoff = 0.2,
      minGSSize = 5,
      readable = TRUE
    ),
    error = function(e) {
      log_msg("enrichGO failed for", gs, ":", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    enrichment <- bind_rows(enrichment, as.data.frame(ego) |> mutate(gene_set = gs))
  }
}

plot_terms <- enrichment |>
  left_join(gene_set_summary, by = "gene_set") |>
  mutate(
    neglog10_BH = -log10(pmax(.data$p.adjust, 1e-300)),
    lineage_label = factor(.data$lineage_label, levels = lineage_labels),
    timing_label = factor(.data$timing_label, levels = c("Early-high", "Late-high")),
    term_label = stringr::str_wrap(.data$Description, width = 32)
  ) |>
  group_by(.data$gene_set) |>
  arrange(.data$p.adjust, .by_group = TRUE) |>
  slice_head(n = 5) |>
  ungroup()

term_levels <- plot_terms |>
  arrange(.data$lineage_label, .data$timing_label, .data$p.adjust) |>
  pull(.data$term_label) |>
  unique() |>
  rev()
plot_terms <- plot_terms |>
  mutate(term_label = factor(.data$term_label, levels = term_levels))

readr::write_csv(gene_metrics, "figures/final_panels/source_data/supplement/S7B_lineage_early_late_GO_gene_metrics.csv")
readr::write_csv(selected_genes, "figures/final_panels/source_data/supplement/S7B_lineage_early_late_GO_selected_genes.csv")
readr::write_csv(gene_set_summary, "figures/final_panels/source_data/supplement/S7B_lineage_early_late_GO_gene_set_summary.csv")
readr::write_csv(enrichment, "figures/final_panels/source_data/supplement/S7B_lineage_early_late_GO_enrichment_all.csv")
readr::write_csv(plot_terms, "figures/final_panels/source_data/supplement/S7B_lineage_early_late_GO_plot_terms.csv")
readr::write_csv(
  tibble(
    background_type = "tradeSeq_gene_universe",
    n_symbol = length(unique(universe)),
    n_entrez = length(enrich_universe),
    top_n_per_lineage_timing = top_n
  ),
  "figures/final_panels/source_data/supplement/S7B_lineage_early_late_GO_background.csv"
)

if (nrow(plot_terms) == 0) {
  stop("No GO terms available for plotting.")
}

panel_labels <- gene_set_summary |>
  mutate(panel_label = paste0(timing_label, "\n", "n=", n_symbol_selected)) |>
  dplyr::select(lineage_label, timing_label, panel_label)

plot_terms <- plot_terms |>
  left_join(panel_labels, by = c("lineage_label", "timing_label"))

gene_set_order <- gene_set_summary |>
  mutate(
    lineage_num = as.integer(str_extract(.data$lineage, "[0-9]+")),
    timing_num = if_else(.data$timing == "early_high", 1L, 2L)
  ) |>
  arrange(.data$lineage_num, .data$timing_num) |>
  mutate(
    x_pos = row_number(),
    x_label = paste0(str_extract(.data$lineage_label, "^L[0-9]"), "\n", str_replace(.data$timing_label, "-high", ""), "\nn=", .data$n_symbol_selected)
  )

plot_terms <- plot_terms |>
  left_join(gene_set_order |> dplyr::select(gene_set, x_pos, x_label), by = "gene_set")

term_order <- plot_terms |>
  group_by(.data$term_label) |>
  summarise(
    first_x = min(.data$x_pos, na.rm = TRUE),
    best_bh = min(.data$p.adjust, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(.data$first_x, .data$best_bh, .data$term_label) |>
  mutate(y_pos = row_number())

plot_terms <- plot_terms |>
  left_join(term_order |> dplyr::select(term_label, y_pos), by = "term_label")

empty_sets <- gene_set_order |>
  filter(!.data$gene_set %in% unique(plot_terms$gene_set)) |>
  mutate(y_pos = max(term_order$y_pos, na.rm = TRUE) * 0.52)

lineage_bands <- tibble(
  xmin = c(0.5, 2.5, 4.5),
  xmax = c(2.5, 4.5, 6.5),
  label = c("L1: NPC-P to MES-I", "L2: NPC-P to MES-V", "L3: NPC-P to OPC-M")
)

p <- ggplot() +
  geom_rect(
    data = lineage_bands,
    aes(xmin = .data$xmin, xmax = .data$xmax, ymin = -Inf, ymax = Inf),
    fill = c("#F7F9FC", "#FFFFFF", "#F7F9FC"),
    color = NA
  ) +
  geom_vline(xintercept = c(2.5, 4.5), linetype = "22", color = "grey68", linewidth = 0.28) +
  geom_point(
    data = plot_terms,
    aes(x = .data$x_pos, y = .data$y_pos, size = .data$Count, fill = .data$neglog10_BH),
    shape = 21, color = "grey22", stroke = 0.18, alpha = 0.96
  ) +
  geom_text(
    data = empty_sets,
    aes(x = .data$x_pos, y = .data$y_pos),
    label = "No enriched\nGO term\n(BH < 0.2)",
    size = 2.0,
    lineheight = 0.9,
    color = "grey45"
  ) +
  annotate(
    "text",
    x = c(1.5, 3.5, 5.5),
    y = max(term_order$y_pos, na.rm = TRUE) + 1.3,
    label = lineage_bands$label,
    size = 2.3,
    fontface = "bold",
    color = "grey15"
  ) +
  scale_x_continuous(
    breaks = gene_set_order$x_pos,
    labels = gene_set_order$x_label,
    limits = c(0.5, 6.5),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = term_order$y_pos,
    labels = term_order$term_label,
    limits = c(0.5, max(term_order$y_pos, na.rm = TRUE) + 1.8),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_gradientn(colors = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5"), name = "-log10(BH)") +
  scale_size_continuous(range = c(2.2, 7.0), name = "Genes") +
  labs(
    title = "S7B. GO enrichment of early- and late-high trajectory driver genes",
    subtitle = "Top 60 dynamic genes per lineage/time group except L3 late-high (n=27); background = 600-gene tradeSeq universe",
    x = NULL,
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 7) +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 7, color = "grey30"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey91", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 5.8, color = "grey15"),
    axis.text.x = element_text(size = 6.4, face = "bold", lineheight = 0.92),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    plot.margin = margin(8, 8, 6, 5.5)
  )

pdf_path <- "figures/final_panels/supplement/S7B_lineage_early_late_GO_dotplot.pdf"
png_path <- "figures/final_panels/supplement/S7B_lineage_early_late_GO_dotplot.png"
ggsave(pdf_path, p, width = 8.4, height = 6.4, useDingbats = FALSE)
ggsave(png_path, p, width = 8.4, height = 6.4, dpi = 400, type = "cairo")

log_msg("Wrote", pdf_path)
log_msg("Gene set sizes:", paste(gene_set_summary$gene_set, gene_set_summary$n_symbol_selected, collapse = "; "))
log_msg("STOP")
