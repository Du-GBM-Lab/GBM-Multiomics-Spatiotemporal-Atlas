#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
})

script_arg <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_arg[grep("^--file=", script_arg)][1])
scenic_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
project_root <- normalizePath(file.path(scenic_dir, "..", "..", ".."), mustWork = TRUE)
final_root <- file.path(project_root, "图片表格", "05_发育时间_TF_ATAC验证", "02_TF_regulon_SCENIC")

candidate_dir <- file.path(final_root, "正文图候选")
main_dir <- file.path(final_root, "正文图")
source_dir <- file.path(final_root, "source_data", "D_四亚型代表regulon_target_GO")
record_dir <- file.path(final_root, "记录")
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(record_dir, recursive = TRUE, showWarnings = FALSE)

grn <- read_csv(file.path(scenic_dir, "outputs", "scenic_grn_edgelist.csv"), show_col_types = FALSE) %>%
  mutate(
    TF = as.character(TF),
    target = as.character(target),
    importance = as.numeric(importance)
  ) %>%
  filter(!is.na(TF), !is.na(target), !is.na(importance))

# Multiple representative regulons per subtype. This restores the original Fig.3D logic:
# subtype-level functional domains are inferred from regulon target sets, not from one TF.
rep_tfs <- tibble(
  subtype = c(
    rep("NPC-P", 3),
    rep("OPC-M", 3),
    rep("MES-V", 3),
    rep("MES-I", 3)
  ),
  TF = c(
    "E2F1", "E2F7", "TFDP1",
    "SOX10", "SOX8", "SREBF1",
    "FOSL1", "ETS1", "BHLHE40",
    "SPI1", "IRF8", "IKZF1"
  ),
  role = c(
    "DNA replication/cell cycle", "mitotic proliferation", "E2F partner",
    "oligodendrocyte/myelin", "oligodendrocyte lineage", "lipid/myelin",
    "AP-1/MES-V candidate", "vascular/adhesion", "hypoxia/stress",
    "myeloid inflammatory", "immune regulon", "lymphoid/immune"
  )
)

available <- intersect(rep_tfs$TF, unique(grn$TF))
missing <- setdiff(rep_tfs$TF, available)
if (length(missing)) message("Missing representative TFs in GRN: ", paste(missing, collapse = ", "))
rep_tfs <- rep_tfs %>% filter(TF %in% available)

target_n_per_tf <- 300
target_tbl <- grn %>%
  inner_join(rep_tfs, by = "TF") %>%
  arrange(subtype, TF, desc(importance)) %>%
  group_by(subtype, TF) %>%
  slice_head(n = target_n_per_tf) %>%
  ungroup()

write_csv(rep_tfs, file.path(source_dir, "representative_TF_selection.csv"))
write_csv(target_tbl, file.path(source_dir, "representative_TF_top_targets.csv"))

universe <- unique(grn$target)
ego_list <- list()
for (sub in unique(rep_tfs$subtype)) {
  genes <- unique(target_tbl$target[target_tbl$subtype == sub])
  ego <- enrichGO(
    gene = genes,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    minGSSize = 5,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    readable = FALSE
  )
  df <- as.data.frame(ego)
  if (nrow(df) > 0) {
    df$subtype <- sub
    df$representative_TFs <- paste(rep_tfs$TF[rep_tfs$subtype == sub], collapse = ", ")
    ego_list[[sub]] <- df
  }
}
ego_all <- bind_rows(ego_list)
write_csv(ego_all, file.path(source_dir, "GO_BP_enrichment_all_terms.csv"))

preferred_patterns <- c(
  "NPC-P" = "DNA replication|cell cycle|mitotic|chromosome|proliferation|axonogenesis|neuron projection",
  "OPC-M" = "myelination|oligodendrocyte|glial cell differentiation|axon ensheathment|lipid|sterol",
  "MES-V" = "extracellular matrix|cell adhesion|migration|vasculature|angiogenesis|hypoxia|wound|glycolytic|glucose|carbohydrate",
  "MES-I" = "immune|inflammatory|leukocyte|cytokine|interferon|innate immune|defense response"
)

select_terms <- function(df, sub) {
  x <- df %>% filter(subtype == sub)
  if (nrow(x) == 0) return(x)
  preferred <- x %>%
    filter(p.adjust <= 0.2, str_detect(tolower(Description), preferred_patterns[sub])) %>%
    arrange(p.adjust, desc(Count))
  fallback <- x %>% filter(p.adjust <= 0.2) %>% arrange(p.adjust, desc(Count))
  bind_rows(preferred, fallback) %>%
    distinct(ID, .keep_all = TRUE) %>%
    slice_head(n = 4)
}

subtype_order <- c("NPC-P", "OPC-M", "MES-V", "MES-I")
plot_terms <- bind_rows(lapply(subtype_order, function(s) select_terms(ego_all, s)))
plot_terms <- plot_terms %>%
  mutate(
    subtype = factor(subtype, levels = subtype_order),
    Description_short = str_replace_all(Description, "regulation of ", "reg. of "),
    Description_short = str_replace_all(Description_short, "positive ", "pos. "),
    Description_short = str_replace_all(Description_short, "negative ", "neg. "),
    Description_short = str_trunc(Description_short, width = 54),
    negLog10FDR = -log10(p.adjust)
  )

term_order <- plot_terms %>%
  arrange(subtype, p.adjust, desc(Count)) %>%
  pull(Description_short) %>%
  unique()
plot_terms$Description_short <- factor(plot_terms$Description_short, levels = rev(term_order))

write_csv(plot_terms, file.path(source_dir, "GO_BP_enrichment_plot_terms.csv"))

subtype_colors <- c("NPC-P" = "#4C78A8", "OPC-M" = "#4B9A8D", "MES-V" = "#C85B4B", "MES-I" = "#8F6BB1")
tf_labels <- rep_tfs %>%
  group_by(subtype) %>%
  summarise(label = paste(TF, collapse = ", "), .groups = "drop")
x_labels <- subtype_order

p <- ggplot(plot_terms, aes(x = subtype, y = Description_short)) +
  geom_point(aes(size = Count, color = negLog10FDR), alpha = 0.96) +
  scale_color_gradient(low = "#F1E7A5", high = "#3B4A89", name = "-log10(FDR)") +
  scale_size_area(max_size = 8.3, name = "target genes") +
  scale_x_discrete(labels = x_labels) +
  labs(
    title = "Subtype representative regulon target-gene GO enrichment",
    subtitle = "Representative regulons: NPC-P (E2F1/E2F7/TFDP1), OPC-M (SOX10/SOX8/SREBF1), MES-V (FOSL1/ETS1/BHLHE40), MES-I (SPI1/IRF8/IKZF1)",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 8.4) +
  theme(
    text = element_text(family = "Arial"),
    panel.grid.major = element_line(color = "#E6E6E6", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8.2, angle = 0, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 6.9, color = "#333333"),
    plot.title = element_text(size = 10.6, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 6.3, hjust = 0.5, color = "#555555"),
    legend.title = element_text(size = 7.6),
    legend.text = element_text(size = 7),
    plot.margin = margin(7, 8, 7, 6)
  )

ggsave(file.path(candidate_dir, "D_四亚型代表regulon_GO功能富集.pdf"), p, width = 7.4, height = 4.7, device = cairo_pdf)
ggsave(file.path(main_dir, "D_四亚型代表regulon_GO功能富集.pdf"), p, width = 7.4, height = 4.7, device = cairo_pdf)
ggsave(file.path(candidate_dir, "D_四亚型代表regulon_GO功能富集.png"), p, width = 7.4, height = 4.7, dpi = 300)

manifest <- tibble(
  panel = "D",
  figure = "D_四亚型代表regulon_GO功能富集.pdf",
  role = "GO BP enrichment for target genes of multiple representative regulons per malignant subtype",
  representative_TFs = paste(paste(rep_tfs$subtype, rep_tfs$TF, sep = "="), collapse = "; "),
  target_definition = paste0("Top ", target_n_per_tf, " SCENIC GRN target genes by edge importance per representative TF; unioned within subtype"),
  claim_boundary = "Functional context for subtype-associated regulon sets; not direct TF-target causality",
  source_data = "source_data/D_四亚型代表regulon_target_GO/GO_BP_enrichment_plot_terms.csv"
)
write_csv(manifest, file.path(record_dir, "PanelD_四亚型代表regulon_GO功能富集_manifest.csv"))
