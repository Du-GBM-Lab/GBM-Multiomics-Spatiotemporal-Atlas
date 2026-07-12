suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(UCell)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
})

# URGENT continuation audit:
# 1) Score non-malignant lineage signatures in the final malignant object.
# 2) Cross-reference revised subtype top50 markers with canonical GBM literature anchors.
#
# This script does not change the Seurat object on disk and does not rerun DE/NMF.

step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir) {
  setwd(file.path(getwd(), step_dir))
}

dir.create(file.path("figures", "audit"), recursive = TRUE, showWarnings = FALSE)

obj_path <- file.path("outputs", "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2")
de_path <- file.path("tables", "10e_one_vs_rest_FindAllMarkers_wilcoxon_downsampled.csv")

out_sig_used <- file.path("tables", "non_malignant_signature_genes_used.csv")
out_per_cell <- file.path("tables", "non_malignant_signature_per_cell_scores.tsv")
out_summary <- file.path("tables", "non_malignant_signature_per_subtype.csv")
out_comp <- file.path("tables", "signature_score_comparison.csv")
out_pdf <- file.path("figures", "audit", "non_malignant_signature_violin.pdf")

out_anchor_detail <- file.path("tables", "canonical_marker_anchoring_detail.csv")
out_anchor_summary <- file.path("tables", "canonical_marker_anchoring_summary.csv")
out_anchor_sanity <- file.path("tables", "canonical_marker_anchoring_sanity_checks.csv")
out_session <- file.path("tables", "urgent_nonmalignant_signature_and_marker_anchor_session_info.txt")

subtype_levels <- paste0("Subtype", 1:4)
subtype_cols <- c(
  "Subtype1" = "#0072B5",
  "Subtype2" = "#E18727",
  "Subtype3" = "#20854E",
  "Subtype4" = "#BC3C29"
)

signature_list <- list(
  Pericyte = c("RGS5", "ACTA2", "TAGLN", "MYH11", "PDGFRB", "KCNJ8", "ABCC9", "HIGD1B", "ADIRF", "NDUFA4L2"),
  vSMC = c("ACTA2", "TAGLN", "MYH11", "ACTG2", "CNN1", "LMOD1"),
  Microglia = c("P2RY12", "P2RY13", "CX3CR1", "TMEM119", "SLC2A5", "CSF1R", "AIF1", "ITGAM"),
  TAM = c("CD163", "MRC1", "MSR1", "CD68", "CD84", "MS4A7", "FCGR1B", "FCGR3A"),
  Endothelial = c("PECAM1", "CDH5", "VWF", "CLDN5", "KDR", "FLT1")
)

obj <- qs2::qs_read(obj_path)
DefaultAssay(obj) <- "RNA"
md <- obj@meta.data
stopifnot("subtype_k4" %in% colnames(md))
md$subtype_k4 <- factor(as.character(md$subtype_k4), levels = subtype_levels)
stopifnot(!any(is.na(md$subtype_k4)))

genes_available <- rownames(obj)
signature_used <- lapply(names(signature_list), function(sig) {
  genes <- signature_list[[sig]]
  tibble(
    signature = sig,
    gene = genes,
    used_in_object = genes %in% genes_available
  )
}) %>% bind_rows()
readr::write_csv(signature_used, out_sig_used)

sig_for_ucell <- lapply(signature_list, function(genes) intersect(genes, genes_available))
if (any(lengths(sig_for_ucell) < 3)) {
  stop("At least one signature has fewer than 3 detected genes: ",
       paste(names(sig_for_ucell)[lengths(sig_for_ucell) < 3], collapse = ", "))
}

obj_scored <- UCell::AddModuleScore_UCell(
  obj,
  features = sig_for_ucell,
  maxRank = 1500,
  ncores = 1
)
score_cols <- paste0(names(sig_for_ucell), "_UCell")
missing_scores <- setdiff(score_cols, colnames(obj_scored@meta.data))
if (length(missing_scores) > 0) {
  stop("Missing UCell score columns: ", paste(missing_scores, collapse = ", "))
}

score_md <- obj_scored@meta.data %>%
  as_tibble(rownames = "cell_id") %>%
  select(cell_id, Pt_number, subtype_k4, all_of(score_cols)) %>%
  rename_with(~ sub("_UCell$", "", .x), all_of(score_cols))

readr::write_tsv(score_md, out_per_cell)

score_long <- score_md %>%
  pivot_longer(cols = all_of(names(sig_for_ucell)), names_to = "signature", values_to = "UCell_score")

sig_summary <- score_long %>%
  group_by(subtype_k4, signature) %>%
  summarise(
    n_cells = n(),
    mean = mean(UCell_score, na.rm = TRUE),
    sd = sd(UCell_score, na.rm = TRUE),
    median = median(UCell_score, na.rm = TRUE),
    q25 = quantile(UCell_score, 0.25, na.rm = TRUE),
    q75 = quantile(UCell_score, 0.75, na.rm = TRUE),
    max = max(UCell_score, na.rm = TRUE),
    pct_gt_0_5 = mean(UCell_score > 0.5, na.rm = TRUE) * 100,
    .groups = "drop"
  )
readr::write_csv(sig_summary, out_summary)

cliff_delta_fast <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  n1 <- length(x)
  n2 <- length(y)
  if (n1 == 0 || n2 == 0) return(NA_real_)
  ranks <- rank(c(x, y), ties.method = "average")
  r1 <- sum(ranks[seq_len(n1)])
  u1 <- r1 - n1 * (n1 + 1) / 2
  (2 * u1 / (n1 * n2)) - 1
}

comparison_plan <- tibble::tribble(
  ~signature, ~target_subtype, ~comparison_group,
  "Pericyte", "Subtype3", "Subtype1+Subtype2",
  "vSMC", "Subtype3", "Subtype1+Subtype2",
  "Endothelial", "Subtype3", "Subtype1+Subtype2",
  "Microglia", "Subtype4", "Subtype1+Subtype2",
  "TAM", "Subtype4", "Subtype1+Subtype2",
  "Endothelial", "Subtype4", "Subtype1+Subtype2"
)

comparison <- comparison_plan %>%
  rowwise() %>%
  mutate(
    n_target = sum(score_md$subtype_k4 == target_subtype),
    n_comparator = sum(score_md$subtype_k4 %in% c("Subtype1", "Subtype2")),
    mean_target = mean(score_md[[signature]][score_md$subtype_k4 == target_subtype], na.rm = TRUE),
    mean_comparator = mean(score_md[[signature]][score_md$subtype_k4 %in% c("Subtype1", "Subtype2")], na.rm = TRUE),
    median_target = median(score_md[[signature]][score_md$subtype_k4 == target_subtype], na.rm = TRUE),
    median_comparator = median(score_md[[signature]][score_md$subtype_k4 %in% c("Subtype1", "Subtype2")], na.rm = TRUE),
    wilcox_p = wilcox.test(
      score_md[[signature]][score_md$subtype_k4 == target_subtype],
      score_md[[signature]][score_md$subtype_k4 %in% c("Subtype1", "Subtype2")]
    )$p.value,
    cliff_delta = cliff_delta_fast(
      score_md[[signature]][score_md$subtype_k4 == target_subtype],
      score_md[[signature]][score_md$subtype_k4 %in% c("Subtype1", "Subtype2")]
    )
  ) %>%
  ungroup() %>%
  mutate(BH_q = p.adjust(wilcox_p, method = "BH")) %>%
  select(signature, target_subtype, comparison_group, n_target, n_comparator,
         mean_target, mean_comparator, median_target, median_comparator,
         wilcox_p, BH_q, cliff_delta)
readr::write_csv(comparison, out_comp)

p <- ggplot(score_long, aes(x = subtype_k4, y = UCell_score, fill = subtype_k4)) +
  geom_violin(scale = "width", linewidth = 0.2, trim = TRUE, color = "grey30") +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.2, fill = "white", alpha = 0.65) +
  facet_wrap(~ signature, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = subtype_cols, guide = "none") +
  labs(x = NULL, y = "UCell score") +
  theme_classic(base_size = 8) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face = "bold"),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 8),
    panel.spacing = unit(4, "mm"),
    plot.margin = margin(4, 4, 4, 4)
  )
ggsave(out_pdf, p, width = 8.8, height = 2.8, useDingbats = FALSE)

de <- readr::read_csv(de_path, show_col_types = FALSE)
stopifnot(all(c("subtype_k4", "gene", "avg_log2FC", "BH_q") %in% colnames(de)))
top50 <- de %>%
  filter(avg_log2FC > 0, BH_q < 0.05) %>%
  group_by(subtype_k4) %>%
  arrange(desc(avg_log2FC), BH_q, .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  slice_head(n = 50) %>%
  ungroup()

anchor_sets <- list(
  Subtype1 = list(
    Proliferation = c("MKI67", "TOP2A", "CDK1", "AURKB", "UBE2C", "CCNB1", "CCNB2", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "CENPF", "BIRC5", "STMN1"),
    NPC = c("STMN2", "ASCL1", "DLL3", "SOX4", "NHLH1", "DLX5", "DLX6", "NEUROD4", "CD24", "LMO3")
  ),
  Subtype2 = list(
    OPC_Myelination = c("MAG", "MOG", "MBP", "PLP1", "OLIG1", "OLIG2", "SOX10", "PDGFRA", "NKX2-2", "NKX2.2", "MOBP", "ERMN", "CNDP1", "CARNS1", "OPALIN", "NKX6-2")
  ),
  Subtype3 = list(
    Classical_MES_ECM = c("SERPINE1", "TIMP1", "CHI3L1", "MGP", "SPP1", "VIM", "FN1", "LOX", "LGALS1", "LGALS3", "A2M", "CD44", "ANXA1", "ANXA2", "COL1A1", "COL1A2", "NDRG1", "BNIP3", "ADM", "HILPDA", "IGFBP3", "PLAU", "PLAUR", "NAMPT", "S100A11", "EMP3"),
    Niche_Pericyte_Coupled = c("RGS5", "ACTA2", "TAGLN", "HIGD1B", "ADIRF", "BMP4", "FRZB", "LRRC32", "PDGFRB", "KCNJ8", "ABCC9")
  ),
  Subtype4 = list(
    Antigen_Presentation = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "B2M", "CIITA", "IFI30", "CTSS", "PSMB9", "HLA-DQA1", "HLA-DQB1", "HLA-DRB5"),
    Malignant_Myeloid_Like = c("SPP1", "C1QA", "C1QB", "C1QC", "TYROBP", "CD14", "FCGR3A", "CCL3", "CCL4", "CCL3L3", "CCL4L2", "MS4A7", "LAPTM5", "CD68", "AIF1", "APOE", "SRGN")
  )
)

source_lookup <- tibble::tribble(
  ~subtype_k4, ~anchor_class, ~source_paper,
  "Subtype1", "Proliferation", "Neftel 2019 Cell; cycling/proliferative GBM programs",
  "Subtype1", "NPC", "Neftel 2019 Cell NPC-like malignant state",
  "Subtype2", "OPC_Myelination", "Neftel 2019 Cell OPC-like state; oligodendrocyte/myelination GBM programs",
  "Subtype3", "Classical_MES_ECM", "Verhaak 2010 Cancer Cell; Wang 2017 Cancer Cell; Neftel 2019 Cell MES-like state",
  "Subtype3", "Niche_Pericyte_Coupled", "Greenwald 2024 Cell vascular-coupled malignant niche programs",
  "Subtype4", "Antigen_Presentation", "Schaettler 2022 Cancer Discovery malignant MHC-II / antigen-presentation program",
  "Subtype4", "Malignant_Myeloid_Like", "Hara 2021 Cancer Cell glioma-immune interaction induced MES/myeloid-like state"
)

anchor_detail <- top50 %>%
  rowwise() %>%
  mutate(
    anchor_class = {
      sets <- anchor_sets[[subtype_k4]]
      hits <- names(sets)[vapply(sets, function(v) gene %in% v, logical(1))]
      if (length(hits) == 0) NA_character_ else paste(hits, collapse = ";")
    },
    literature_anchor = !is.na(anchor_class)
  ) %>%
  ungroup() %>%
  left_join(source_lookup, by = c("subtype_k4", "anchor_class")) %>%
  mutate(source_paper = ifelse(is.na(source_paper) & !is.na(anchor_class) & str_detect(anchor_class, ";"), "Multiple anchor classes", source_paper)) %>%
  select(subtype_k4, rank, gene, avg_log2FC, BH_q, literature_anchor, anchor_class, source_paper, pct.1, pct.2)

readr::write_csv(anchor_detail, out_anchor_detail)

anchor_summary <- anchor_detail %>%
  group_by(subtype_k4) %>%
  summarise(
    total_markers_top50 = n(),
    n_anchored = sum(literature_anchor),
    pct_anchored = n_anchored / total_markers_top50 * 100,
    anchor_source_breakdown = paste(
      names(table(anchor_class[literature_anchor])),
      as.integer(table(anchor_class[literature_anchor])),
      sep = "=",
      collapse = ";"
    ),
    anchored_genes = paste(gene[literature_anchor], collapse = ";"),
    .groups = "drop"
  )
readr::write_csv(anchor_summary, out_anchor_summary)

anchor_sanity <- tibble::tibble(
  metric = c(
    "n_cells",
    "n_signatures",
    "min_signature_genes_detected",
    "n_top50_marker_rows",
    "n_subtypes_with_50_markers",
    "Subtype3_classical_MES_ECM_hits",
    "Subtype3_niche_pericyte_hits",
    "Subtype4_antigen_presentation_hits",
    "Subtype4_malignant_myeloid_like_hits"
  ),
  value = c(
    nrow(score_md),
    length(sig_for_ucell),
    min(lengths(sig_for_ucell)),
    nrow(anchor_detail),
    length(unique(anchor_detail$subtype_k4)),
    sum(anchor_detail$subtype_k4 == "Subtype3" & anchor_detail$anchor_class == "Classical_MES_ECM"),
    sum(anchor_detail$subtype_k4 == "Subtype3" & anchor_detail$anchor_class == "Niche_Pericyte_Coupled"),
    sum(anchor_detail$subtype_k4 == "Subtype4" & anchor_detail$anchor_class == "Antigen_Presentation"),
    sum(anchor_detail$subtype_k4 == "Subtype4" & anchor_detail$anchor_class == "Malignant_Myeloid_Like")
  )
)
readr::write_csv(anchor_sanity, out_anchor_sanity)

writeLines(capture.output(sessionInfo()), out_session)

cat("Non-malignant signature per subtype:\n")
print(sig_summary)
cat("\nSignature comparisons:\n")
print(comparison)
cat("\nCanonical marker anchoring summary:\n")
print(anchor_summary)
cat("\nOutputs:\n")
cat(out_summary, "\n", out_comp, "\n", out_pdf, "\n", out_anchor_detail, "\n", out_anchor_summary, "\n", sep = "")
