suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(scales)
  library(tibble)
})

step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir && dir.exists(step_dir)) {
  setwd(step_dir)
}

out_dir <- "figures/supplement"
src_dir <- "figures/source_data/supplement"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(src_dir, showWarnings = FALSE, recursive = TRUE)

subtype_levels <- paste0("Subtype", 1:4)
subtype_colors <- c(
  "Subtype1" = "#0072B5",
  "Subtype2" = "#E18727",
  "Subtype3" = "#20854E",
  "Subtype4" = "#BC3C29"
)

theme_supp <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.35, color = "grey40"),
      axis.ticks = element_line(linewidth = 0.3, color = "grey40"),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold", color = "grey15"),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold")
    )
}

save_pdf <- function(plot, filename, width, height) {
  ggsave(file.path(out_dir, filename), plot, width = width, height = height,
         useDingbats = FALSE)
}

pick_file <- function(pattern) {
  hit <- Sys.glob(pattern)
  if (length(hit) == 0) stop("Missing file pattern: ", pattern)
  hit[[1]]
}

# -------------------------------------------------------------------------
# S2: QC and DoubletFinder summary
# -------------------------------------------------------------------------
qc <- readr::read_csv("../02_scRNA_QC/tables/QC_filtering_summary_by_sample.csv", show_col_types = FALSE)
df <- readr::read_csv("../02_scRNA_QC/tables/DoubletFinder_summary_by_sample.csv", show_col_types = FALSE)
ret <- readr::read_csv("../02_scRNA_QC/tables/celltype_retention_before_after_QC.csv", show_col_types = FALSE)

sample_order <- qc %>%
  filter(qc_status == "Pass_pre_doublet") %>%
  arrange(desc(percent), sample) %>%
  pull(sample)

qc <- qc %>%
  mutate(sample = factor(sample, levels = sample_order),
         qc_status = factor(qc_status, levels = c("Pass_pre_doublet", "Fail_MT20", "Fail_MAD")))
df <- df %>%
  mutate(sample = factor(sample, levels = sample_order),
         doubletfinder_class = factor(doubletfinder_class, levels = c("Singlet", "Doublet", "Skipped_low_cells")))

p_qc <- ggplot(qc, aes(sample, percent, fill = qc_status)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = c(Pass_pre_doublet = "#4D9221", Fail_MT20 = "#C51B7D", Fail_MAD = "#BDBDBD")) +
  scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Cells (%)", title = "QC filtering") +
  theme_supp() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = "bottom")

p_df <- ggplot(df, aes(sample, percent, fill = doubletfinder_class)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = c(Singlet = "#4D9221", Doublet = "#D6604D", Skipped_low_cells = "#BDBDBD")) +
  scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Cells (%)", title = "DoubletFinder calls") +
  theme_supp() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = "bottom")

p_ret <- ret %>%
  arrange(retain_rate) %>%
  mutate(annotation = factor(annotation, levels = annotation)) %>%
  ggplot(aes(retain_rate, annotation)) +
  geom_col(fill = "#4393C3", width = 0.65) +
  scale_x_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.04))) +
  labs(x = "Retention after QC (%)", y = NULL, title = "Cell-type retention") +
  theme_supp() +
  theme(axis.text.y = element_text(size = 6))

readr::write_csv(qc, file.path(src_dir, "S2_QC_filtering_summary_by_sample.csv"))
readr::write_csv(df, file.path(src_dir, "S2_DoubletFinder_summary_by_sample.csv"))
readr::write_csv(ret, file.path(src_dir, "S2_celltype_retention_before_after_QC.csv"))
save_pdf((p_qc / p_df) | p_ret, "S2_QC_DoubletFinder_summary.pdf", 9.5, 6.2)

# -------------------------------------------------------------------------
# S3: inferCNV malignant calling summary
# -------------------------------------------------------------------------
infer_sum <- readr::read_csv(
  pick_file("../04_inferCNV_*/tables/infercnv_immune_reference_call_summary_by_sample.csv"),
  show_col_types = FALSE
)
infer_sum <- infer_sum %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    call_group = case_when(
      infercnv_call == "malignant_like_CNV_high_confidence" ~ "High-confidence malignant",
      infercnv_call == "non_malignant_reference" ~ "Reference",
      grepl("^non_malignant_like", infercnv_call) ~ "Non-malignant-like",
      grepl("^malignant_like", infercnv_call) ~ "Other malignant-like",
      TRUE ~ "Non-malignant-like"
    ),
    call_group = factor(call_group, levels = c("High-confidence malignant", "Other malignant-like", "Non-malignant-like", "Reference"))
  )

p_infer_stack <- infer_sum %>%
  group_by(sample, call_group) %>%
  summarise(percent = sum(percent), n_cells = sum(n_cells), .groups = "drop") %>%
  ggplot(aes(sample, percent, fill = call_group)) +
  geom_col(width = 0.82, color = "white", linewidth = 0.12) +
  scale_fill_manual(values = c(
    "High-confidence malignant" = "#B2182B",
    "Other malignant-like" = "#EF8A62",
    "Non-malignant-like" = "#92C5DE",
    "Reference" = "#D9D9D9"
  )) +
  scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Cells (%)", title = "Immune-reference inferCNV calls") +
  theme_supp() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = "bottom")

p_infer_hc <- infer_sum %>%
  filter(infercnv_call == "malignant_like_CNV_high_confidence") %>%
  ggplot(aes(sample, percent)) +
  geom_col(fill = "#B2182B", width = 0.75) +
  scale_y_continuous(labels = label_percent(scale = 1), expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "High-confidence malignant (%)", title = "High-confidence malignant fraction") +
  theme_supp() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))

readr::write_csv(infer_sum, file.path(src_dir, "S3_inferCNV_immune_reference_call_summary.csv"))
save_pdf(p_infer_stack / p_infer_hc, "S3_inferCNV_malignant_calling_summary.pdf", 8.6, 6.2)

# -------------------------------------------------------------------------
# S4: subtype CNV / doublet audit
# -------------------------------------------------------------------------
cnv <- readr::read_csv("tables/urgent_reaudit_cnv_by_subtype.csv", show_col_types = FALSE) %>%
  mutate(subtype_k4 = factor(subtype_k4, levels = subtype_levels))

p_cnv1 <- ggplot(cnv, aes(subtype_k4, cnv_burden_z_median, color = subtype_k4)) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_errorbar(aes(ymin = cnv_burden_z_q25, ymax = cnv_burden_z_q75), width = 0.18, linewidth = 0.55) +
  geom_point(size = 2.7) +
  scale_color_manual(values = subtype_colors) +
  labs(x = NULL, y = "CNV burden z-score", title = "CNV burden") +
  theme_supp() + theme(legend.position = "none")

p_cnv2 <- ggplot(cnv, aes(subtype_k4, cnv_correlation_ref_z_median, color = subtype_k4)) +
  geom_hline(yintercept = 3, linetype = "dashed", color = "grey45", linewidth = 0.35) +
  geom_errorbar(aes(ymin = cnv_correlation_ref_z_q25, ymax = cnv_correlation_ref_z_q75), width = 0.18, linewidth = 0.55) +
  geom_point(size = 2.7) +
  scale_color_manual(values = subtype_colors) +
  labs(x = NULL, y = "CNV correlation z-score", title = "CNV correlation") +
  theme_supp() + theme(legend.position = "none")

p_cnv3 <- ggplot(cnv, aes(subtype_k4, doublet_rate_pct, fill = subtype_k4)) +
  geom_col(width = 0.65) +
  scale_fill_manual(values = subtype_colors) +
  labs(x = NULL, y = "Doublet rate (%)", title = "Residual doublet rate") +
  theme_supp() + theme(legend.position = "none")

readr::write_csv(cnv, file.path(src_dir, "S4_subtype_CNV_doublet_audit.csv"))
save_pdf(p_cnv1 + p_cnv2 + p_cnv3, "S4_subtype_CNV_doublet_audit.pdf", 8.5, 3.0)

# -------------------------------------------------------------------------
# S5: non-malignant signature scoring
# -------------------------------------------------------------------------
sig_cell <- readr::read_tsv("tables/non_malignant_signature_per_cell_scores.tsv", show_col_types = FALSE) %>%
  mutate(subtype_k4 = factor(subtype_k4, levels = subtype_levels)) %>%
  pivot_longer(c(Pericyte, vSMC, Microglia, TAM, Endothelial), names_to = "signature", values_to = "score") %>%
  mutate(signature = factor(signature, levels = c("Pericyte", "vSMC", "Microglia", "TAM", "Endothelial")))

set.seed(42)
sig_plot <- sig_cell %>%
  group_by(subtype_k4, signature) %>%
  group_modify(~ dplyr::slice_sample(.x, n = min(2500, nrow(.x)))) %>%
  ungroup()

p_sig <- ggplot(sig_plot, aes(subtype_k4, score, fill = subtype_k4)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, color = "grey35") +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.25, fill = "white", alpha = 0.75) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey35", linewidth = 0.25) +
  facet_wrap(~signature, nrow = 1) +
  scale_fill_manual(values = subtype_colors) +
  labs(x = NULL, y = "UCell score") +
  theme_supp() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

readr::write_csv(sig_plot, file.path(src_dir, "S5_non_malignant_signature_plot_downsampled.csv"))
save_pdf(p_sig, "S5_non_malignant_signature_violin.pdf", 9.2, 3.0)

# -------------------------------------------------------------------------
# S6: per-patient NMF cophenetic curves and chosen k
# -------------------------------------------------------------------------
cons <- readr::read_csv("tables/10a_per_patient_consensus_metrics.csv", show_col_types = FALSE)
chosen <- readr::read_csv("tables/10a_per_patient_chosen_k.csv", show_col_types = FALSE)
cons <- cons %>% left_join(chosen %>% select(patient, chosen_k), by = "patient")

cons_plot <- cons %>%
  mutate(
    k = as.integer(k),
    chosen = k == chosen_k,
    highlight = case_when(
      patient == "Pt9" ~ "Pt9",
      chosen_k != 4 ~ "non-k4",
      TRUE ~ "other"
    )
  )

cophenetic_summary <- cons_plot %>%
  group_by(k) %>%
  summarise(
    median_cophenetic = median(cophenetic, na.rm = TRUE),
    q25 = quantile(cophenetic, 0.25, na.rm = TRUE),
    q75 = quantile(cophenetic, 0.75, na.rm = TRUE),
    min_cophenetic = min(cophenetic, na.rm = TRUE),
    max_cophenetic = max(cophenetic, na.rm = TRUE),
    .groups = "drop"
  )

label_points <- cons_plot %>%
  filter(chosen & chosen_k != 4) %>%
  mutate(label = paste0(patient, " k=", chosen_k))

p_coph <- ggplot() +
  geom_ribbon(
    data = cophenetic_summary,
    aes(k, ymin = q25, ymax = q75),
    fill = "grey86", alpha = 0.8
  ) +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.35, color = "#B2182B") +
  geom_line(
    data = cons_plot,
    aes(k, cophenetic, group = patient),
    color = "grey62", alpha = 0.38, linewidth = 0.32
  ) +
  geom_point(
    data = cons_plot,
    aes(k, cophenetic),
    color = "grey55", alpha = 0.55, size = 1.15
  ) +
  geom_line(
    data = cophenetic_summary,
    aes(k, median_cophenetic),
    color = "grey10", linewidth = 0.9
  ) +
  geom_point(
    data = cophenetic_summary,
    aes(k, median_cophenetic),
    shape = 21, fill = "white", color = "grey10", stroke = 0.45, size = 2.2
  ) +
  geom_point(
    data = label_points,
    aes(k, cophenetic),
    shape = 21, fill = "#B2182B", color = "white", stroke = 0.35, size = 2.4
  ) +
  geom_text(
    data = label_points,
    aes(k, cophenetic, label = label),
    nudge_x = 0.10, nudge_y = 0.004,
    size = 2.35, hjust = 0, color = "#B2182B", fontface = "bold"
  ) +
  scale_x_continuous(breaks = 4:7, limits = c(3.9, 7.45), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(limits = c(0.90, 1.005), breaks = seq(0.90, 1.00, 0.025)) +
  labs(x = "NMF k", y = "Cophenetic correlation", title = "Per-patient NMF stability") +
  theme_supp() +
  theme(
    panel.grid.major.y = element_line(color = "grey91", linewidth = 0.25),
    axis.text = element_text(size = 7.2),
    legend.position = "none"
  )

p_chosen <- chosen %>%
  count(chosen_k) %>%
  ggplot(aes(factor(chosen_k), n)) +
  geom_col(fill = "#4D4D4D", width = 0.56) +
  geom_text(aes(label = n), vjust = -0.35, size = 2.6, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)), breaks = scales::pretty_breaks(n = 4)) +
  labs(x = "Chosen k", y = "Patients", title = "Chosen k distribution") +
  theme_supp() +
  theme(axis.text = element_text(size = 7.2))

readr::write_csv(cons, file.path(src_dir, "S6_NMF_cophenetic_metrics.csv"))
readr::write_csv(cons_plot, file.path(src_dir, "S6_NMF_cophenetic_curve_plot_data.csv"))
readr::write_csv(cophenetic_summary, file.path(src_dir, "S6_NMF_cophenetic_summary.csv"))
save_pdf(p_coph | p_chosen + plot_layout(widths = c(3.1, 1)), "S6_NMF_cophenetic_chosen_k.pdf", 7.2, 3.4)

# -------------------------------------------------------------------------
# S7: metaprogram clustering diagnostics
# -------------------------------------------------------------------------
clsel <- readr::read_csv("tables/10b_metaprogram_cluster_selection.csv", show_col_types = FALSE)
mpsum <- readr::read_csv("tables/10b_metaprogram_summary.csv", show_col_types = FALSE)

p_mp_k <- ggplot(clsel, aes(candidate_k, mean_silhouette)) +
  geom_line(color = "grey35", linewidth = 0.45) +
  geom_point(aes(fill = selected), shape = 21, size = 2.4, color = "grey20") +
  scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "white")) +
  scale_x_continuous(breaks = clsel$candidate_k) +
  labs(x = "Candidate metaprogram clusters", y = "Mean silhouette", title = "Metaprogram cluster-count selection") +
  theme_supp() + theme(legend.position = "none")

p_mp_rec <- mpsum %>%
  mutate(metaprogram_id = factor(metaprogram_id, levels = metaprogram_id)) %>%
  ggplot(aes(metaprogram_id, n_patients_unique, fill = is_recurrent)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 6, linetype = "dashed", linewidth = 0.35, color = "grey35") +
  scale_fill_manual(values = c(`TRUE` = "#4D9221", `FALSE` = "#D9D9D9")) +
  labs(x = NULL, y = "Unique patients", title = "Metaprogram recurrence") +
  theme_supp() + theme(legend.position = "none")

readr::write_csv(clsel, file.path(src_dir, "S7_metaprogram_cluster_selection.csv"))
readr::write_csv(mpsum, file.path(src_dir, "S7_metaprogram_summary.csv"))
save_pdf(p_mp_k | p_mp_rec, "S7_NMF_metaprogram_diagnostics.pdf", 8.2, 3.2)

# -------------------------------------------------------------------------
# S9: S3 vs S4 marker DE volcano
# -------------------------------------------------------------------------
de <- readr::read_csv("tables/10d_S3_vs_S4_FindMarkers_all.csv", show_col_types = FALSE) %>%
  mutate(
    sig = BH_q < 0.05 & abs(avg_log2FC) > 0.5,
    direction = case_when(
      sig & avg_log2FC > 0 ~ "Subtype3-up",
      sig & avg_log2FC < 0 ~ "Subtype4-up",
      TRUE ~ "Not passing"
    ),
    neg_log10_q = pmin(-log10(BH_q + 1e-300), 80),
    direction = factor(direction, levels = c("Not passing", "Subtype4-up", "Subtype3-up"))
  )

top_label <- de %>%
  filter(sig) %>%
  group_by(direction) %>%
  arrange(BH_q, desc(abs(avg_log2FC)), .by_group = TRUE) %>%
  slice_head(n = 7) %>%
  ungroup()

volcano_counts <- de %>%
  count(direction) %>%
  mutate(label = paste0(direction, ": ", n))

p_volcano <- ggplot(de, aes(avg_log2FC, neg_log10_q)) +
  annotate("rect", xmin = -0.5, xmax = 0.5, ymin = -Inf, ymax = Inf,
           fill = "grey96", color = NA) +
  geom_point(
    data = de %>% filter(direction == "Not passing"),
    color = "grey82", size = 0.28, alpha = 0.42, stroke = 0
  ) +
  geom_point(
    data = de %>% filter(direction != "Not passing"),
    aes(color = direction),
    size = 0.62, alpha = 0.72, stroke = 0
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.32, color = "grey38") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.32, color = "grey38") +
  scale_color_manual(
    values = c("Subtype3-up" = "#20854E", "Subtype4-up" = "#BC3C29", "Not passing" = "grey78"),
    breaks = c("Subtype4-up", "Subtype3-up")
  ) +
  ggrepel::geom_text_repel(
    data = top_label,
    aes(label = gene),
    size = 2.45,
    max.overlaps = Inf,
    box.padding = 0.28,
    point.padding = 0.12,
    min.segment.length = 0,
    segment.size = 0.22,
    segment.color = "grey45"
  ) +
  annotate("text", x = min(de$avg_log2FC, na.rm = TRUE) * 0.92, y = 77,
           label = "Subtype4-up", hjust = 0, size = 3.0, fontface = "bold", color = "#BC3C29") +
  annotate("text", x = max(de$avg_log2FC, na.rm = TRUE) * 0.92, y = 77,
           label = "Subtype3-up", hjust = 1, size = 3.0, fontface = "bold", color = "#20854E") +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.04))) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.04))) +
  labs(
    x = expression(log[2]~"fold change (Subtype3 / Subtype4)"),
    y = expression(-log[10]~"BH q"),
    color = NULL,
    title = "Subtype3 vs Subtype4 marker genes"
  ) +
  theme_supp(base_size = 8.5) +
  theme(
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.25),
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 10.5, face = "bold"),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )

readr::write_csv(top_label, file.path(src_dir, "S9_S3_vs_S4_DE_top_labeled_genes.csv"))
readr::write_csv(volcano_counts, file.path(src_dir, "S9_S3_vs_S4_DE_direction_counts.csv"))
save_pdf(p_volcano, "S9_S3_vs_S4_marker_DE_volcano.pdf", 5.8, 4.8)

# -------------------------------------------------------------------------
# S10: S3 vs S4 GSEA diverging plot
# -------------------------------------------------------------------------
gsea34 <- readr::read_csv("tables/10e_S3_vs_S4_GSEA_combined.csv", show_col_types = FALSE) %>%
  filter(p.adjust < 0.05) %>%
  mutate(
    direction_clean = ifelse(direction == "Subtype3_enriched", "Subtype3", "Subtype4"),
    NES_signed = ifelse(direction_clean == "Subtype3", abs(NES), -abs(NES)),
    term_label = Description %>%
      str_remove("^HALLMARK_") %>%
      str_remove("^REACTOME_") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence(),
    neg_log10_q = pmin(-log10(p.adjust + 1e-300), 10)
  )

top_gsea34 <- bind_rows(
  gsea34 %>% filter(direction_clean == "Subtype3") %>% arrange(p.adjust, desc(abs(NES))) %>% slice_head(n = 10),
  gsea34 %>% filter(direction_clean == "Subtype4") %>% arrange(p.adjust, desc(abs(NES))) %>% slice_head(n = 10)
) %>%
  mutate(
    term_label_wrapped = str_wrap(term_label, width = 34)
  )

gsea_label_levels <- top_gsea34 %>%
  arrange(NES_signed) %>%
  pull(term_label_wrapped) %>%
  rev()

top_gsea34 <- top_gsea34 %>%
  mutate(
    term_label_wrapped = factor(term_label_wrapped, levels = unique(gsea_label_levels))
  )

stopifnot(!any(is.na(top_gsea34$term_label_wrapped)))

p_gsea34 <- ggplot(top_gsea34, aes(NES_signed, term_label_wrapped, fill = direction_clean)) +
  geom_vline(xintercept = 0, color = "grey35", linewidth = 0.38) +
  geom_col(aes(alpha = neg_log10_q), width = 0.62, color = "white", linewidth = 0.18) +
  geom_text(
    aes(label = sprintf("%.2f", abs(NES)),
        hjust = ifelse(NES_signed > 0, -0.15, 1.15)),
    size = 2.25, color = "grey20"
  ) +
  scale_fill_manual(values = c(Subtype3 = "#20854E", Subtype4 = "#BC3C29")) +
  scale_alpha_continuous(range = c(0.72, 1), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  labs(x = "Signed normalized enrichment score (NES)", y = NULL, fill = NULL,
       title = "Subtype3 vs Subtype4 enriched pathways") +
  theme_supp(base_size = 8.5) +
  theme(
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
    axis.text.y = element_text(size = 7.2, lineheight = 0.88),
    axis.text.x = element_text(size = 8),
    axis.title.x = element_text(size = 9),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 10.5, face = "bold")
  ) +
  coord_cartesian(clip = "off")

readr::write_csv(top_gsea34, file.path(src_dir, "S10_S3_vs_S4_GSEA_top_terms.csv"))
save_pdf(p_gsea34, "S10_S3_vs_S4_GSEA_diverging.pdf", 7.2, 5.2)

# -------------------------------------------------------------------------
# S11: original manuscript marker remapping
# -------------------------------------------------------------------------
remap <- readr::read_csv("tables/manuscript_subtype_remapping.csv", show_col_types = FALSE)
remap_long <- remap %>%
  select(original_subtype, original_marker, starts_with("mean_expr_S"), starts_with("pct_expressed_S")) %>%
  pivot_longer(starts_with("mean_expr_S"), names_to = "subtype", values_to = "mean_expr") %>%
  mutate(subtype_k4 = paste0("Subtype", str_extract(subtype, "[1-4]"))) %>%
  select(-subtype) %>%
  left_join(
    remap %>%
      select(original_marker, starts_with("pct_expressed_S")) %>%
      pivot_longer(starts_with("pct_expressed_S"), names_to = "subtype", values_to = "pct_expr") %>%
      mutate(subtype_k4 = paste0("Subtype", str_extract(subtype, "[1-4]"))) %>%
      select(original_marker, subtype_k4, pct_expr),
    by = c("original_marker", "subtype_k4")
  ) %>%
  mutate(
    subtype_k4 = factor(subtype_k4, levels = subtype_levels),
    marker_label = paste0(original_marker, "\n(", original_subtype, ")"),
    marker_label = factor(marker_label, levels = unique(marker_label))
  )

p_remap <- ggplot(remap_long, aes(subtype_k4, marker_label)) +
  geom_point(aes(size = pct_expr, fill = mean_expr), shape = 21, color = "grey25", stroke = 0.3) +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "Mean expr.") +
  scale_size_continuous(range = c(1.8, 7), name = "% expressed") +
  labs(x = NULL, y = NULL, title = "Original manuscript marker remapping") +
  theme_supp() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

readr::write_csv(remap_long, file.path(src_dir, "S11_original_marker_remapping_plot_data.csv"))
save_pdf(p_remap, "S11_original_marker_remapping.pdf", 5.8, 3.4)

writeLines(capture.output(sessionInfo()), file.path(src_dir, "S2_S11_supplement_batch_session_info.txt"))

cat("Supplement batch completed.\n")
cat("Output directory:", out_dir, "\n")
print(list.files(out_dir, pattern = "^S[0-9].*\\.pdf$", full.names = FALSE))
