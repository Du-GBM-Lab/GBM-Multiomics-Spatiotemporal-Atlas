suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir && dir.exists(step_dir)) {
  setwd(step_dir)
}

out_dir <- "figures/supplement"
src_dir <- "figures/source_data/supplement"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(src_dir, showWarnings = FALSE, recursive = TRUE)

subtype_colors <- c(
  "Subtype1" = "#0072B5",
  "Subtype2" = "#E18727",
  "Subtype3" = "#20854E",
  "Subtype4" = "#BC3C29"
)

mp_colors <- c(
  "MP01" = "#6A51A3",
  "MP02" = "#20854E",
  "MP03" = "#E18727",
  "MP04" = "#BC3C29",
  "MP05" = "#7F7F7F",
  "MP06" = "#0072B5",
  "MP07" = "#BDBDBD"
)

theme_supp <- function(base_size = 8.5) {
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

# -------------------------------------------------------------------------
# S8: S3 vs S4 within-patient metaprogram consistency
# -------------------------------------------------------------------------
scores <- readr::read_tsv("tables/10c_per_cell_metaprogram_scores.tsv", show_col_types = FALSE)
mp_test <- readr::read_csv("tables/10c_S3_vs_S4_metaprogram_test.csv", show_col_types = FALSE)

mps <- c("MP01", "MP02", "MP03", "MP04", "MP06")
mp_order <- c("MP02", "MP06", "MP04", "MP03", "MP01")

patient_delta <- scores %>%
  filter(subtype_k4 %in% c("Subtype3", "Subtype4")) %>%
  select(Pt_number, subtype_k4, all_of(mps)) %>%
  pivot_longer(all_of(mps), names_to = "metaprogram", values_to = "score") %>%
  group_by(Pt_number, subtype_k4, metaprogram) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), n_cells = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = subtype_k4,
    values_from = c(mean_score, n_cells),
    names_sep = "__"
  ) %>%
  filter(!is.na(mean_score__Subtype3), !is.na(mean_score__Subtype4)) %>%
  mutate(
    delta_S3_minus_S4 = mean_score__Subtype3 - mean_score__Subtype4,
    sign_delta = sign(delta_S3_minus_S4)
  ) %>%
  left_join(
    mp_test %>%
      transmute(
        metaprogram,
        global_delta = delta_mean_S3_minus_S4,
        global_sign = sign(delta_mean_S3_minus_S4),
        within_patient_consistency_pct
      ),
    by = "metaprogram"
  ) %>%
  mutate(
    consistent_with_global = sign_delta == global_sign,
    metaprogram = factor(metaprogram, levels = mp_order)
  )

patient_order <- patient_delta %>%
  filter(metaprogram %in% c("MP02", "MP04")) %>%
  select(Pt_number, metaprogram, delta_S3_minus_S4) %>%
  pivot_wider(names_from = metaprogram, values_from = delta_S3_minus_S4) %>%
  mutate(order_score = MP02 - MP04) %>%
  arrange(desc(order_score), Pt_number) %>%
  pull(Pt_number)

patient_delta <- patient_delta %>%
  mutate(Pt_number = factor(Pt_number, levels = rev(patient_order)))

mp_label <- mp_test %>%
  filter(metaprogram %in% mps) %>%
  mutate(
    metaprogram = factor(metaprogram, levels = mp_order),
    label = sprintf("%s\n%.0f%%", metaprogram, within_patient_consistency_pct)
  ) %>%
  arrange(metaprogram) %>%
  select(metaprogram, label)

mp_label_vec <- setNames(as.character(mp_label$label), as.character(mp_label$metaprogram))

p_s8 <- ggplot(patient_delta, aes(metaprogram, Pt_number, fill = delta_S3_minus_S4)) +
  geom_tile(color = "white", linewidth = 0.28) +
  geom_point(
    data = patient_delta %>% filter(!consistent_with_global),
    shape = 4, size = 1.1, stroke = 0.35, color = "grey15"
  ) +
  scale_x_discrete(labels = mp_label_vec, position = "top") +
  scale_fill_gradient2(
    low = "#BC3C29",
    mid = "white",
    high = "#20854E",
    midpoint = 0,
    limits = c(-0.18, 0.18),
    oob = scales::squish,
    name = "Mean score\nS3 - S4"
  ) +
  labs(
    x = NULL,
    y = "Patient",
    title = "Within-patient metaprogram score differences"
  ) +
  theme_supp(base_size = 8.5) +
  theme(
    axis.text.x = element_text(size = 7.4, face = "bold", lineheight = 0.9),
    axis.text.y = element_text(size = 6.5),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey78", linewidth = 0.35),
    legend.position = "right",
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7)
  )

readr::write_csv(patient_delta, file.path(src_dir, "S8_S3_S4_within_patient_metaprogram_delta.csv"))
readr::write_csv(mp_label, file.path(src_dir, "S8_S3_S4_within_patient_consistency_summary.csv"))
ggsave(file.path(out_dir, "S8_S3_S4_within_patient_metaprogram_consistency.pdf"),
       p_s8, width = 5.8, height = 5.2, useDingbats = FALSE)

# -------------------------------------------------------------------------
# S13: NMF program Jaccard clustering heatmap with dendrogram
# -------------------------------------------------------------------------
jaccard <- readRDS("outputs/nmf/10b_jaccard_matrix.rds")
hc <- readRDS("outputs/nmf/10b_program_dendrogram.rds")
assignment <- readr::read_csv("tables/10b_metaprogram_assignment.csv", show_col_types = FALSE)
mp_summary <- readr::read_csv("tables/10b_metaprogram_summary.csv", show_col_types = FALSE)

stopifnot(is.matrix(jaccard))
stopifnot(setequal(rownames(jaccard), colnames(jaccard)))
stopifnot(setequal(rownames(jaccard), assignment$program_id))

assignment <- assignment %>%
  left_join(
    mp_summary %>% select(metaprogram_id, is_recurrent, is_cycling, n_patients_unique),
    by = "metaprogram_id"
  ) %>%
  mutate(
    metaprogram_id = factor(metaprogram_id, levels = names(mp_colors))
  )

ann_df <- assignment %>%
  select(program_id, metaprogram_id) %>%
  tibble::column_to_rownames("program_id") %>%
  .[rownames(jaccard), , drop = FALSE]

mp_label_map <- c(
  "MP01" = "MP01: IGFBP / stress",
  "MP02" = "MP02: ECM organization",
  "MP03" = "MP03: Myelination",
  "MP04" = "MP04: MHC-II antigen presentation",
  "MP05" = "MP05: Cell cycle",
  "MP06" = "MP06: Angiogenesis signaling",
  "MP07" = "MP07: Minor program"
)

top_ann <- HeatmapAnnotation(
  Metaprogram = ann_df$metaprogram_id,
  col = list(Metaprogram = mp_colors),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 7, fontface = "bold"),
  simple_anno_size = unit(2.8, "mm"),
  annotation_legend_param = list(
    Metaprogram = list(
      title = "Metaprogram",
      at = names(mp_colors),
      labels = unname(mp_label_map[names(mp_colors)]),
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      grid_height = unit(3, "mm"),
      grid_width = unit(3, "mm")
    )
  )
)

# Off-diagonal Jaccard similarities are mostly below 0.35; using the full 0-1
# range makes biologically relevant blocks look washed out.  The source matrix
# remains raw, while the display scale is capped for contrast.
jaccard_col <- circlize::colorRamp2(
  c(0, 0.05, 0.10, 0.20, 0.35),
  c("#FFFFFF", "#E0F0F7", "#A6D5E8", "#4A9BC7", "#08519C")
)

ordered_programs <- hc$labels[hc$order]
ordered_mp <- assignment$metaprogram_id[match(ordered_programs, assignment$program_id)]
block_rle <- rle(as.character(ordered_mp))
block_end <- cumsum(block_rle$lengths)
block_start <- block_end - block_rle$lengths + 1
block_df <- tibble::tibble(
  metaprogram_id = block_rle$values,
  start = block_start,
  end = block_end,
  n_programs_ordered = block_rle$lengths
) %>%
  filter(!is.na(metaprogram_id))

ht <- Heatmap(
  jaccard,
  name = "Jaccard",
  col = jaccard_col,
  cluster_rows = hc,
  cluster_columns = hc,
  show_row_names = FALSE,
  show_column_names = FALSE,
  show_row_dend = TRUE,
  show_column_dend = TRUE,
  row_dend_width = unit(9, "mm"),
  column_dend_height = unit(9, "mm"),
  top_annotation = top_ann,
  border = TRUE,
  rect_gp = gpar(col = NA),
  heatmap_legend_param = list(
    at = c(0, 0.10, 0.20, 0.35),
    labels = c("0", "0.10", "0.20", ">=0.35"),
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7),
    legend_height = unit(28, "mm"),
    grid_width = unit(3, "mm")
  )
)

pdf(file.path(out_dir, "S13_NMF_program_jaccard_clustering.pdf"),
    width = 6.2, height = 6.2, useDingbats = FALSE)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
     padding = unit(c(3, 3, 3, 3), "mm"))
decorate_heatmap_body("Jaccard", {
  n_mat <- nrow(jaccard)
  for (i in seq_len(nrow(block_df))) {
    s <- block_df$start[i]
    e <- block_df$end[i]
    block_n <- e - s + 1
    x_mid <- (s + e - 1) / (2 * n_mat)
    y_mid <- 1 - (s + e - 1) / (2 * n_mat)
    wh <- block_n / n_mat

    grid.rect(
      x = unit(x_mid, "npc"),
      y = unit(y_mid, "npc"),
      width = unit(wh, "npc"),
      height = unit(wh, "npc"),
      gp = gpar(col = "grey15", fill = NA, lwd = 0.75, lty = 2)
    )

  }
})
dev.off()

readr::write_csv(
  assignment %>% arrange(metaprogram_id, patient, program_id),
  file.path(src_dir, "S13_NMF_program_metaprogram_assignment.csv")
)
readr::write_csv(
  block_df %>% mutate(label = mp_label_map[metaprogram_id]),
  file.path(src_dir, "S13_NMF_program_jaccard_diagonal_blocks.csv")
)
write.csv(jaccard, file.path(src_dir, "S13_NMF_program_jaccard_matrix.csv"))

writeLines(capture.output(sessionInfo()), file.path(src_dir, "S8_S13_required_supplement_session_info.txt"))

cat("Required supplement additions completed.\n")
cat("Outputs:\n")
cat("- figures/supplement/S8_S3_S4_within_patient_metaprogram_consistency.pdf\n")
cat("- figures/supplement/S13_NMF_program_jaccard_clustering.pdf\n")
