#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

module_dir <- "07_细胞通讯"
archive_dir <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/06_细胞通讯"
main_dir <- file.path(archive_dir, "正文图")
supp_dir <- file.path(archive_dir, "补充图")
source_dir <- file.path(archive_dir, "source_data")

dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

comm_path <- file.path(module_dir, "tables", "CellChat纯TME_全部通讯结果.csv")
stopifnot(file.exists(comm_path))

subtype_levels <- c(
  "Malignant_NPC-P", "Malignant_OPC-M", "Malignant_MES-V", "Malignant_MES-I"
)
tme_levels <- c(
  "Macrophages", "Microglial", "Monocytes", "cDCs", "pDCs",
  "T cells", "NK cells", "B cells", "Endothelial", "Mural cells"
)

subtype_labels <- c(
  "Malignant_NPC-P" = "NPC-P",
  "Malignant_OPC-M" = "OPC-M",
  "Malignant_MES-V" = "MES-V",
  "Malignant_MES-I" = "MES-I"
)
tme_labels <- c(
  "Macrophages" = "Macrophages",
  "Microglial" = "Microglia",
  "Monocytes" = "Monocytes",
  "cDCs" = "cDCs",
  "pDCs" = "pDCs",
  "T cells" = "T cells",
  "NK cells" = "NK cells",
  "B cells" = "B cells",
  "Endothelial" = "Endothelial",
  "Mural cells" = "Mural"
)

comm <- read_csv(comm_path, show_col_types = FALSE)

malignant_to_tme <- comm %>%
  filter(source %in% subtype_levels, target %in% tme_levels) %>%
  mutate(
    subtype = factor(unname(subtype_labels[source]), levels = unname(subtype_labels[subtype_levels])),
    target_label = factor(unname(tme_labels[target]), levels = unname(tme_labels[tme_levels])),
    lr_pair = interaction_name_2,
    lr_short = paste0(ligand, " -> ", receptor),
    is_plau_plaur = ligand == "PLAU" & receptor == "PLAUR"
  )

top_lr <- malignant_to_tme %>%
  filter(source == "Malignant_MES-V") %>%
  group_by(lr_pair, lr_short, pathway_name) %>%
  summarise(
    mesv_total_prob = sum(prob, na.rm = TRUE),
    n_edges = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mesv_total_prob)) %>%
  mutate(mesv_rank = row_number()) %>%
  slice_head(n = 35)

selected_lr <- top_lr

display_order <- selected_lr %>%
  arrange(desc(mesv_total_prob)) %>%
  pull(lr_short)

plot_data <- malignant_to_tme %>%
  semi_join(selected_lr, by = "lr_pair") %>%
  group_by(subtype, target_label, lr_pair, lr_short, pathway_name, ligand, receptor, is_plau_plaur) %>%
  summarise(
    prob = sum(prob, na.rm = TRUE),
    min_pval = min(pval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    selected_lr %>% select(lr_pair, mesv_lr_prob = mesv_total_prob, mesv_lr_rank = mesv_rank),
    by = "lr_pair"
  ) %>%
  mutate(
    display_rank_top_to_bottom = match(as.character(lr_short), display_order),
    neg_log10_p = pmin(-log10(pmax(min_pval, 1e-300)), 50),
    lr_short = factor(lr_short, levels = rev(display_order))
  )

write_csv(plot_data, file.path(source_dir, "四亚型TME配受体散点图_source.csv"))

p <- ggplot(plot_data, aes(x = target_label, y = lr_short)) +
  geom_rect(
    data = data.frame(
      subtype = factor("MES-V", levels = unname(subtype_labels[subtype_levels])),
      xmin = 0.5,
      xmax = length(tme_levels) + 0.5,
      ymin = match("PLAU -> PLAUR", levels(plot_data$lr_short)) - 0.36,
      ymax = match("PLAU -> PLAUR", levels(plot_data$lr_short)) + 0.36
    ),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA,
    color = "#BC3C29",
    linewidth = 0.28
  ) +
  geom_point(
    aes(size = prob, color = prob),
    alpha = 0.9
  ) +
  facet_grid(. ~ subtype, scales = "free_x", space = "free_x") +
  scale_color_gradientn(
    colors = c("#7FB7D9", "#F4A261", "#B2182B"),
    trans = "sqrt",
    labels = scientific
  ) +
  scale_size_continuous(range = c(0.32, 1.8), labels = scientific) +
  labs(
    x = NULL,
    y = NULL,
    color = "Communication\nprobability",
    size = "Communication\nprobability"
  ) +
  theme_bw(base_size = 5.1) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.1),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey55", linewidth = 0.2),
    strip.text = element_text(face = "bold", size = 4.8),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 3.9),
    axis.text.y = element_text(size = 3.8),
    legend.position = "right",
    legend.key.height = unit(0.14, "cm"),
    legend.key.width = unit(0.12, "cm"),
    legend.text = element_text(size = 3.8),
    legend.title = element_text(size = 4),
    panel.spacing.x = unit(0.04, "lines"),
    plot.margin = margin(1.5, 3.5, 1.5, 1.5)
  )

ggsave(
  file.path(main_dir, "待审核_四亚型TME配受体散点图.pdf"),
  p,
  width = 5.45,
  height = 3.25,
  useDingbats = FALSE
)

subtype_diff <- malignant_to_tme %>%
  semi_join(selected_lr, by = "lr_pair") %>%
  group_by(subtype, lr_pair, lr_short, pathway_name, ligand, receptor, is_plau_plaur) %>%
  summarise(total_prob = sum(prob, na.rm = TRUE), .groups = "drop") %>%
  complete(
    subtype = factor(unname(subtype_labels[subtype_levels]), levels = unname(subtype_labels[subtype_levels])),
    lr_short = factor(display_order, levels = rev(display_order)),
    fill = list(total_prob = 0)
  ) %>%
  group_by(lr_short) %>%
  mutate(
    mean_prob = mean(total_prob, na.rm = TRUE),
    sd_prob = sd(total_prob, na.rm = TRUE),
    z_score = if_else(sd_prob > 0, (total_prob - mean_prob) / sd_prob, 0),
    z_score = pmax(pmin(z_score, 2), -2)
  ) %>%
  ungroup()

specificity_wide_for_order <- subtype_diff %>%
  select(lr_short, subtype, z_score) %>%
  pivot_wider(names_from = subtype, values_from = z_score, values_fill = 0)

specificity_order_mat <- specificity_wide_for_order %>%
  select(all_of(unname(subtype_labels[subtype_levels]))) %>%
  as.matrix()
rownames(specificity_order_mat) <- specificity_wide_for_order$lr_short

if (nrow(specificity_order_mat) > 1) {
  clustered_lr <- rownames(specificity_order_mat)[hclust(dist(specificity_order_mat), method = "ward.D2")$order]
} else {
  clustered_lr <- rownames(specificity_order_mat)
}

subtype_diff <- subtype_diff %>%
  mutate(lr_short = factor(as.character(lr_short), levels = clustered_lr))

write_csv(subtype_diff, file.path(source_dir, "四亚型配受体差异热图_source.csv"))

p_diff <- ggplot(subtype_diff, aes(x = subtype, y = lr_short, fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.18) +
  scale_fill_gradient2(
    low = "#3B6EA8",
    mid = "white",
    high = "#BC3C29",
    midpoint = 0,
    limits = c(-2, 2),
    name = "Row z-score"
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", size = 7),
    axis.text.y = element_text(size = 5.8),
    legend.position = "right",
    plot.margin = margin(5, 8, 5, 5)
  )

ggsave(
  file.path(supp_dir, "四亚型配受体差异热图.pdf"),
  p_diff,
  width = 4.2,
  height = 6.4,
  useDingbats = FALSE
)

specificity_dot <- subtype_diff %>%
  group_by(lr_short) %>%
  mutate(
    other_mean_prob = vapply(seq_along(total_prob), function(i) {
      mean(total_prob[-i], na.rm = TRUE)
    }, numeric(1)),
    log2_enrichment_vs_others = log2((total_prob + 1e-6) / (other_mean_prob + 1e-6)),
    positive_enrichment = pmax(log2_enrichment_vs_others, 0),
    positive_enrichment = pmin(positive_enrichment, 6)
  ) %>%
  ungroup()

specificity_wide <- specificity_dot %>%
  select(lr_short, subtype, log2_enrichment_vs_others) %>%
  pivot_wider(names_from = subtype, values_from = log2_enrichment_vs_others, values_fill = 0)

specificity_mat <- specificity_wide %>%
  select(all_of(unname(subtype_labels[subtype_levels]))) %>%
  as.matrix()
rownames(specificity_mat) <- specificity_wide$lr_short

specificity_dot <- specificity_dot %>%
  mutate(lr_short = factor(as.character(lr_short), levels = clustered_lr))

write_csv(specificity_dot, file.path(source_dir, "四亚型配受体差异散点图_source.csv"))

p_specificity <- ggplot(specificity_dot, aes(x = subtype, y = lr_short)) +
  geom_point(
    aes(color = total_prob, size = positive_enrichment),
    alpha = 0.9
  ) +
  geom_rect(
    data = data.frame(
      subtype = factor("MES-V", levels = unname(subtype_labels[subtype_levels])),
      xmin = which(levels(specificity_dot$subtype) == "MES-V") - 0.42,
      xmax = which(levels(specificity_dot$subtype) == "MES-V") + 0.42,
      ymin = match("PLAU -> PLAUR", levels(specificity_dot$lr_short)) - 0.36,
      ymax = match("PLAU -> PLAUR", levels(specificity_dot$lr_short)) + 0.36
    ),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = NA,
    color = "#BC3C29",
    linewidth = 0.28
  ) +
  scale_color_gradientn(
    colors = c("#7FB7D9", "#F4A261", "#B2182B"),
    trans = "sqrt",
    labels = scientific,
    name = "Total\nprobability"
  ) +
  scale_size_continuous(
    range = c(0.35, 2.6),
    limits = c(0, 6),
    name = "Subtype\nenrichment"
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 5.8) +
  theme(
    panel.grid.major = element_line(color = "grey88", linewidth = 0.1),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold", size = 5.1),
    axis.text.y = element_text(size = 4.55),
    legend.position = "right",
    legend.key.height = unit(0.16, "cm"),
    legend.key.width = unit(0.13, "cm"),
    legend.text = element_text(size = 3.8),
    legend.title = element_text(size = 4),
    plot.margin = margin(1.5, 3.5, 1.5, 1.5)
  )

ggsave(
  file.path(main_dir, "待审核_四亚型配受体差异散点图.pdf"),
  p_specificity,
  width = 4.05,
  height = 3.25,
  useDingbats = FALSE
)
