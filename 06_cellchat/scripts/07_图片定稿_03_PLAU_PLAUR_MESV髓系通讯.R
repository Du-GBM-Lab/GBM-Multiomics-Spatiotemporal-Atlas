#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

module_dir <- "07_细胞通讯"
figures_dir <- file.path(module_dir, "figures")
source_data_dir <- file.path(figures_dir, "source_data")
archive_dir <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/06_细胞通讯"
archive_main_dir <- file.path(archive_dir, "正文图")
archive_source_dir <- file.path(archive_dir, "source_data")

dir.create(archive_main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(archive_source_dir, recursive = TRUE, showWarnings = FALSE)

comm_path <- file.path(module_dir, "tables", "PLAU_PLAUR_纯TME方向表达支撑表.csv")
comm <- read_csv(comm_path, show_col_types = FALSE) %>%
  filter(ligand == "PLAU", receptor == "PLAUR")

focus_groups <- c(
  "Malignant_NPC-P", "Malignant_OPC-M", "Malignant_MES-V", "Malignant_MES-I",
  "Macrophages", "Microglial"
)

plot_df <- expand_grid(source = focus_groups, target = focus_groups) %>%
  left_join(comm, by = c("source", "target")) %>%
  mutate(
    prob = replace_na(prob, 0),
    source_PLAU_mean = replace_na(source_PLAU_mean, NA_real_),
    source_PLAU_pct = replace_na(source_PLAU_pct, NA_real_),
    target_PLAUR_mean = replace_na(target_PLAUR_mean, NA_real_),
    target_PLAUR_pct = replace_na(target_PLAUR_pct, NA_real_),
    source = factor(source, levels = rev(focus_groups)),
    target = factor(target, levels = focus_groups),
    prob_label = sprintf("%.2f", prob * 1e4),
    pair_class = case_when(
      as.character(source) == "Malignant_MES-V" & as.character(target) %in% c("Macrophages", "Microglial") ~ "MES-V to myeloid",
      as.character(source) %in% c("Macrophages", "Microglial") & as.character(target) == "Malignant_MES-V" ~ "myeloid to MES-V",
      as.character(source) %in% c("Macrophages", "Microglial") & as.character(target) %in% c("Macrophages", "Microglial") ~ "myeloid internal",
      TRUE ~ "other"
    ),
    source_label = recode(
      as.character(source),
      "Malignant_NPC-P" = "NPC-P",
      "Malignant_OPC-M" = "OPC-M",
      "Malignant_MES-V" = "MES-V",
      "Malignant_MES-I" = "MES-I",
      "Microglial" = "Microglia",
      .default = as.character(source)
    ),
    target_label = recode(
      as.character(target),
      "Malignant_NPC-P" = "NPC-P",
      "Malignant_OPC-M" = "OPC-M",
      "Malignant_MES-V" = "MES-V",
      "Malignant_MES-I" = "MES-I",
      "Microglial" = "Microglia",
      .default = as.character(target)
    ),
    source_label = factor(source_label, levels = rev(c("NPC-P", "OPC-M", "MES-V", "MES-I", "Macrophages", "Microglia"))),
    target_label = factor(target_label, levels = c("NPC-P", "OPC-M", "MES-V", "MES-I", "Macrophages", "Microglia"))
  )

source_out <- file.path(source_data_dir, "PLAU_PLAUR_MESV髓系通讯_source.csv")
archive_source_out <- file.path(archive_source_dir, "PLAU_PLAUR_MESV髓系通讯_source.csv")
write_csv(plot_df, source_out)
write_csv(plot_df, archive_source_out)

highlight_df <- plot_df %>%
  filter(pair_class %in% c("MES-V to myeloid", "myeloid to MES-V"))

p <- ggplot(plot_df, aes(target_label, source_label)) +
  geom_tile(aes(fill = prob), color = "grey55", linewidth = 0.22) +
  geom_tile(
    data = highlight_df,
    aes(target_label, source_label),
    fill = NA,
    color = "#B2182B",
    linewidth = 0.75
  ) +
  geom_text(aes(label = prob_label), size = 2.7, color = "black") +
  scale_fill_gradient(
    low = "white",
    high = "#B2182B",
    name = "Communication\nprobability",
    labels = label_scientific(digits = 2)
  ) +
  labs(
    title = "PLAU-PLAUR communication in MES-V-myeloid pairs",
    subtitle = "Values indicate communication probability x 10,000",
    x = "Receiver",
    y = "Sender"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(size = 8, color = "grey25", hjust = 0),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title = element_text(size = 9),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(fill = NA, color = "grey35", linewidth = 0.45),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  ) +
  coord_fixed()

out_pdf <- file.path(figures_dir, "[正文]待审核_PLAU_PLAUR_MESV髓系通讯.pdf")
archive_pdf <- file.path(archive_main_dir, "待审核_PLAU_PLAUR_MESV髓系通讯.pdf")

ggsave(out_pdf, p, width = 5.6, height = 4.6, useDingbats = FALSE)
ggsave(archive_pdf, p, width = 5.6, height = 4.6, useDingbats = FALSE)

file.copy(
  file.path(figures_dir, "[正文]CellChat纯TME_PLAU_malignant_to_TME_bubble.pdf"),
  file.path(archive_main_dir, "待审核_CellChat纯TME_PLAU_malignant_to_TME_bubble.pdf"),
  overwrite = TRUE
)
file.copy(
  file.path(figures_dir, "[正文]CellChat纯TME_PLAU_TME_to_malignant_bubble.pdf"),
  file.path(archive_main_dir, "待审核_CellChat纯TME_PLAU_TME_to_malignant_bubble.pdf"),
  overwrite = TRUE
)

message("Updated focused PLAU-PLAUR review figure: ", archive_pdf)
