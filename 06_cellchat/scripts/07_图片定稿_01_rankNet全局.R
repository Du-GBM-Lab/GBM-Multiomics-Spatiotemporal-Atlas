#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

module_dir <- "07_细胞通讯"
figures_dir <- file.path(module_dir, "figures")
source_data_dir <- file.path(figures_dir, "source_data")
archive_dir <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/06_细胞通讯"
archive_main_dir <- file.path(archive_dir, "正文候选")
archive_source_dir <- file.path(archive_dir, "source_data")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(archive_main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(archive_source_dir, recursive = TRUE, showWarnings = FALSE)

rank_path <- file.path(source_data_dir, "CellChat纯TME_global_pathway_rank.csv")
rank_tbl <- read_csv(rank_path, show_col_types = FALSE) %>%
  arrange(rank) %>%
  mutate(
    pathway_name = as.character(pathway_name),
    is_PLAU = pathway_name == "PLAU",
    pathway_plot = factor(pathway_name, levels = rev(pathway_name))
  )

plau <- rank_tbl %>% filter(pathway_name == "PLAU")
if (nrow(plau) != 1) {
  stop("Expected exactly one PLAU row in global pathway rank table.")
}

plot_tbl <- rank_tbl %>% slice_head(n = 30)

p <- ggplot(plot_tbl, aes(x = total_prob, y = pathway_plot)) +
  geom_col(aes(fill = is_PLAU), width = 0.74) +
  geom_text(
    data = plau,
    aes(
      x = total_prob,
      y = pathway_plot,
      label = paste0("PLAU, rank ", rank)
    ),
    inherit.aes = FALSE,
    hjust = -0.08,
    size = 3.6,
    fontface = "bold",
    color = "#B2182B"
  ) +
  scale_fill_manual(values = c("FALSE" = "grey72", "TRUE" = "#B2182B"), guide = "none") +
  scale_x_continuous(
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "Global signaling pathway ranking",
    x = "Total communication probability",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.text.y = element_text(color = "black", size = 9),
    axis.text.x = element_text(color = "black", size = 9),
    axis.title.x = element_text(size = 10),
    axis.line = element_line(linewidth = 0.45),
    axis.ticks = element_line(linewidth = 0.45),
    plot.margin = margin(5.5, 28, 5.5, 5.5)
  ) +
  coord_cartesian(clip = "off")

main_pdf <- file.path(figures_dir, "[正文定稿]PLAU通路全局排序.pdf")
archive_pdf <- file.path(archive_main_dir, "PLAU通路全局排序.pdf")
source_csv <- file.path(source_data_dir, "PLAU通路全局排序_source.csv")
archive_source_csv <- file.path(archive_source_dir, "PLAU通路全局排序_source.csv")

ggsave(main_pdf, p, width = 4.3, height = 5.8, useDingbats = FALSE)
ggsave(archive_pdf, p, width = 4.3, height = 5.8, useDingbats = FALSE)

write_csv(plot_tbl, source_csv)
write_csv(plot_tbl, archive_source_csv)

message("Finalized global rankNet figure: ", archive_pdf)
