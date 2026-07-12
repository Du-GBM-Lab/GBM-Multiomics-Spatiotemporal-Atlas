suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(tidyr)
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
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# S12: AMS vs UCell concordance
ams <- readr::read_csv("tables/04_AMS_UCell_spearman.csv", show_col_types = FALSE) %>%
  mutate(state = factor(state, levels = state))

p_ams <- ggplot(ams, aes(state, spearman)) +
  geom_col(fill = "#4D4D4D", width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", spearman)), vjust = -0.35, size = 2.5) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Spearman correlation", title = "AMS and UCell score concordance") +
  theme_supp()

readr::write_csv(ams, file.path(src_dir, "S12_AMS_UCell_spearman.csv"))
ggsave(file.path(out_dir, "S12_AMS_UCell_concordance.pdf"), p_ams,
       width = 4.8, height = 3.0, useDingbats = FALSE)

# S13: fine cluster to subtype mapping
map_path <- if (file.exists("figures/source_data/F1_cluster_subtype_mapping.csv")) {
  "figures/source_data/F1_cluster_subtype_mapping.csv"
} else {
  "tables/03c_cluster_to_subtype_k4_theta5.csv"
}
mapping <- readr::read_csv(map_path, show_col_types = FALSE)
if ("fine_cluster" %in% names(mapping)) {
  mapping <- mapping %>% rename(cluster = fine_cluster)
}
mapping <- mapping %>%
  mutate(
    cluster = factor(as.character(cluster), levels = as.character(sort(as.numeric(as.character(cluster))))),
    subtype = factor(subtype, levels = subtype_levels)
  )

p_map <- ggplot(mapping, aes(cluster, 1, fill = subtype)) +
  geom_tile(color = "white", linewidth = 0.45, height = 0.7) +
  geom_text(aes(label = cluster), size = 2.2, color = "white", fontface = "bold") +
  scale_fill_manual(values = subtype_colors) +
  labs(x = "Fine cluster", y = NULL, fill = "Subtype",
       title = "Fine cluster to final subtype mapping") +
  theme_supp() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")

readr::write_csv(mapping, file.path(src_dir, "S13_fine_cluster_to_subtype_mapping.csv"))
ggsave(file.path(out_dir, "S13_fine_cluster_to_subtype_mapping.pdf"), p_map,
       width = 8.0, height = 2.0, useDingbats = FALSE)

writeLines(capture.output(sessionInfo()), file.path(src_dir, "S12_S13_supplement_extra_session_info.txt"))

cat("S12/S13 supplement extras completed.\n")
print(list.files(out_dir, pattern = "^S1[23].*\\.pdf$", full.names = FALSE))
