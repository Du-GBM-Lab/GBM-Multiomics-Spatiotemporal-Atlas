suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(scales)
})

step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir) {
  setwd(file.path(getwd(), step_dir))
}

ora_path <- "figures/source_data/F6_ORA_source.csv"
out_pdf <- "figures/F6_GO_BP_per_subtype.pdf"
out_source <- "figures/source_data/F6_GO_BP_per_subtype_source.csv"
out_myelin_audit <- "figures/source_data/F6_Subtype2_myelin_term_audit.csv"
out_session <- "figures/source_data/F6_GO_BP_publication_v2_session_info.txt"

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/source_data", showWarnings = FALSE, recursive = TRUE)

subtype_levels <- paste0("Subtype", 1:4)
subtype_colors <- c(
  "Subtype1" = "#0072B5",
  "Subtype2" = "#E18727",
  "Subtype3" = "#20854E",
  "Subtype4" = "#BC3C29"
)

ora <- readr::read_csv(ora_path, show_col_types = FALSE) %>%
  mutate(
    subtype_k4 = factor(subtype_k4, levels = subtype_levels),
    p.adjust = suppressWarnings(as.numeric(p.adjust)),
    Count = suppressWarnings(as.numeric(Count)),
    neg_log10_BH_q = -log10(p.adjust + 1e-300)
  )

go <- ora %>%
  filter(database == "GO_BP", !is.na(p.adjust), p.adjust < 0.05) %>%
  arrange(subtype_k4, p.adjust, desc(neg_log10_BH_q))

keyword_regex <- "(?i)myelin|ensheathment|oligodendrocyte"
s2_ranked <- go %>%
  filter(subtype_k4 == "Subtype2") %>%
  arrange(p.adjust, desc(neg_log10_BH_q)) %>%
  mutate(numeric_rank = row_number())

s2_myelin_audit <- s2_ranked %>%
  filter(grepl(keyword_regex, Description)) %>%
  select(subtype_k4, ID, Description, numeric_rank, p.adjust, Count, neg_log10_BH_q)
readr::write_csv(s2_myelin_audit, out_myelin_audit)

select_terms_for_subtype <- function(st) {
  x <- go %>%
    filter(subtype_k4 == st) %>%
    arrange(p.adjust, desc(neg_log10_BH_q)) %>%
    mutate(numeric_rank = row_number())

  if (st != "Subtype2") {
    return(x %>% slice_head(n = 10) %>% mutate(manual_promotion = FALSE, selection_rule = "top10_by_BH_q"))
  }

  base <- x %>% slice_head(n = 8) %>% mutate(manual_promotion = FALSE, selection_rule = "top8_by_BH_q")
  promoted <- x %>%
    filter(numeric_rank <= 20, grepl(keyword_regex, Description)) %>%
    filter(!ID %in% base$ID) %>%
    mutate(
      promotion_priority = case_when(
        Description == "ensheathment of neurons" ~ 1L,
        Description == "myelination" ~ 2L,
        Description == "axon ensheathment" ~ 3L,
        grepl("oligodendrocyte", Description, ignore.case = TRUE) ~ 4L,
        TRUE ~ 9L
      )
    ) %>%
    arrange(promotion_priority, p.adjust, desc(neg_log10_BH_q)) %>%
    slice_head(n = 2) %>%
    select(-promotion_priority) %>%
    mutate(manual_promotion = TRUE, selection_rule = "manual_promote_myelin_related_from_top20")

  bind_rows(base, promoted) %>%
    arrange(manual_promotion, p.adjust, desc(neg_log10_BH_q))
}

plot_tbl <- bind_rows(lapply(subtype_levels, select_terms_for_subtype)) %>%
  mutate(
    subtype_k4 = factor(subtype_k4, levels = subtype_levels),
    Description_wrapped = stringr::str_wrap(Description, width = 36),
    display_order = ave(p.adjust, subtype_k4, FUN = function(v) rank(v, ties.method = "first"))
  ) %>%
  group_by(subtype_k4) %>%
  arrange(desc(neg_log10_BH_q), .by_group = TRUE) %>%
  mutate(y_id = row_number()) %>%
  ungroup()

readr::write_csv(plot_tbl, out_source)

plot_one_subtype <- function(st) {
  df <- plot_tbl %>%
    filter(subtype_k4 == st) %>%
    arrange(neg_log10_BH_q) %>%
    mutate(
      y = row_number(),
      y_lab = factor(Description_wrapped, levels = Description_wrapped),
      zebra = y %% 2 == 0
    )

  x_max <- max(df$neg_log10_BH_q, na.rm = TRUE)
  x_pad <- max(1.2, x_max * 0.18)
  strip_y <- nrow(df) + 1.2
  strip_half_height <- 0.34

  ggplot(df, aes(x = neg_log10_BH_q, y = y_lab)) +
    geom_rect(
      data = df %>% filter(zebra),
      aes(ymin = as.numeric(y_lab) - 0.5, ymax = as.numeric(y_lab) + 0.5),
      xmin = -Inf, xmax = Inf,
      inherit.aes = FALSE,
      fill = "#F5F5F5", color = NA
    ) +
    geom_col(
      fill = subtype_colors[[st]],
      color = scales::alpha(subtype_colors[[st]], 0.65),
      linewidth = 0.45,
      width = 0.62
    ) +
    geom_text(
      aes(label = paste0("n=", Count), x = neg_log10_BH_q + x_pad * 0.18),
      hjust = 0, size = 4.4, color = "grey30"
    ) +
    annotate(
      "rect",
      xmin = 0, xmax = x_max + x_pad,
      ymin = strip_y - strip_half_height,
      ymax = strip_y + strip_half_height,
      fill = subtype_colors[[st]], color = NA
    ) +
    annotate(
      "text",
      x = (x_max + x_pad) / 2,
      y = strip_y,
      label = sub("Subtype", "Subtype ", st),
      color = "white", fontface = "bold", size = 5.8
    ) +
    scale_x_continuous(
      limits = c(0, x_max + x_pad),
      breaks = scales::pretty_breaks(n = 5),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 17) +
    theme(
      axis.text.y = element_text(size = 16, color = "black"),
      axis.text.x = element_text(size = 14, color = "black"),
      axis.title = element_blank(),
      axis.line = element_line(linewidth = 0.5, color = "grey50"),
      axis.ticks = element_line(linewidth = 0.45, color = "grey50"),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      plot.margin = margin(t = 22, r = 16, b = 10, l = 6)
    )
}

plots <- lapply(subtype_levels, plot_one_subtype)
combined <- (plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]]) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(4, 4, 4, 4)
    )
  ) &
  theme(plot.margin = margin(4, 4, 4, 4))

combined <- combined + plot_layout(guides = "collect") &
  theme(legend.position = "none")

pdf(out_pdf, width = 12.5, height = 10.5, useDingbats = FALSE)
print(combined)
grid::grid.text(
  expression(-log[10]("BH q")),
  x = 0.53,
  y = unit(0.02, "npc"),
  gp = grid::gpar(fontsize = 17)
)
dev.off()

writeLines(capture.output(sessionInfo()), out_session)

cat("F6 GO BP publication v2 completed.\n")
cat("Subtype2 myelin-related hits in top20:\n")
print(s2_myelin_audit %>% filter(numeric_rank <= 20))
cat("Manual promotions used:\n")
print(plot_tbl %>% filter(manual_promotion) %>% select(subtype_k4, Description, numeric_rank, p.adjust, Count))
cat("Output PDF:", out_pdf, "\n")
cat("Source:", out_source, "\n")
