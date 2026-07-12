# =============================================================================
# F8  GSEA subtype x pathway dot-matrix heatmap
#     Dense 20 x 4 dot heatmap
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(ggh4x)
  library(tibble)
  library(grid)
})

step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir && dir.exists(step_dir)) {
  setwd(step_dir)
}

gsea_path   <- "tables/10e_one_vs_rest_GSEA_combined_wilcoxon_downsampled.csv"
out_pdf     <- "figures/F8_GSEA_subtype_pathway_dotplot.pdf"
out_source  <- "figures/source_data/F8_GSEA_subtype_pathway_dotplot_source.csv"
out_session <- "figures/source_data/F8_GSEA_subtype_pathway_dotplot_session_info.txt"

dir.create("figures",             showWarnings = FALSE, recursive = TRUE)
dir.create("figures/source_data", showWarnings = FALSE, recursive = TRUE)

subtype_levels <- paste0("Subtype", 1:4)
subtype_names <- c(
  "Subtype1" = "Proliferative-NPC",
  "Subtype2" = "OPC-Myelination",
  "Subtype3" = "Vascular-niche MES",
  "Subtype4" = "MES-Antigen-presenting"
)
# 列标题折两行,避免放大字号后被窄 panel 截断(行标题竖排仍用单行)
col_names <- c(
  "Subtype1" = "Proliferative-\nNPC",
  "Subtype2" = "OPC-\nMyelination",
  "Subtype3" = "Vascular-niche\nMES",
  "Subtype4" = "MES-Antigen-\npresenting"
)
subtype_colors <- c(
  "Subtype1" = "#0072B5",
  "Subtype2" = "#E18727",
  "Subtype3" = "#20854E",
  "Subtype4" = "#BC3C29"
)

selected_terms <- tibble::tribble(
  ~subtype_k4, ~ID, ~display_label, ~display_order,
  "Subtype1", "HALLMARK_G2M_CHECKPOINT",                              "G2M checkpoint",                     1,
  "Subtype1", "HALLMARK_E2F_TARGETS",                                 "E2F targets",                        2,
  "Subtype1", "REACTOME_CELL_CYCLE_CHECKPOINTS",                      "Cell cycle checkpoints",             3,
  "Subtype1", "REACTOME_CELL_CYCLE",                                  "Cell cycle",                         4,
  "Subtype1", "REACTOME_S_PHASE",                                     "S phase",                            5,

  "Subtype2", "REACTOME_NEURONAL_SYSTEM",                             "Neuronal system",                    1,
  "Subtype2", "REACTOME_PROTEIN_PROTEIN_INTERACTIONS_AT_SYNAPSES",    "Synaptic interactions",              2,
  "Subtype2", "REACTOME_ION_CHANNEL_TRANSPORT",                       "Ion channel transport",              3,
  "Subtype2", "REACTOME_NEUREXINS_AND_NEUROLIGINS",                   "Neurexins / neuroligins",            4,
  "Subtype2", "REACTOME_METABOLISM_OF_LIPIDS",                        "Lipid metabolism",                   5,

  "Subtype3", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",           "Epithelial-mesenchymal transition",  1,
  "Subtype3", "REACTOME_COLLAGEN_FORMATION",                          "Collagen formation",                 2,
  "Subtype3", "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",           "ECM organization",                   3,
  "Subtype3", "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",          "Integrin cell-surface interactions", 4,
  "Subtype3", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",                   "Oxidative phosphorylation",          5,

  "Subtype4", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",                     "TNFA / NFKB signaling",              1,
  "Subtype4", "HALLMARK_INTERFERON_GAMMA_RESPONSE",                   "Interferon-gamma response",          2,
  "Subtype4", "HALLMARK_INFLAMMATORY_RESPONSE",                       "Inflammatory response",              3,
  "Subtype4", "HALLMARK_IL6_JAK_STAT3_SIGNALING",                     "IL6 / JAK / STAT3 signaling",        4,
  "Subtype4", "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",         "Cytokine signaling",                 5
)

gsea <- readr::read_csv(gsea_path, show_col_types = FALSE) %>%
  mutate(
    subtype_k4 = as.character(subtype_k4),
    NES        = suppressWarnings(as.numeric(NES)),
    p.adjust   = suppressWarnings(as.numeric(p.adjust))
  ) %>%
  distinct(ID, subtype_k4, .keep_all = TRUE)

matrix_data <- selected_terms %>%
  rename(owner_subtype = subtype_k4) %>%
  tidyr::crossing(subtype_k4 = subtype_levels) %>%
  left_join(gsea %>% select(ID, subtype_k4, NES, p.adjust),
            by = c("ID", "subtype_k4")) %>%
  mutate(
    significant   = !is.na(p.adjust) & p.adjust < 0.05 & !is.na(NES) & NES > 0,
    is_own        = (subtype_k4 == owner_subtype),
    neg_log10_q   = ifelse(significant, pmin(-log10(p.adjust + 1e-300), 10), NA_real_),
    NES_show      = ifelse(significant, NES, NA_real_),
    owner_subtype = factor(owner_subtype, levels = subtype_levels),
    subtype_k4    = factor(subtype_k4,    levels = subtype_levels)
  ) %>%
  arrange(owner_subtype, display_order)

missing_own <- matrix_data %>%
  filter(is_own & !significant) %>%
  transmute(pair = paste(owner_subtype, ID, sep = "::"))
if (nrow(missing_own) > 0) {
  stop("Curated GSEA terms missing/non-significant in own subtype: ",
       paste(missing_own$pair, collapse = ", "))
}

matrix_data <- matrix_data %>%
  mutate(display_label = str_wrap(display_label, width = 36))

y_levels <- matrix_data %>%
  filter(subtype_k4 == subtype_levels[1]) %>%
  arrange(owner_subtype, desc(display_order)) %>%
  pull(display_label)

matrix_data <- matrix_data %>%
  mutate(
    display_label = factor(display_label, levels = y_levels),
    owner_strip = factor(
      subtype_names[as.character(owner_subtype)],
      levels = unname(subtype_names[subtype_levels])
    ),
    subtype_strip = factor(
      col_names[as.character(subtype_k4)],
      levels = unname(col_names[subtype_levels])
    ),
    x_pos = 1L
  )

readr::write_csv(matrix_data, out_source)

strip_des <- strip_themed(
  background_x = elem_list_rect(
    fill  = unname(subtype_colors[subtype_levels]),
    color = NA
  ),
  text_x = elem_list_text(
    color = "white", face = "bold", size = 17, lineheight = 0.95
  ),
  background_y = elem_list_rect(
    fill  = unname(subtype_colors[subtype_levels]),
    color = NA
  ),
  text_y = elem_list_text(
    color = "white", face = "bold", size = 16, lineheight = 0.95
  )
)

nes_min <- min(matrix_data$NES_show, na.rm = TRUE)
nes_max <- max(matrix_data$NES_show, na.rm = TRUE)

p <- ggplot(matrix_data, aes(x = x_pos, y = display_label)) +
  geom_point(
    data    = matrix_data %>% filter(significant),
    aes(size = neg_log10_q, fill = NES_show),
    shape   = 21,
    stroke  = 0.5,
    color   = "grey20"
  ) +
  geom_text(
    data = matrix_data %>% filter(is_own),
    aes(label = sprintf("%.1f", NES_show),
        color = ifelse(NES_show > 2.2, "white", "grey15")),
    size     = 5,
    fontface = "bold"
  ) +
  scale_color_identity() +
  scale_fill_gradient(
    low    = "#FFE9B8",
    high   = "#8E1818",
    name   = "NES",
    breaks = pretty(c(nes_min, nes_max), n = 4),
    guide  = guide_colorbar(
      barwidth        = unit(5, "mm"),
      barheight       = unit(34, "mm"),
      ticks.colour    = NA,
      frame.colour    = "grey55",
      frame.linewidth = 0.5,
      title.position  = "top"
    )
  ) +
  scale_size_continuous(
    range  = c(3.5, 11),
    limits = c(0, 10),
    breaks = c(2, 5, 8, 10),
    labels = c("2", "5", "8", ">=10"),
    name   = expression(-log[10](BH~italic(q))),
    guide  = guide_legend(
      override.aes  = list(fill = "grey75", color = "grey25", stroke = 0.3),
      title.position = "top"
    )
  ) +
  scale_x_continuous(breaks = NULL, expand = expansion(add = 0.55)) +
  scale_y_discrete(expand = expansion(add = 0.35)) +
  facet_grid2(
    rows   = vars(owner_strip),
    cols   = vars(subtype_strip),
    scales = "free_y",
    space  = "free_y",
    switch = "y",
    strip  = strip_des
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 19) +
  theme(
    axis.text.x         = element_blank(),
    axis.ticks          = element_blank(),
    axis.line           = element_blank(),
    axis.text.y         = element_text(size = 18, color = "black",
                                       lineheight = 0.9),
    panel.background    = element_rect(fill = "white", color = NA),
    panel.grid          = element_blank(),
    panel.border        = element_rect(fill = NA, color = "grey88", linewidth = 0.4),
    panel.spacing.x     = unit(2, "pt"),
    panel.spacing.y     = unit(2, "pt"),
    plot.background     = element_rect(fill = "white", color = NA),
    strip.placement     = "outside",
    strip.text.x        = element_text(margin = margin(t = 7, b = 7)),
    strip.text.y.left   = element_text(angle = 90, margin = margin(r = 7, l = 4)),
    legend.position     = "right",
    legend.title        = element_text(size = 14, face = "bold"),
    legend.text         = element_text(size = 12.5),
    legend.key.size     = unit(6, "mm"),
    legend.box.spacing  = unit(4, "mm"),
    legend.spacing.y    = unit(6, "mm"),
    plot.margin         = margin(t = 6, r = 6, b = 6, l = 6)
  )

ggsave(out_pdf, p, width = 12.5, height = 7.2, useDingbats = FALSE)

writeLines(capture.output(sessionInfo()), out_session)

cat("F8 GSEA dot-matrix completed.\n")
cat("Output PDF: ", out_pdf, "\n")
cat("Source CSV: ", out_source, "\n\n")

cat("Matrix coverage (significant hits per subtype column):\n")
print(
  matrix_data %>%
    group_by(subtype_k4) %>%
    summarise(
      n_sig   = sum(significant),
      n_own   = sum(is_own & significant),
      n_cross = sum(!is_own & significant),
      .groups = "drop"
    ) %>%
    as.data.frame()
)

cat("\nSelected terms (curated):\n")
print(
  matrix_data %>%
    filter(is_own) %>%
    arrange(owner_subtype, display_order) %>%
    select(owner_subtype, ID, display_label, NES, p.adjust, neg_log10_q) %>%
    as.data.frame()
)
