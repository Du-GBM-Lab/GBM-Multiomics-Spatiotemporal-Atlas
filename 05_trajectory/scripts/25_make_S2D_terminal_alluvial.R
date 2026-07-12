options(stringsAsFactors = FALSE)

## ======================================================================
## Figure S2D | Slingshot lineage -> Monocle3 terminal correspondence
## Alluvial panel, redesigned to journal-figure standards.
## ======================================================================

if (basename(getwd()) == "06_恶性细胞拟时序") {
  project_root <- normalizePath(file.path(getwd(), ".."), winslash = "/", mustWork = TRUE)
} else {
  project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
setwd(file.path(project_root, "06_恶性细胞拟时序"))

dir.create("figures/final_panels/supplement", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/final_panels/source_data/supplement", recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggalluvial)
})

ARROW <- "→"
DOT <- "·"

lineage_levels <- c("L1", "L2", "L3")
lineage_path <- c(
  L1 = paste("NPC-P", ARROW, "MES-I"),
  L2 = paste("NPC-P", ARROW, "MES-V"),
  L3 = paste("NPC-P", ARROW, "OPC-M")
)

lineage_pal <- c(L1 = "#BC3C29", L2 = "#20854E", L3 = "#E18727")
band_grey <- "#E8E7E4"

terminal_order <- c("Y_26", "Y_70", "Y_72", "Y_41", "Y_98", "Y_43", "Y_6", "Y_31")
terminal_note <- c(
  Y_26 = "MES-I",
  Y_70 = "MES-mix",
  Y_72 = "MES-mix",
  Y_41 = "NPC/MES",
  Y_98 = "MES-V",
  Y_43 = "NPC leaf",
  Y_6 = "NPC leaf",
  Y_31 = "OPC-M"
)

flow_raw <- read_csv("figures/source_data/03_panel_D.csv", show_col_types = FALSE)
terminal_detection <- read_csv("tables/monocle3_terminal_detection.csv", show_col_types = FALSE)

terminal_summary <- terminal_detection |>
  group_by(monocle3_terminal, terminal_total_cells) |>
  slice_max(terminal_fraction, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(
    monocle3_terminal,
    terminal_total_cells,
    subtype_short,
    terminal_fraction,
    terminal_note = terminal_note[monocle3_terminal]
  ) |>
  mutate(monocle3_terminal = factor(monocle3_terminal, levels = terminal_order)) |>
  arrange(monocle3_terminal)

flow <- flow_raw |>
  transmute(
    lineage = recode(
      slingshot_assigned_lineage,
      lineage1 = "L1",
      lineage2 = "L2",
      lineage3 = "L3"
    ),
    terminal = monocle3_terminal,
    n_cells = n_cells
  ) |>
  left_join(terminal_summary, by = c("terminal" = "monocle3_terminal")) |>
  mutate(
    lineage = factor(lineage, levels = lineage_levels),
    terminal = factor(terminal, levels = terminal_order)
  )

total_cells <- sum(flow$n_cells)
lineage_n <- tapply(flow$n_cells, flow$lineage, sum)
terminal_n <- tapply(flow$n_cells, flow$terminal, sum)

fmt_n <- function(x) format(as.integer(x), big.mark = ",", trim = TRUE)

lineage_detail <- setNames(
  paste0(lineage_path[names(lineage_n)], "   ", DOT, "   ", fmt_n(lineage_n), " cells"),
  names(lineage_n)
)

terminal_frac_pct <- round(terminal_summary$terminal_fraction * 100)
terminal_label_full <- setNames(
  paste0(
    names(terminal_n), "    ",
    terminal_note[names(terminal_n)], "   ", DOT, "   ",
    fmt_n(terminal_n), " cells (", terminal_frac_pct, "%)"
  ),
  names(terminal_n)
)

p <- ggplot(flow, aes(axis1 = lineage, axis2 = terminal, y = n_cells)) +
  geom_alluvium(
    aes(fill = lineage),
    width = 0.15,
    alpha = 0.55,
    knot.pos = 0.34,
    color = NA
  ) +
  geom_stratum(
    aes(fill = after_stat(ifelse(x == 1, as.character(stratum), "TERM"))),
    width = 0.15,
    color = "white",
    linewidth = 0.7
  ) +
  geom_text(
    stat = "stratum",
    na.rm = TRUE,
    fontface = "bold",
    size = 2.85,
    color = "grey12",
    hjust = 1,
    nudge_x = -0.115,
    nudge_y = total_cells * 0.016,
    aes(label = after_stat(ifelse(x == 1, as.character(stratum), NA)))
  ) +
  geom_text(
    stat = "stratum",
    na.rm = TRUE,
    size = 2.35,
    color = "grey42",
    hjust = 1,
    nudge_x = -0.115,
    nudge_y = -total_cells * 0.017,
    aes(label = after_stat(ifelse(x == 1, lineage_detail[as.character(stratum)], NA)))
  ) +
  geom_text(
    stat = "stratum",
    na.rm = TRUE,
    size = 2.40,
    color = "grey18",
    hjust = 0,
    nudge_x = 0.115,
    aes(label = after_stat(ifelse(x == 2, terminal_label_full[as.character(stratum)], NA)))
  ) +
  annotate(
    "text",
    x = 1,
    y = total_cells * 1.05,
    label = "Slingshot lineage",
    fontface = "bold",
    size = 3.05,
    color = "grey20"
  ) +
  annotate(
    "text",
    x = 2,
    y = total_cells * 1.05,
    label = "Monocle3 terminal",
    fontface = "bold",
    size = 3.05,
    color = "grey20"
  ) +
  scale_fill_manual(values = c(lineage_pal, TERM = band_grey), guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.10))) +
  coord_cartesian(xlim = c(0, 3), clip = "off") +
  labs(
    title = paste("Slingshot lineage", ARROW, "Monocle3 terminal correspondence"),
    subtitle = paste0(
      "Ribbon width proportional to cell number; terminals ordered by dominant subtype context (n = ",
      fmt_n(total_cells),
      " malignant cells)"
    )
  ) +
  theme_void(base_size = 8) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(face = "bold", size = 9.5, hjust = 0, margin = margin(b = 1.5)),
    plot.subtitle = element_text(size = 7, color = "grey38", margin = margin(b = 10)),
    plot.title.position = "plot",
    legend.position = "none",
    plot.margin = margin(t = 9, r = 16, b = 9, l = 14)
  )

out_pdf <- "figures/final_panels/supplement/S2D_lineage_terminal_alluvial.pdf"
out_png <- "figures/final_panels/supplement/S2D_lineage_terminal_alluvial.png"

ggsave(out_pdf, p, width = 7.6, height = 5.6, device = cairo_pdf, bg = "white")
ggsave(out_png, p, width = 7.6, height = 5.6, dpi = 600, bg = "white")

flow_summary <- flow |>
  group_by(lineage) |>
  mutate(
    lineage_total = sum(n_cells),
    lineage_fraction = n_cells / lineage_total
  ) |>
  ungroup() |>
  arrange(lineage, desc(n_cells))

readr::write_csv(
  flow_summary,
  "figures/final_panels/source_data/supplement/S2D_lineage_terminal_alluvial.csv"
)
readr::write_csv(
  terminal_summary,
  "figures/final_panels/source_data/supplement/S2D_terminal_order_summary.csv"
)

message("Wrote ", out_pdf)
