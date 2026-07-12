# =============================================================================
# F7  Patient x malignant-subtype composition + clinical annotations
#     Circular (circos-style) plot - v2
#
# Required packages (install once if missing):
#   install.packages(c("qs2", "dplyr", "tidyr", "readr", "tibble", "circlize"))
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install("ComplexHeatmap")
# =============================================================================

suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
  library(gridBase)
})

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir) {
  setwd(file.path(getwd(), step_dir))
}

obj_path     <- "outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"
out_pdf      <- "figures/F7_patient_composition_clinical_circular.pdf"
out_comp     <- "figures/source_data/F7_patient_subtype_composition.csv"
out_clinical <- "figures/source_data/F7_patient_clinical_annotation.csv"
out_sanity   <- "figures/source_data/F7_patient_composition_clinical_sanity_checks.csv"
out_session  <- "figures/source_data/F7_patient_composition_clinical_session_info.txt"

dir.create("figures",             showWarnings = FALSE, recursive = TRUE)
dir.create("figures/source_data", showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Subtype config  (descriptive names only)
# -----------------------------------------------------------------------------
subtype_levels <- c(
  "Proliferative-NPC",
  "OPC-Myelination",
  "Vascular-niche MES",
  "MES-Antigen-presenting"
)

subtype_map <- c(
  "Subtype1" = "Proliferative-NPC",
  "Subtype2" = "OPC-Myelination",
  "Subtype3" = "Vascular-niche MES",
  "Subtype4" = "MES-Antigen-presenting"
)

subtype_colors <- c(
  "Proliferative-NPC"      = "#0072B5",
  "OPC-Myelination"        = "#E18727",
  "Vascular-niche MES"     = "#20854E",
  "MES-Antigen-presenting" = "#BC3C29"
)

# -----------------------------------------------------------------------------
# Clinical config
# -----------------------------------------------------------------------------
clinical_cols <- c("diagnosis", "WHO", "gender", "IDH", "age",
                   "treatment_1", "treatment_3")

clinical_label_map <- c(
  diagnosis   = "Diagnosis",
  WHO         = "WHO grade",
  gender      = "Sex",
  IDH         = "IDH",
  age         = "Age",
  treatment_1 = "Treatment",
  treatment_3 = "Response"
)

clinical_palettes <- list(
  diagnosis   = c(recurrence   = "#7B3294",
                  ND           = "#A6DBA0"),
  WHO         = c(four         = "#D6604D",
                  three        = "#4393C3"),
  gender      = c(male         = "#4D4D4D",
                  female       = "#BDBDBD"),
  IDH         = c(WT           = "#B2182B",
                  mutant       = "#2166AC"),
  age         = c(below60      = "#92C5DE",
                  above60      = "#F4A582",
                  unkown       = "#D9D9D9",
                  unknown      = "#D9D9D9"),
  treatment_1 = c(TMZ          = "#1B9E77",
                  untreated    = "#D9D9D9",
                  PD1          = "#7570B3"),
  treatment_3 = c(Rec          = "#7B3294",
                  Untreated    = "#D9D9D9",
                  nonresponder = "#E7298A",
                  responder    = "#66A61E")
)

NA_COLOR <- "#EEEEEE"

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
obj <- qs2::qs_read(obj_path)
md  <- obj@meta.data

required_cols <- c("Pt_number", "subtype_k4", clinical_cols)
missing_cols  <- setdiff(required_cols, colnames(md))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

md <- md %>%
  mutate(
    Pt_number     = as.character(Pt_number),
    subtype_label = factor(subtype_map[as.character(subtype_k4)],
                           levels = subtype_levels)
  )

# -----------------------------------------------------------------------------
# Composition table
# -----------------------------------------------------------------------------
comp <- md %>%
  count(Pt_number, subtype_label, name = "n_cells") %>%
  complete(
    Pt_number,
    subtype_label = factor(subtype_levels, levels = subtype_levels),
    fill = list(n_cells = 0)
  ) %>%
  group_by(Pt_number) %>%
  mutate(
    patient_total_cells = sum(n_cells),
    fraction = ifelse(patient_total_cells > 0, n_cells / patient_total_cells, NA_real_)
  ) %>%
  ungroup()

patient_summary <- comp %>%
  group_by(Pt_number) %>%
  summarise(
    patient_total_cells = first(patient_total_cells),
    S1_frac = fraction[subtype_label == "Proliferative-NPC"],
    S2_frac = fraction[subtype_label == "OPC-Myelination"],
    S3_frac = fraction[subtype_label == "Vascular-niche MES"],
    S4_frac = fraction[subtype_label == "MES-Antigen-presenting"],
    .groups = "drop"
  )

# ===== ORDERING =====
#   primary  : Proliferative-NPC fraction -- DESCENDING
#   secondary: OPC-Myelination fraction   -- ASCENDING
#   tiebreak : Pt_number alphabetical
patient_order <- patient_summary %>%
  arrange(desc(S1_frac), S2_frac, Pt_number) %>%
  pull(Pt_number)

comp <- comp %>%
  mutate(
    Pt_number     = factor(Pt_number, levels = patient_order),
    subtype_label = factor(subtype_label, levels = subtype_levels)
  )

# -----------------------------------------------------------------------------
# Clinical lookup
# -----------------------------------------------------------------------------
clinical <- md %>%
  select(Pt_number, all_of(clinical_cols)) %>%
  distinct() %>%
  group_by(Pt_number) %>%
  summarise(
    across(all_of(clinical_cols),
           ~ paste(unique(as.character(.x)), collapse = " | ")),
    .groups = "drop"
  ) %>%
  mutate(Pt_number = factor(Pt_number, levels = patient_order)) %>%
  arrange(Pt_number)

clinical_vec <- lapply(clinical_cols, function(v) {
  setNames(clinical[[v]], as.character(clinical$Pt_number))
})
names(clinical_vec) <- clinical_cols

readr::write_csv(comp     %>% mutate(Pt_number = as.character(Pt_number)), out_comp)
readr::write_csv(clinical %>% mutate(Pt_number = as.character(Pt_number)), out_clinical)

# =============================================================================
# Circular plot
# =============================================================================
n_pt    <- length(patient_order)
GAP_DEG <- 16

gap_after       <- rep(0.28, n_pt)
gap_after[n_pt] <- GAP_DEG

pdf(out_pdf, width = 13, height = 9.5, useDingbats = FALSE)
grid.newpage()

circle_size <- unit(9.15, "in")
pushViewport(viewport(
  x      = unit(0.3, "in"),
  y      = unit(0.5, "npc"),
  width  = circle_size,
  height = circle_size,
  just   = c("left", "center")
))
par(omi = gridOMI(), mar = c(0, 0, 0, 0))

circos.clear()
circos.par(
  cell.padding            = c(0, 0, 0, 0),
  gap.after               = gap_after,
  start.degree            = 90 - GAP_DEG / 2,
  track.margin            = c(0.003, 0.003),
  points.overflow.warning = FALSE
)

circos.initialize(
  factors = factor(patient_order, levels = patient_order),
  xlim    = c(0, 1)
)

circos.track(
  ylim         = c(0, 1),
  bg.border    = NA,
  track.height = 0.055,
  panel.fun    = function(x, y) {
    circos.text(0.5, 0.5, CELL_META$sector.index,
                facing     = "clockwise",
                niceFacing = TRUE,
                cex        = 0.76,
                font       = 2,
                col        = "grey10")
  }
)

circos.track(
  ylim         = c(0, 1),
  bg.border    = "grey75",
  bg.col       = NA,
  track.height = 0.32,
  panel.fun    = function(x, y) {
    pt <- CELL_META$sector.index
    pd <- comp[comp$Pt_number == pt, ]
    pd <- pd[order(pd$subtype_label), ]
    cum <- 0
    for (i in seq_len(nrow(pd))) {
      nxt <- cum + pd$fraction[i]
      circos.rect(0, cum, 1, nxt,
                  col    = subtype_colors[as.character(pd$subtype_label[i])],
                  border = "white",
                  lwd    = 0.45)
      cum <- nxt
    }
  }
)

circos.par(track.margin = c(0.002, 0.010))

for (i in seq_along(clinical_cols)) {
  v   <- clinical_cols[i]
  pal <- clinical_palettes[[v]]
  vec <- clinical_vec[[v]]

  local({
    pal_local <- pal
    vec_local <- vec
    circos.track(
      ylim         = c(0, 1),
      bg.border    = NA,
      track.height = 0.041,
      panel.fun    = function(x, y) {
        pt  <- CELL_META$sector.index
        val <- as.character(vec_local[pt])
        col <- pal_local[val]
        if (length(col) == 0 || is.na(col)) col <- NA_COLOR
        circos.rect(0, 0, 1, 1, col = col, border = "white", lwd = 0.3)
      }
    )
  })

  if (i == 1) {
    circos.par(track.margin = c(0.002, 0.002))
  }
}

r_comp_bot <- get.cell.meta.data("cell.bottom.radius",
                                 sector.index = patient_order[1],
                                 track.index  = 2)
r_comp_top <- get.cell.meta.data("cell.top.radius",
                                 sector.index = patient_order[1],
                                 track.index  = 2)
r_comp_mid <- (r_comp_bot + r_comp_top) / 2

text(0, r_comp_top, "100%", cex = 0.72, font = 2, col = "grey25")
text(0, r_comp_mid, "50%",  cex = 0.68, font = 2, col = "grey45")
text(0, r_comp_bot, "0%",   cex = 0.72, font = 2, col = "grey25")

for (i in seq_along(clinical_cols)) {
  track_idx <- 2 + i
  r_top <- get.cell.meta.data("cell.top.radius",
                              sector.index = patient_order[1],
                              track.index  = track_idx)
  r_bot <- get.cell.meta.data("cell.bottom.radius",
                              sector.index = patient_order[1],
                              track.index  = track_idx)
  r <- (r_top + r_bot) / 2
  text(0.18, r, clinical_label_map[[clinical_cols[i]]],
       cex = 0.62, font = 2, col = "grey10", adj = c(0, 0.5))
}

text(0,  0.075, sprintf("%d patients", n_pt),
     cex = 1.45, font = 2, col = "grey5")
text(0, -0.015, sprintf("%s cells", formatC(nrow(md), big.mark = ",")),
     cex = 0.92, font = 2, col = "grey30")
text(0, -0.095, "4 malignant subtypes",
     cex = 0.78, font = 2, col = "grey40")

circos.clear()
upViewport()

build_legend <- function(labels, colors, title,
                         title_size = 11, label_size = 9,
                         grid_size  = 3.6) {
  Legend(
    labels      = labels,
    legend_gp   = gpar(fill = colors),
    title       = title,
    title_gp    = gpar(fontsize = title_size, fontface = "bold",
                       col = "grey5"),
    labels_gp   = gpar(fontsize = label_size, fontface = "bold",
                       col = "grey15"),
    grid_height = unit(grid_size, "mm"),
    grid_width  = unit(grid_size, "mm"),
    title_gap   = unit(1.6, "mm"),
    background  = "white",
    border      = FALSE
  )
}

lgd_subtype <- build_legend(
  labels     = subtype_levels,
  colors     = unname(subtype_colors[subtype_levels]),
  title      = "Malignant cell subtype",
  title_size = 12,
  label_size = 10,
  grid_size  = 4.2
)

dedup_pal <- function(pal) {
  if (all(c("unkown", "unknown") %in% names(pal)) &&
      pal[["unkown"]] == pal[["unknown"]]) {
    pal <- pal[names(pal) != "unkown"]
  }
  pal
}

clinical_lgds <- lapply(clinical_cols, function(v) {
  pal <- dedup_pal(clinical_palettes[[v]])
  build_legend(
    labels = names(pal),
    colors = unname(pal),
    title  = clinical_label_map[[v]]
  )
})

clinical_packed <- packLegend(
  list       = clinical_lgds,
  direction  = "vertical",
  row_gap    = unit(3.2, "mm"),
  column_gap = unit(5, "mm"),
  max_height = unit(8, "in")
)

all_lgd <- packLegend(
  lgd_subtype, clinical_packed,
  direction = "vertical",
  row_gap   = unit(7, "mm")
)

draw(all_lgd,
     x    = unit(9.4, "in"),
     y    = unit(0.5, "npc"),
     just = c("left", "center"))

dev.off()

n_all4 <- {
  tally <- comp %>%
    group_by(Pt_number) %>%
    summarise(present = sum(fraction > 0, na.rm = TRUE), .groups = "drop")
  sum(tally$present == 4)
}

sanity_tbl <- tibble::tibble(
  metric = c(
    "n_cells", "n_patients", "clinical_columns",
    "missing_clinical_values", "patient_order_rule",
    "patients_with_all_4_subtypes",
    "min_patient_cells", "max_patient_cells"
  ),
  value = c(
    nrow(md),
    length(patient_order),
    paste(clinical_cols, collapse = ";"),
    sum(is.na(md[, clinical_cols]) | md[, clinical_cols] == ""),
    "Proliferative-NPC_DESC_then_OPC-Myelination_ASC_then_Pt_number",
    n_all4,
    min(patient_summary$patient_total_cells),
    max(patient_summary$patient_total_cells)
  )
)
readr::write_csv(sanity_tbl, out_sanity)
writeLines(capture.output(sessionInfo()), out_session)

cat("F7 circular composition + clinical plot completed.\n")
cat("Patients:        ", length(patient_order), "\n")
cat("Output PDF:      ", out_pdf, "\n")
cat("Composition CSV: ", out_comp, "\n")
cat("Clinical CSV:    ", out_clinical, "\n")
cat("Patient order (Proliferative-NPC DESC, OPC-Myelination ASC):\n")
cat(paste(patient_order, collapse = ", "), "\n")
cat("Subtype totals:\n")
print(comp %>% group_by(subtype_label) %>%
        summarise(n = sum(n_cells), .groups = "drop"))

