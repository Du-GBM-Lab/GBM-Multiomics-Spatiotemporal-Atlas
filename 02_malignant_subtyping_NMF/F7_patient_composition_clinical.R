suppressPackageStartupMessages({
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir) {
  setwd(file.path(getwd(), step_dir))
}

obj_path <- "outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"
out_pdf <- "figures/F7_patient_composition_clinical.pdf"
out_comp <- "figures/source_data/F7_patient_subtype_composition.csv"
out_clinical <- "figures/source_data/F7_patient_clinical_annotation.csv"
out_sanity <- "figures/source_data/F7_patient_composition_clinical_sanity_checks.csv"
out_session <- "figures/source_data/F7_patient_composition_clinical_session_info.txt"

dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/source_data", showWarnings = FALSE, recursive = TRUE)

subtype_levels <- paste0("Subtype", 1:4)
subtype_colors <- c(
  "Subtype1" = "#0072B5",
  "Subtype2" = "#E18727",
  "Subtype3" = "#20854E",
  "Subtype4" = "#BC3C29"
)
subtype_legend_labels <- c(
  "Subtype1" = "Subtype1: Proliferative-NPC",
  "Subtype2" = "Subtype2: OPC-Myelination",
  "Subtype3" = "Subtype3: Vascular-niche MES",
  "Subtype4" = "Subtype4: MES-Antigen-presenting"
)

clinical_cols <- c("diagnosis", "WHO", "gender", "IDH", "age", "treatment_1", "treatment_3")
clinical_label_map <- c(
  "diagnosis" = "Diagnosis",
  "WHO" = "WHO",
  "gender" = "Gender",
  "IDH" = "IDH",
  "age" = "Age",
  "treatment_1" = "Treatment 1",
  "treatment_3" = "Treatment 3"
)

clinical_palettes <- list(
  diagnosis = c("recurrence" = "#7B3294", "ND" = "#A6DBA0"),
  WHO = c("four" = "#D6604D", "three" = "#4393C3"),
  gender = c("male" = "#4D4D4D", "female" = "#BDBDBD"),
  IDH = c("WT" = "#B2182B", "mutant" = "#2166AC"),
  age = c("below60" = "#92C5DE", "above60" = "#F4A582", "unkown" = "#D9D9D9", "unknown" = "#D9D9D9"),
  treatment_1 = c("TMZ" = "#1B9E77", "untreated" = "#D9D9D9", "PD1" = "#7570B3"),
  treatment_3 = c("Rec" = "#7B3294", "Untreated" = "#D9D9D9", "nonresponder" = "#E7298A", "responder" = "#66A61E")
)

obj <- qs2::qs_read(obj_path)
md <- obj@meta.data
required_cols <- c("Pt_number", "subtype_k4", clinical_cols)
missing_cols <- setdiff(required_cols, colnames(md))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

md <- md %>%
  mutate(
    Pt_number = as.character(Pt_number),
    subtype_k4 = factor(as.character(subtype_k4), levels = subtype_levels)
  )

comp <- md %>%
  count(Pt_number, subtype_k4, name = "n_cells") %>%
  complete(Pt_number, subtype_k4 = factor(subtype_levels, levels = subtype_levels), fill = list(n_cells = 0)) %>%
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
    dominant_subtype = subtype_k4[which.max(fraction)],
    dominant_fraction = max(fraction),
    Subtype1_fraction = fraction[subtype_k4 == "Subtype1"],
    Subtype2_fraction = fraction[subtype_k4 == "Subtype2"],
    Subtype3_fraction = fraction[subtype_k4 == "Subtype3"],
    Subtype4_fraction = fraction[subtype_k4 == "Subtype4"],
    .groups = "drop"
  )

patient_order <- patient_summary %>%
  arrange(dominant_subtype, desc(dominant_fraction), desc(patient_total_cells), Pt_number) %>%
  pull(Pt_number)

comp <- comp %>%
  mutate(
    Pt_number = factor(Pt_number, levels = patient_order),
    subtype_k4 = factor(subtype_k4, levels = subtype_levels)
  )

clinical <- md %>%
  select(Pt_number, all_of(clinical_cols)) %>%
  distinct() %>%
  group_by(Pt_number) %>%
  summarise(
    across(all_of(clinical_cols), ~ paste(unique(as.character(.x)), collapse = " | ")),
    .groups = "drop"
  ) %>%
  mutate(Pt_number = factor(Pt_number, levels = patient_order)) %>%
  arrange(Pt_number)

readr::write_csv(comp %>% mutate(Pt_number = as.character(Pt_number)), out_comp)
readr::write_csv(clinical %>% mutate(Pt_number = as.character(Pt_number)), out_clinical)

clinical_long <- clinical %>%
  pivot_longer(cols = all_of(clinical_cols), names_to = "clinical_variable", values_to = "clinical_value") %>%
  mutate(
    clinical_label = clinical_label_map[clinical_variable],
    clinical_variable = factor(clinical_variable, levels = rev(clinical_cols), labels = rev(clinical_label_map[clinical_cols])),
    clinical_value = as.character(clinical_value),
    clinical_key = paste0(clinical_label, ": ", clinical_value)
  )

clinical_colors <- unlist(lapply(names(clinical_palettes), function(v) {
  vals <- clinical_palettes[[v]]
  stats::setNames(unname(vals), paste0(clinical_label_map[[v]], ": ", names(vals)))
}))
extra_values <- setdiff(unique(clinical_long$clinical_key), names(clinical_colors))
fallback_cols <- setNames(rep("#F0F0F0", length(extra_values)), extra_values)
clinical_colors <- c(clinical_colors, fallback_cols)

p_ann <- ggplot(clinical_long, aes(x = Pt_number, y = clinical_variable, fill = clinical_key)) +
  geom_tile(color = "white", linewidth = 0.25, height = 0.9) +
  scale_fill_manual(values = clinical_colors, name = NULL, na.value = "#F0F0F0") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 7, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 6.5),
    legend.key.size = unit(2.8, "mm"),
    plot.margin = margin(t = 2, r = 4, b = 0, l = 2)
  ) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(linewidth = 0)))

p_bar <- ggplot(comp, aes(x = Pt_number, y = fraction, fill = subtype_k4)) +
  geom_col(width = 0.82, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = subtype_colors, labels = subtype_legend_labels, name = "Subtype") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.01))) +
  labs(x = NULL, y = "Fraction of malignant cells") +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.5, color = "black"),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.title.y = element_text(size = 8),
    axis.line = element_line(linewidth = 0.35, color = "grey40"),
    axis.ticks = element_line(linewidth = 0.3, color = "grey40"),
    legend.position = "bottom",
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, face = "bold"),
    legend.key.size = unit(3, "mm"),
    plot.margin = margin(t = 1, r = 4, b = 2, l = 2)
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

combined <- p_ann / p_bar + plot_layout(heights = c(0.43, 1), guides = "collect") &
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(out_pdf, combined, width = 9.2, height = 5.2, useDingbats = FALSE)

sanity_tbl <- tibble::tibble(
  metric = c(
    "n_cells",
    "n_patients",
    "clinical_columns",
    "missing_clinical_values",
    "patient_order_rule",
    "patients_with_all_4_subtypes",
    "min_patient_cells",
    "max_patient_cells"
  ),
  value = c(
    nrow(md),
    length(patient_order),
    paste(clinical_cols, collapse = ";"),
    sum(is.na(md[, clinical_cols]) | md[, clinical_cols] == ""),
    "dominant_subtype_then_dominant_fraction_then_cell_count",
    sum(patient_summary %>% rowwise() %>% mutate(n_present = sum(c_across(ends_with("_fraction")) > 0)) %>% pull(n_present) == 4),
    min(patient_summary$patient_total_cells),
    max(patient_summary$patient_total_cells)
  )
)
readr::write_csv(sanity_tbl, out_sanity)
writeLines(capture.output(sessionInfo()), out_session)

cat("F7 patient composition + clinical annotation completed.\n")
cat("Patients:", length(patient_order), "\n")
cat("Patient order:\n")
cat(paste(patient_order, collapse = ", "), "\n")
cat("Clinical columns:", paste(clinical_cols, collapse = ", "), "\n")
cat("Subtype cell totals:\n")
print(comp %>% group_by(subtype_k4) %>% summarise(n = sum(n_cells), .groups = "drop"))
cat("Output PDF:", out_pdf, "\n")
cat("Composition source:", out_comp, "\n")
cat("Clinical source:", out_clinical, "\n")
