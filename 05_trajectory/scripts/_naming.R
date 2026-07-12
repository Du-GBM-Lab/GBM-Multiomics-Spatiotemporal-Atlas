subtype_naming_mapping <- tibble::tribble(
  ~subtype_k4, ~subtype_label_original, ~subtype_label_final, ~abbreviation, ~color,
  "Subtype1", "Proliferative-NPC", "Proliferative-NPC-like", "NPC-P", "#0072B5",
  "Subtype2", "OPC-Myelination", "OPC-Myelination-like", "OPC-M", "#E18727",
  "Subtype3", "Vascular-niche MES", "Vascular-niche MES-like", "MES-V", "#20854E",
  "Subtype4", "MES-Antigen-presenting", "MES-Immune-interacting-like", "MES-I", "#BC3C29"
)

recode_subtype <- function(metadata) {
  stopifnot("subtype_k4" %in% colnames(metadata))
  if ("subtype_label_final" %in% colnames(metadata)) {
    metadata <- metadata |>
      dplyr::rename(subtype_label_source_metadata = "subtype_label_final")
  } else {
    metadata$subtype_label_source_metadata <- NA_character_
  }
  metadata |>
    dplyr::mutate(
      subtype_k4 = factor(as.character(.data$subtype_k4), levels = subtype_naming_mapping$subtype_k4)
    ) |>
    dplyr::left_join(
      subtype_naming_mapping |>
        dplyr::mutate(subtype_k4 = factor(.data$subtype_k4, levels = subtype_naming_mapping$subtype_k4)),
      by = "subtype_k4"
    ) |>
    dplyr::mutate(
      subtype_label_final = .data$subtype_label_final,
      subtype_short = factor(.data$abbreviation, levels = subtype_naming_mapping$abbreviation)
    )
}

write_subtype_naming_mapping <- function(path = file.path("tables", "subtype_naming_mapping.csv")) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(
    subtype_naming_mapping |>
      dplyr::select(subtype_k4, subtype_label_original, subtype_label_final, abbreviation),
    path
  )
}
