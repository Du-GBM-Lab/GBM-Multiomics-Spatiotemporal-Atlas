suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(stringr)
})

# Input:
#   outputs/GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2
#   tables/10e_one_vs_rest_FindAllMarkers_wilcoxon_downsampled.csv
#   tables/audit_manuscript_extracted_text.txt
#
# Output:
#   tables/manuscript_subtype_remapping.csv
#   tables/new_subtype_to_original_mapping.csv
#   tables/subtype4_manuscript_traceability.csv
#   tables/audit_manuscript_subtype_remapping_session_info.txt
#
# This audit does not change the subtype framework or figure outputs.
# It documents how the first-submission marker framework maps onto the
# revised high-confidence malignant-cell subtype framework.

step_dir <- "05_恶性细胞分亚群与Neftel对照"
if (basename(getwd()) != step_dir) {
  setwd(file.path(getwd(), step_dir))
}

obj_path <- file.path("outputs", "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2")
markers_path <- file.path("tables", "10e_one_vs_rest_FindAllMarkers_wilcoxon_downsampled.csv")
manuscript_text_path <- file.path("tables", "audit_manuscript_extracted_text.txt")

out_remap <- file.path("tables", "manuscript_subtype_remapping.csv")
out_reverse <- file.path("tables", "new_subtype_to_original_mapping.csv")
out_s4 <- file.path("tables", "subtype4_manuscript_traceability.csv")
out_session <- file.path("tables", "audit_manuscript_subtype_remapping_session_info.txt")

subtype_levels <- paste0("Subtype", 1:4)

original_markers <- tibble::tribble(
  ~original_subtype, ~original_marker, ~prior_expected_new_subtype,
  "Astro-sub", "GAP43", NA_character_,
  "Prolif-sub", "TOP2A", "Subtype1",
  "Mesen-sub", "COL1A1", "Subtype3",
  "Oligo-sub", "MAG", "Subtype2"
)

markers <- readr::read_csv(markers_path, show_col_types = FALSE)
stopifnot(all(c("subtype_k4", "gene", "avg_log2FC", "BH_q") %in% colnames(markers)))

markers_ranked <- markers %>%
  group_by(subtype_k4) %>%
  arrange(desc(avg_log2FC), BH_q, .by_group = TRUE) %>%
  mutate(rank_in_subtype = row_number()) %>%
  ungroup()

obj <- qs2::qs_read(obj_path)
DefaultAssay(obj) <- "RNA"
md <- obj@meta.data
stopifnot("subtype_k4" %in% colnames(md))
md$subtype_k4 <- factor(as.character(md$subtype_k4), levels = subtype_levels)
stopifnot(!any(is.na(md$subtype_k4)))

expr <- GetAssayData(obj, assay = "RNA", layer = "data")
genes <- original_markers$original_marker
stopifnot(all(genes %in% rownames(expr)))

mean_mat <- matrix(NA_real_, nrow = length(genes), ncol = length(subtype_levels),
                   dimnames = list(genes, subtype_levels))
pct_mat <- matrix(NA_real_, nrow = length(genes), ncol = length(subtype_levels),
                  dimnames = list(genes, subtype_levels))

for (g in genes) {
  vals <- expr[g, , drop = TRUE]
  for (s in subtype_levels) {
    cells <- which(md$subtype_k4 == s)
    vals_s <- vals[cells]
    mean_mat[g, s] <- Matrix::mean(vals_s)
    pct_mat[g, s] <- mean(vals_s > 0) * 100
  }
}

rank_lookup <- function(subtype, gene, col) {
  row <- markers_ranked %>% filter(subtype_k4 == subtype, gene == !!gene)
  if (nrow(row) == 0) return(NA)
  row[[col]][1]
}

classify_mapping <- function(original_subtype, gene, top_subtype, expected_subtype,
                             specificity_score, bh_q, pct_values, mean_values) {
  if (all(pct_values < 10)) return("lost")
  sorted_means <- sort(mean_values, decreasing = TRUE)
  if (length(sorted_means) >= 2 && sorted_means[2] > 0 && sorted_means[1] < 1.5 * sorted_means[2]) {
    return("split")
  }
  if (is.na(specificity_score) || specificity_score < 0.4) return("absorbed")
  if (!is.na(expected_subtype)) {
    if (!is.na(bh_q) && bh_q < 0.05 && top_subtype == expected_subtype) return("one_to_one")
    if (top_subtype != expected_subtype) return("remapped")
  }
  if (original_subtype == "Astro-sub" && top_subtype == "Subtype1") return("remapped")
  "remapped"
}

remap <- original_markers %>%
  rowwise() %>%
  mutate(
    mean_expr_S1 = mean_mat[original_marker, "Subtype1"],
    mean_expr_S2 = mean_mat[original_marker, "Subtype2"],
    mean_expr_S3 = mean_mat[original_marker, "Subtype3"],
    mean_expr_S4 = mean_mat[original_marker, "Subtype4"],
    pct_expressed_S1 = pct_mat[original_marker, "Subtype1"],
    pct_expressed_S2 = pct_mat[original_marker, "Subtype2"],
    pct_expressed_S3 = pct_mat[original_marker, "Subtype3"],
    pct_expressed_S4 = pct_mat[original_marker, "Subtype4"],
    new_top_subtype = subtype_levels[which.max(c(mean_expr_S1, mean_expr_S2, mean_expr_S3, mean_expr_S4))],
    specificity_score = max(c(mean_expr_S1, mean_expr_S2, mean_expr_S3, mean_expr_S4), na.rm = TRUE) /
      sum(c(mean_expr_S1, mean_expr_S2, mean_expr_S3, mean_expr_S4), na.rm = TRUE),
    new_top_subtype_rank = rank_lookup(new_top_subtype, original_marker, "rank_in_subtype"),
    log2FC_top_vs_rest = rank_lookup(new_top_subtype, original_marker, "avg_log2FC"),
    BH_q = rank_lookup(new_top_subtype, original_marker, "BH_q"),
    mapping_status = classify_mapping(
      original_subtype,
      original_marker,
      new_top_subtype,
      prior_expected_new_subtype,
      specificity_score,
      BH_q,
      c(pct_expressed_S1, pct_expressed_S2, pct_expressed_S3, pct_expressed_S4),
      c(mean_expr_S1, mean_expr_S2, mean_expr_S3, mean_expr_S4)
    ),
    narrative_note = case_when(
      original_marker == "GAP43" & new_top_subtype == "Subtype1" ~
        "GAP43 was an original Astro-sub representative marker, but in the revised framework it maps to Subtype1. Treat this as remapping/absorption of the original Astro-like signal into the revised Subtype1 rather than as evidence that an independent Astro-sub persists.",
      mapping_status == "one_to_one" ~
        paste0(original_subtype, " representative marker maps cleanly to ", new_top_subtype, "."),
      mapping_status == "split" ~
        "Marker is shared across multiple revised subtypes; avoid using it as a single-subtype-defining marker without qualification.",
      mapping_status == "absorbed" ~
        "Marker is no longer sufficiently subtype-specific under the revised framework.",
      TRUE ~
        "Marker is remapped under the revised framework; use narrative reconciliation."
    )
  ) %>%
  ungroup() %>%
  select(
    original_subtype, original_marker, new_top_subtype, new_top_subtype_rank,
    log2FC_top_vs_rest, BH_q,
    mean_expr_S1, mean_expr_S2, mean_expr_S3, mean_expr_S4,
    pct_expressed_S1, pct_expressed_S2, pct_expressed_S3, pct_expressed_S4,
    specificity_score, mapping_status, narrative_note
  )

readr::write_csv(remap, out_remap)

top30_long <- markers_ranked %>%
  filter(avg_log2FC > 0, BH_q < 0.05) %>%
  group_by(subtype_k4) %>%
  slice_head(n = 30) %>%
  ungroup()

reverse <- top30_long %>%
  group_by(subtype_k4) %>%
  summarise(
    top30_marker_genes = paste(gene, collapse = ";"),
    GAP43_rank = rank_in_subtype[match("GAP43", gene)],
    TOP2A_rank = rank_in_subtype[match("TOP2A", gene)],
    COL1A1_rank = rank_in_subtype[match("COL1A1", gene)],
    MAG_rank = rank_in_subtype[match("MAG", gene)],
    original_representative_markers_in_top30 = paste(intersect(gene, original_markers$original_marker), collapse = ";"),
    n_original_representative_markers_in_top30 = length(intersect(gene, original_markers$original_marker)),
    other_original_manuscript_molecules_in_top30 = paste(
      intersect(gene, c("PLAU", "PLAUR", "VIM", "CDH2", "CDH1", "MMP9", "MMP14", "FOSL1", "CD44", "FN1")),
      collapse = ";"
    ),
    .groups = "drop"
  ) %>%
  right_join(tibble(subtype_k4 = subtype_levels), by = "subtype_k4") %>%
  arrange(factor(subtype_k4, levels = subtype_levels)) %>%
  mutate(
    original_framework_interpretation = case_when(
      subtype_k4 == "Subtype1" ~ "Revised Subtype1 captures the proliferative/NPC-flavored program and is expected to contain TOP2A; GAP43 may be remapped here from the original Astro-sub framework.",
      subtype_k4 == "Subtype2" ~ "Revised Subtype2 maps to the original Oligo-sub/myelination axis if MAG or myelin markers are enriched.",
      subtype_k4 == "Subtype3" ~ "Revised Subtype3 maps to the original Mesen-sub ECM axis if COL1A1/ECM markers are enriched.",
      subtype_k4 == "Subtype4" ~ "Revised Subtype4 has no direct original representative wet-lab marker counterpart; interpret as a newly resolved MES antigen-presentation substate.",
      TRUE ~ NA_character_
    )
  )

readr::write_csv(reverse, out_reverse)

text_lines <- readLines(manuscript_text_path, encoding = "UTF-8", warn = FALSE)
trace_terms <- c(
  "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "B2M",
  "antigen presentation", "antigen-presenting", "MHC class II", "MHC-II",
  "immune-like", "immunoreactive"
)

traceability <- lapply(trace_terms, function(term) {
  hit_idx <- which(str_detect(text_lines, regex(term, ignore_case = TRUE)))
  if (length(hit_idx) == 0) {
    return(tibble(
      keyword = term,
      mentioned_in_original_manuscript = FALSE,
      paragraph_number = NA_integer_,
      context_excerpt = NA_character_
    ))
  }
  lapply(hit_idx, function(i) {
    line <- text_lines[i]
    paragraph_number <- suppressWarnings(as.integer(str_match(line, "\\[word/document\\.xml:(\\d+)\\]")[, 2]))
    plain <- str_replace(line, "^\\[[^\\]]+\\]\\t", "")
    loc <- str_locate(str_to_lower(plain), fixed(str_to_lower(term)))[1, ]
    if (any(is.na(loc))) {
      context <- substr(plain, 1, min(nchar(plain), 180))
    } else {
      start <- max(1, loc[1] - 50)
      end <- min(nchar(plain), loc[2] + 50)
      context <- substr(plain, start, end)
    }
    tibble(
      keyword = term,
      mentioned_in_original_manuscript = TRUE,
      paragraph_number = paragraph_number,
      context_excerpt = context
    )
  }) %>% bind_rows()
}) %>% bind_rows()

readr::write_csv(traceability, out_s4)
writeLines(capture.output(sessionInfo()), out_session)

cat("manuscript_subtype_remapping.csv:\n")
print(remap %>% select(original_subtype, original_marker, new_top_subtype, new_top_subtype_rank, log2FC_top_vs_rest, BH_q, specificity_score, mapping_status))
cat("\nnew_subtype_to_original_mapping.csv:\n")
print(reverse %>% select(subtype_k4, original_representative_markers_in_top30, n_original_representative_markers_in_top30, other_original_manuscript_molecules_in_top30))
cat("\nsubtype4_manuscript_traceability.csv hits:\n")
print(traceability %>% filter(mentioned_in_original_manuscript))
cat("\nOutputs:\n")
cat(out_remap, "\n")
cat(out_reverse, "\n")
cat(out_s4, "\n")
