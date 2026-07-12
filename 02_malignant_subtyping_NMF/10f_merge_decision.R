# 05_恶性细胞分亚群与Neftel对照/10f_merge_decision.R
# Rule v2/v2.2 merge decision for Subtype3 vs Subtype4.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(dplyr)
  library(readr)
})

params <- list(
  tables_dir = file.path("05_恶性细胞分亚群与Neftel对照", "tables"),
  evidence1_file = file.path("05_恶性细胞分亚群与Neftel对照", "tables", "10c_evidence1_decision.csv"),
  evidence2_file = file.path("05_恶性细胞分亚群与Neftel对照", "tables", "10d_evidence2_decision.csv"),
  evidence3_file = file.path("05_恶性细胞分亚群与Neftel对照", "tables", "10e_evidence3_decision.csv")
)

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  if (is.na(x) || x == "") return(NA)
  toupper(as.character(x)) %in% c("TRUE", "T", "YES", "Y", "1")
}

e1 <- read_csv(params$evidence1_file, show_col_types = FALSE)
e2 <- read_csv(params$evidence2_file, show_col_types = FALSE)
e3 <- read_csv(params$evidence3_file, show_col_types = FALSE)

evidence1_satisfied <- as_logical_safe(e1$evidence_satisfied[1])
evidence2_statistical <- as_logical_safe(e2$statistical_threshold_passed[1])
evidence2_visual <- as_logical_safe(e2$visually_separable[1])
evidence2_satisfied <- if (is.na(evidence2_visual)) {
  NA
} else {
  evidence2_statistical && evidence2_visual
}
evidence3_satisfied <- as_logical_safe(e3$evidence3_satisfied[1])

n_satisfied_confirmed <- sum(c(evidence1_satisfied, evidence2_satisfied, evidence3_satisfied), na.rm = TRUE)
n_pending <- sum(is.na(c(evidence1_satisfied, evidence2_satisfied, evidence3_satisfied)))

decision_status <- if (n_pending > 0) {
  "pending_visual_heatmap_review"
} else if (n_satisfied_confirmed >= 2) {
  "keep_S3_S4_separate"
} else if (n_satisfied_confirmed == 1) {
  "keep_separate_with_limitations"
} else {
  "merge_S3_S4_k3"
}

decision_if_visual_true <- {
  e2_if <- evidence2_statistical && TRUE
  n_if <- sum(c(evidence1_satisfied, e2_if, evidence3_satisfied), na.rm = TRUE)
  if (n_if >= 2) {
    "keep_S3_S4_separate_k4_final"
  } else if (n_if == 1) {
    "keep_separate_with_limitations"
  } else {
    "merge_S3_S4_k3"
  }
}

decision_if_visual_false <- {
  e2_if <- FALSE
  n_if <- sum(c(evidence1_satisfied, e2_if, evidence3_satisfied), na.rm = TRUE)
  if (n_if >= 2) {
    "keep_S3_S4_separate_k4_final"
  } else if (n_if == 1) {
    "keep_separate_with_limitations"
  } else {
    "merge_S3_S4_k3"
  }
}

decision <- tibble(
  rule_version = "Rule_v2.2",
  evidence1_NMF_program_level = evidence1_satisfied,
  evidence1_detail = paste0(
    "PERMANOVA_p=", e1$PERMANOVA_p[1],
    "; passing_metaprograms=", e1$metaprograms_passing_univariate[1]
  ),
  evidence2_marker_DE_statistical = evidence2_statistical,
  evidence2_marker_heatmap_visually_separable = evidence2_visual,
  evidence2_marker_DE_satisfied = evidence2_satisfied,
  evidence2_detail = paste0(
    "n_passing=", e2$n_passing[1],
    "; n_S3_up=", e2$n_S3_up[1],
    "; n_S4_up=", e2$n_S4_up[1]
  ),
  evidence3_pathway_functional = evidence3_satisfied,
  evidence3_detail = paste0(
    "S3_sig=", e3$n_S3_enriched_sig_pathways[1],
    "; S4_sig=", e3$n_S4_enriched_sig_pathways[1]
  ),
  n_satisfied_confirmed = n_satisfied_confirmed,
  n_pending = n_pending,
  decision_status = decision_status,
  decision_if_heatmap_visually_separable_TRUE = decision_if_visual_true,
  decision_if_heatmap_visually_separable_FALSE = decision_if_visual_false,
  recommended_next_action = ifelse(
    is.na(evidence2_visual),
    "Manually review 10d_S3_vs_S4_marker_heatmap_draft.pdf and set visually_separable TRUE/FALSE in 10d_evidence2_decision.csv.",
    "Proceed to final figure design."
  )
)

write_csv(decision, file.path(params$tables_dir, "10f_merge_decision.csv"))

subtype_candidates <- tibble::tribble(
  ~subtype_k4, ~subtype_label_candidate, ~label_status, ~basis,
  "Subtype1", "Proliferative-NPC subtype", "candidate_pending_heatmap_10f_finalization", "Cycling axis; NPC2 mixed signal; downsampled GSEA G2M/E2F/MYC/DNA repair.",
  "Subtype2", "OPC-Myelination subtype", "candidate_pending_heatmap_10f_finalization", "MP03 myelination program; OPC signal; neuronal-system/myelination-related pathway overview.",
  "Subtype3", "MES-ECM subtype", "candidate_pending_heatmap_10f_finalization", "MP02 ECM program and MP06 angiogenesis-signaling program; S3 GSEA EMT/OXPHOS/myogenesis/hypoxia.",
  "Subtype4", "MES-Antigen-presenting subtype", "candidate_pending_heatmap_10f_finalization", "MP04 MHC-II antigen-presentation program; S4 GSEA TNFA/IFN/allograft/complement."
)
write_csv(subtype_candidates, file.path(params$tables_dir, "10f_subtype_label_candidates.csv"))

writeLines(capture.output(sessionInfo()), file.path(params$tables_dir, "10f_session_info.txt"))
print(decision)
