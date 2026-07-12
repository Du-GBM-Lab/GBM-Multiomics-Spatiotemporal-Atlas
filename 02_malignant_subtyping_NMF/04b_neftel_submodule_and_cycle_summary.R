# 05_恶性细胞分亚群与Neftel对照/04b_neftel_submodule_and_cycle_summary.R
# Summarize Neftel submodules separately, especially MES1/MES2 and cell-cycle axes.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(qs2)
  library(Seurat)
  library(dplyr)
  library(tidyr)
})

proj <- "05_恶性细胞分亚群与Neftel对照"
in_obj <- file.path(proj, "outputs", "GBM.malignant.subtyped.neftel_scored.qs2")
out_tab <- file.path(proj, "tables")
dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)

obj <- qs2::qs_read(in_obj)
md <- obj@meta.data

stopifnot("subtype_k4" %in% colnames(md))

required <- c(
  "AMS_MES1", "AMS_MES2", "AMS_AC", "AMS_OPC", "AMS_NPC1", "AMS_NPC2", "AMS_G1S", "AMS_G2M",
  "MES1_UCell", "MES2_UCell", "AC_UCell", "OPC_UCell", "NPC1_UCell", "NPC2_UCell", "G1S_UCell", "G2M_UCell"
)
missing <- setdiff(required, colnames(md))
if (length(missing) > 0) {
  stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

ams_modules <- c("AMS_MES1", "AMS_MES2", "AMS_AC", "AMS_OPC", "AMS_NPC1", "AMS_NPC2", "AMS_G1S", "AMS_G2M")
ucell_modules <- c("MES1_UCell", "MES2_UCell", "AC_UCell", "OPC_UCell", "NPC1_UCell", "NPC2_UCell", "G1S_UCell", "G2M_UCell")

summary_by_subtype <- md |>
  group_by(subtype_k4) |>
  summarise(
    n_cells = n(),
    across(all_of(ams_modules), list(mean = mean, median = median), .names = "{.col}_{.fn}"),
    across(all_of(ucell_modules), list(mean = mean, median = median), .names = "{.col}_{.fn}"),
    .groups = "drop"
  )

write.csv(
  summary_by_subtype,
  file.path(out_tab, "04b_neftel_submodule_scores_by_subtype.csv"),
  row.names = FALSE
)

# Dominant non-cycling submodule without merging MES1/MES2 or NPC1/NPC2.
ams_submodule_cols <- c("AMS_MES1", "AMS_MES2", "AMS_AC", "AMS_OPC", "AMS_NPC1", "AMS_NPC2")
ucell_submodule_cols <- c("MES1_UCell", "MES2_UCell", "AC_UCell", "OPC_UCell", "NPC1_UCell", "NPC2_UCell")

ams_names <- c("MES1", "MES2", "AC", "OPC", "NPC1", "NPC2")
ucell_names <- c("MES1", "MES2", "AC", "OPC", "NPC1", "NPC2")

md$neftel_submodule_AMS <- ams_names[apply(md[, ams_submodule_cols], 1, which.max)]
md$neftel_submodule_UCell <- ucell_names[apply(md[, ucell_submodule_cols], 1, which.max)]
md$neftel_mes_submodule_AMS <- ifelse(md$AMS_MES1 >= md$AMS_MES2, "MES1", "MES2")
md$neftel_mes_submodule_UCell <- ifelse(md$MES1_UCell >= md$MES2_UCell, "MES1", "MES2")

submodule_counts_ams <- as.data.frame.matrix(table(md$subtype_k4, md$neftel_submodule_AMS))
submodule_counts_ucell <- as.data.frame.matrix(table(md$subtype_k4, md$neftel_submodule_UCell))
submodule_pct_ams <- prop.table(as.matrix(table(md$subtype_k4, md$neftel_submodule_AMS)), margin = 1) * 100
submodule_pct_ucell <- prop.table(as.matrix(table(md$subtype_k4, md$neftel_submodule_UCell)), margin = 1) * 100

write.csv(submodule_counts_ams, file.path(out_tab, "04b_subtype_x_neftel_submodule_AMS_counts.csv"))
write.csv(submodule_counts_ucell, file.path(out_tab, "04b_subtype_x_neftel_submodule_UCell_counts.csv"))
write.csv(round(submodule_pct_ams, 2), file.path(out_tab, "04b_subtype_x_neftel_submodule_AMS_rowpct.csv"))
write.csv(round(submodule_pct_ucell, 2), file.path(out_tab, "04b_subtype_x_neftel_submodule_UCell_rowpct.csv"))

# MES1/MES2 balance inside each subtype.
mes_balance <- md |>
  group_by(subtype_k4) |>
  summarise(
    n_cells = n(),
    AMS_MES1_gt_MES2_pct = round(mean(AMS_MES1 >= AMS_MES2) * 100, 2),
    UCell_MES1_gt_MES2_pct = round(mean(MES1_UCell >= MES2_UCell) * 100, 2),
    AMS_MES1_mean = mean(AMS_MES1),
    AMS_MES2_mean = mean(AMS_MES2),
    UCell_MES1_mean = mean(MES1_UCell),
    UCell_MES2_mean = mean(MES2_UCell),
    .groups = "drop"
  )

write.csv(mes_balance, file.path(out_tab, "04b_MES1_MES2_balance_by_subtype.csv"), row.names = FALSE)

# Cycling: use both Neftel G1S/G2M and Seurat Phase if present.
cycle_summary <- md |>
  group_by(subtype_k4) |>
  summarise(
    n_cells = n(),
    AMS_G1S_mean = mean(AMS_G1S),
    AMS_G2M_mean = mean(AMS_G2M),
    AMS_cycling_mean = mean(AMS_cycling),
    UCell_G1S_mean = mean(G1S_UCell),
    UCell_G2M_mean = mean(G2M_UCell),
    UCell_cycling_mean = mean(UCell_cycling),
    .groups = "drop"
  )

write.csv(cycle_summary, file.path(out_tab, "04b_cycling_scores_by_subtype.csv"), row.names = FALSE)

if ("Phase" %in% colnames(md)) {
  phase_counts <- as.data.frame.matrix(table(md$subtype_k4, md$Phase))
  phase_pct <- prop.table(as.matrix(table(md$subtype_k4, md$Phase)), margin = 1) * 100
  write.csv(phase_counts, file.path(out_tab, "04b_subtype_x_seurat_phase_counts.csv"))
  write.csv(round(phase_pct, 2), file.path(out_tab, "04b_subtype_x_seurat_phase_rowpct.csv"))
} else {
  phase_pct <- NULL
}

# Persist the explicit submodule labels.
obj@meta.data$neftel_submodule_AMS <- md$neftel_submodule_AMS
obj@meta.data$neftel_submodule_UCell <- md$neftel_submodule_UCell
obj@meta.data$neftel_mes_submodule_AMS <- md$neftel_mes_submodule_AMS
obj@meta.data$neftel_mes_submodule_UCell <- md$neftel_mes_submodule_UCell
qs2::qs_save(obj, file.path(proj, "outputs", "GBM.malignant.subtyped.neftel_scored.submodule_labeled.qs2"))

cat("\nMES1/MES2 balance by subtype:\n")
print(mes_balance)

cat("\nDominant Neftel submodule, UCell row %:\n")
print(round(submodule_pct_ucell, 2))

cat("\nCycling scores by subtype:\n")
print(cycle_summary)

if (!is.null(phase_pct)) {
  cat("\nSeurat Phase row %:\n")
  print(round(phase_pct, 2))
}

cat("\nDone.\n")
