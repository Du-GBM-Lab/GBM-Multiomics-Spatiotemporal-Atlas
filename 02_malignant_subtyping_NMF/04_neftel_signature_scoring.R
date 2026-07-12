# 05_恶性细胞分亚群与Neftel对照/04_neftel_signature_scoring.R
# Neftel 2019 meta-module scoring on malignant cells.
# Gene sets are read from the user-provided official Table S2 Excel file.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(Seurat)
  library(qs2)
  library(readxl)
  library(UCell)
  library(dplyr)
  library(tidyr)
})

set.seed(42)

proj <- "05_恶性细胞分亚群与Neftel对照"
signature_xlsx <- "<DATA_ROOT>/zetora/storage/7SGY8BYX/亚型marker.xlsx"

in_obj <- file.path(proj, "outputs", "GBM.malignant.subtyped.theta5.k4.qs2")
out_obj <- file.path(proj, "outputs", "GBM.malignant.subtyped.neftel_scored.qs2")
out_tab <- file.path(proj, "tables")

dir.create(out_tab, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

read_neftel_signatures <- function(path) {
  if (!file.exists(path)) {
    stop("Neftel signature Excel not found: ", path, call. = FALSE)
  }

  raw <- readxl::read_excel(path, sheet = "Table S2", skip = 3, col_types = "text")
  names(raw) <- gsub("\\s+", "", names(raw))
  names(raw) <- gsub("/", "", names(raw))

  expected <- c("MES2", "MES1", "AC", "OPC", "NPC1", "NPC2", "G1S", "G2M")
  missing <- setdiff(expected, names(raw))
  if (length(missing) > 0) {
    stop("Missing signature columns in Excel: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  sigs <- lapply(expected, function(module) {
    genes <- raw[[module]]
    genes <- trimws(as.character(genes))
    genes <- genes[!is.na(genes) & genes != ""]
    unique(genes)
  })
  names(sigs) <- expected
  sigs[c("MES1", "MES2", "AC", "OPC", "NPC1", "NPC2", "G1S", "G2M")]
}

safe_col <- function(md, candidates) {
  hit <- candidates[candidates %in% colnames(md)]
  if (length(hit) == 0) {
    stop("None of these metadata columns exist: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  hit[1]
}

msg("Reading Neftel signatures from:", signature_xlsx)
neftel_sigs_raw <- read_neftel_signatures(signature_xlsx)

raw_gene_count <- data.frame(
  module = names(neftel_sigs_raw),
  n_genes_in_excel = vapply(neftel_sigs_raw, length, integer(1)),
  stringsAsFactors = FALSE
)

msg("Loading object:", in_obj)
obj <- qs2::qs_read(in_obj)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

if (!"data" %in% Layers(obj[["RNA"]])) {
  msg("RNA data layer not found; running NormalizeData.")
  obj <- NormalizeData(obj, verbose = FALSE)
}

stopifnot("subtype_k4" %in% colnames(obj@meta.data))

neftel_sigs <- lapply(neftel_sigs_raw, function(g) intersect(g, rownames(obj)))
gene_used <- raw_gene_count |>
  mutate(
    n_genes_used = vapply(neftel_sigs[module], length, integer(1)),
    recovery_pct = round(100 * n_genes_used / n_genes_in_excel, 1)
  )

msg("Genes recovered per module:")
print(gene_used)
if (any(gene_used$n_genes_used < 25)) {
  stop("At least one Neftel module has fewer than 25 detected genes. Check gene symbols.", call. = FALSE)
}

gene_detail <- bind_rows(lapply(names(neftel_sigs_raw), function(module) {
  data.frame(
    module = module,
    gene = neftel_sigs_raw[[module]],
    used_in_object = neftel_sigs_raw[[module]] %in% rownames(obj),
    stringsAsFactors = FALSE
  )
}))

write.csv(gene_used, file.path(out_tab, "04_neftel_module_genes_used.csv"), row.names = FALSE)
write.csv(gene_detail, file.path(out_tab, "04_neftel_module_gene_detail.csv"), row.names = FALSE)

msg("Running AddModuleScore.")
obj <- AddModuleScore(
  obj,
  features = neftel_sigs,
  name = "AMS_tmp_",
  ctrl = 100,
  seed = 42
)

for (i in seq_along(neftel_sigs)) {
  obj@meta.data[[paste0("AMS_", names(neftel_sigs)[i])]] <- obj@meta.data[[paste0("AMS_tmp_", i)]]
  obj@meta.data[[paste0("AMS_tmp_", i)]] <- NULL
}

msg("Running UCell.")
obj <- AddModuleScore_UCell(
  obj,
  features = neftel_sigs,
  ncores = 4
)

md <- obj@meta.data

ucell_cols_expected <- paste0(names(neftel_sigs), "_UCell")
missing_ucell <- setdiff(ucell_cols_expected, colnames(md))
if (length(missing_ucell) > 0) {
  stop("UCell output columns missing: ", paste(missing_ucell, collapse = ", "), call. = FALSE)
}

md$AMS_MES <- pmax(md$AMS_MES1, md$AMS_MES2)
md$AMS_NPC <- pmax(md$AMS_NPC1, md$AMS_NPC2)
md$AMS_cycling <- md$AMS_G1S + md$AMS_G2M

md$UCell_MES <- pmax(md$MES1_UCell, md$MES2_UCell)
md$UCell_AC <- md$AC_UCell
md$UCell_OPC <- md$OPC_UCell
md$UCell_NPC <- pmax(md$NPC1_UCell, md$NPC2_UCell)
md$UCell_cycling <- md$G1S_UCell + md$G2M_UCell

states <- c("MES", "AC", "OPC", "NPC")
md$neftel_state_AMS <- states[apply(md[, c("AMS_MES", "AMS_AC", "AMS_OPC", "AMS_NPC")], 1, which.max)]
md$neftel_state_UCell <- states[apply(md[, c("UCell_MES", "UCell_AC", "UCell_OPC", "UCell_NPC")], 1, which.max)]
md$neftel_state_concordant <- md$neftel_state_AMS == md$neftel_state_UCell

obj@meta.data <- md

msg("Saving scored object:", out_obj)
qs2::qs_save(obj, out_obj)

state_dist_ams <- as.data.frame(table(md$neftel_state_AMS))
colnames(state_dist_ams) <- c("state", "n_cells")
state_dist_ams$pct <- round(100 * state_dist_ams$n_cells / sum(state_dist_ams$n_cells), 2)

state_dist_ucell <- as.data.frame(table(md$neftel_state_UCell))
colnames(state_dist_ucell) <- c("state", "n_cells")
state_dist_ucell$pct <- round(100 * state_dist_ucell$n_cells / sum(state_dist_ucell$n_cells), 2)

concordance <- mean(md$neftel_state_concordant) * 100
cor_check <- sapply(states, function(s) {
  cor(
    md[[paste0("AMS_", s)]],
    md[[paste0("UCell_", s)]],
    method = "spearman",
    use = "complete.obs"
  )
})
cor_check_df <- data.frame(state = names(cor_check), spearman = as.numeric(cor_check), row.names = NULL)

subtype_x_ams <- as.data.frame.matrix(table(md$subtype_k4, md$neftel_state_AMS))
subtype_x_ucell <- as.data.frame.matrix(table(md$subtype_k4, md$neftel_state_UCell))

subtype_x_ams_pct <- prop.table(as.matrix(table(md$subtype_k4, md$neftel_state_AMS)), margin = 1) * 100
subtype_x_ucell_pct <- prop.table(as.matrix(table(md$subtype_k4, md$neftel_state_UCell)), margin = 1) * 100

write.csv(state_dist_ams, file.path(out_tab, "04_neftel_state_dist_AMS.csv"), row.names = FALSE)
write.csv(state_dist_ucell, file.path(out_tab, "04_neftel_state_dist_UCell.csv"), row.names = FALSE)
write.csv(subtype_x_ams, file.path(out_tab, "04_subtype_x_state_AMS_counts.csv"))
write.csv(subtype_x_ucell, file.path(out_tab, "04_subtype_x_state_UCell_counts.csv"))
write.csv(round(subtype_x_ams_pct, 2), file.path(out_tab, "04_subtype_x_state_AMS_rowpct.csv"))
write.csv(round(subtype_x_ucell_pct, 2), file.path(out_tab, "04_subtype_x_state_UCell_rowpct.csv"))
write.csv(cor_check_df, file.path(out_tab, "04_AMS_UCell_spearman.csv"), row.names = FALSE)

summary_df <- data.frame(
  metric = c("n_cells", "AMS_UCell_state_concordance_pct"),
  value = c(ncol(obj), round(concordance, 2))
)
write.csv(summary_df, file.path(out_tab, "04_neftel_scoring_summary.csv"), row.names = FALSE)

cat("\nCells:", ncol(obj), "\n")
cat("\nAMS state distribution:\n")
print(state_dist_ams)
cat("\nUCell state distribution:\n")
print(state_dist_ucell)
cat("\nAMS-UCell concordance:", round(concordance, 1), "%\n")
cat("\nAMS-UCell Spearman per state:\n")
print(round(cor_check, 3))
cat("\nSubtype x Neftel state, AMS row %:\n")
print(round(subtype_x_ams_pct, 2))
cat("\nSubtype x Neftel state, UCell row %:\n")
print(round(subtype_x_ucell_pct, 2))

msg("Done.")
