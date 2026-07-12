# R4 04: build new four-subtype gene sets from the reanalysis marker table.
# Output is used for bulk CGGA ssGSEA audit. No survival or clustering here.

suppressPackageStartupMessages({
  library(dplyr)
})

root <- getwd()
marker_file <- file.path(
  dirname(root),
  "05_恶性细胞分亚群与Neftel对照",
  "tables",
  "10e_one_vs_rest_FindAllMarkers_wilcoxon_downsampled.csv"
)
out_dir <- file.path(root, "data", "processed")
tab_dir <- file.path(root, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

SUBTYPE_LEVELS <- c("NPC_P", "OPC_M", "MES_V", "MES_I")
SUBTYPE_MAP <- c(
  "Subtype1" = "NPC_P",
  "Subtype2" = "OPC_M",
  "Subtype3" = "MES_V",
  "Subtype4" = "MES_I"
)
FORBIDDEN <- c("PLAUR", "FOSL1", "FOSL2")

PADJ_MAX <- 0.05
LFC_MIN <- 0.50
PCT1_MIN <- 0.25
DELTA_PCT <- 0.10
TOP_N <- 100
DROP_SHARED <- TRUE

cat("marker_file:", marker_file, "\n")
stopifnot(file.exists(marker_file))

fam <- read.csv(marker_file, check.names = FALSE)
required <- c("gene", "cluster", "avg_log2FC", "p_val_adj", "pct.1", "pct.2")
missing_required <- setdiff(required, colnames(fam))
cat("missing required columns:", paste(missing_required, collapse = ", "), "\n")
stopifnot(length(missing_required) == 0)
stopifnot(all(names(SUBTYPE_MAP) %in% unique(fam$cluster)))

fam$subtype_r4 <- unname(SUBTYPE_MAP[as.character(fam$cluster)])

sel <- fam %>%
  filter(
    subtype_r4 %in% SUBTYPE_LEVELS,
    p_val_adj < PADJ_MAX,
    avg_log2FC >= LFC_MIN,
    pct.1 >= PCT1_MIN,
    (pct.1 - pct.2) >= DELTA_PCT,
    !toupper(gene) %in% FORBIDDEN
  ) %>%
  group_by(subtype_r4) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  distinct(gene, .keep_all = TRUE) %>%
  slice_head(n = TOP_N) %>%
  ungroup()

shared <- names(which(table(sel$gene) > 1))
cat("cross-subtype shared genes before dropping:", length(shared), "\n")
if (DROP_SHARED && length(shared)) {
  sel <- sel %>% filter(!gene %in% shared)
}

gene_sets <- split(sel$gene, factor(sel$subtype_r4, levels = SUBTYPE_LEVELS))
gene_sets <- lapply(gene_sets, function(x) unique(as.character(x[!is.na(x) & x != ""])))

forbidden_hits <- intersect(toupper(unlist(gene_sets)), FORBIDDEN)
stopifnot(length(forbidden_hits) == 0)

size_tbl <- data.frame(
  subtype = names(gene_sets),
  n_genes = lengths(gene_sets),
  stringsAsFactors = FALSE
)
param_tbl <- data.frame(
  parameter = c("PADJ_MAX", "LFC_MIN", "PCT1_MIN", "DELTA_PCT", "TOP_N", "DROP_SHARED"),
  value = c(PADJ_MAX, LFC_MIN, PCT1_MIN, DELTA_PCT, TOP_N, DROP_SHARED),
  stringsAsFactors = FALSE
)

cat("\n== gene set sizes ==\n")
print(size_tbl)
cat("\n== first genes ==\n")
print(lapply(gene_sets, head, 12))

saveRDS(gene_sets, file.path(out_dir, "sc_subtype_markers.rds"))
write.csv(sel, file.path(tab_dir, "STOP2_四亚型gene_set入选基因.csv"), row.names = FALSE)
write.csv(size_tbl, file.path(tab_dir, "STOP2_四亚型gene_set大小.csv"), row.names = FALSE)
write.csv(param_tbl, file.path(tab_dir, "STOP2_四亚型gene_set参数.csv"), row.names = FALSE)

cat("\nWrote:", file.path(out_dir, "sc_subtype_markers.rds"), "\n")
cat("R4 boundary checked: PLAUR/FOSL1/FOSL2 absent.\n")
