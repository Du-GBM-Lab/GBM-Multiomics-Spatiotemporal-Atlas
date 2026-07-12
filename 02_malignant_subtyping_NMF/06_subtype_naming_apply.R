# 05_恶性细胞分亚群与Neftel对照/06_subtype_naming_apply.R
# Apply biological subtype names and export final UMAP panels.
# This only adds label metadata; it does not change clusters or subtype assignment.

suppressPackageStartupMessages({
  .libPaths(c("<DATA_ROOT>/环境/稳稳的r包", .libPaths()))
  library(Seurat)
  library(qs2)
  library(ggplot2)
  library(dplyr)
})

set.seed(42)

proj <- "05_恶性细胞分亚群与Neftel对照"
in_obj <- file.path(proj, "outputs", "GBM.malignant.subtyped.umap_candidates.qs2")
out_obj <- file.path(proj, "outputs", "GBM.malignant.subtyped.named.qs2")
fig_dir <- file.path(proj, "figures")
tab_dir <- file.path(proj, "tables")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

name_map <- c(
  "Subtype1" = "NPC-Cycling",
  "Subtype2" = "OPC-like",
  "Subtype3" = "MES-Perivascular",
  "Subtype4" = "MES-Inflammatory"
)

named_levels <- c("NPC-Cycling", "OPC-like", "MES-Perivascular", "MES-Inflammatory")
named_palette <- c(
  "NPC-Cycling" = "#00468B",
  "OPC-like" = "#ED0000",
  "MES-Perivascular" = "#42B540",
  "MES-Inflammatory" = "#0099B4"
)
named_palette_alpha <- grDevices::adjustcolor(named_palette, alpha.f = 0.62)

msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

theme_umap_final <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.key.size = unit(0.42, "cm"),
      plot.title = element_text(size = 12, face = "bold"),
      axis.title = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.4)
    )
}

msg("Loading object:", in_obj)
obj <- qs2::qs_read(in_obj)
DefaultAssay(obj) <- "RNA"
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

stopifnot("subtype_k4" %in% colnames(obj@meta.data))
stopifnot("umap_closer4" %in% Reductions(obj))

obj$subtype_k4 <- factor(as.character(obj$subtype_k4), levels = names(name_map))
obj$subtype_named <- factor(
  unname(name_map[as.character(obj$subtype_k4)]),
  levels = named_levels
)

if (any(is.na(obj$subtype_named))) {
  stop("Subtype naming generated NA values. Check subtype_k4 labels.", call. = FALSE)
}

rename_tab <- table(subtype_k4 = obj$subtype_k4, subtype_named = obj$subtype_named)
expected_nonzero <- cbind(names(name_map), unname(name_map))
for (i in seq_len(nrow(expected_nonzero))) {
  old <- expected_nonzero[i, 1]
  new <- expected_nonzero[i, 2]
  if (rename_tab[old, new] != sum(obj$subtype_k4 == old)) {
    stop("Rename sanity failed for ", old, " -> ", new, call. = FALSE)
  }
}
if (sum(rename_tab[row(rename_tab) != match(colnames(rename_tab), unname(name_map))[col(rename_tab)]]) != 0) {
  stop("Rename sanity failed: off-diagonal mappings detected.", call. = FALSE)
}

cat("\n== Rename sanity ==\n")
print(rename_tab)

p_main <- DimPlot(
  obj,
  reduction = "umap_closer4",
  group.by = "subtype_named",
  pt.size = 0.30,
  label = TRUE,
  label.size = 5,
  repel = TRUE
) +
  scale_color_manual(values = named_palette_alpha, drop = FALSE) +
  ggtitle("Malignant subtypes (n=28,213 cells, 24 patients)") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_umap_final()

ggsave(
  file.path(fig_dir, "06_final_main_UMAP_named.pdf"),
  p_main,
  width = 7.5,
  height = 5.5,
  device = cairo_pdf
)

p_nolabel <- DimPlot(
  obj,
  reduction = "umap_closer4",
  group.by = "subtype_named",
  pt.size = 0.30,
  label = FALSE
) +
  scale_color_manual(values = named_palette_alpha, drop = FALSE) +
  ggtitle("Malignant subtypes") +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme_umap_final()

ggsave(
  file.path(fig_dir, "06_final_main_UMAP_named_nolabel.pdf"),
  p_nolabel,
  width = 7.5,
  height = 5.5,
  device = cairo_pdf
)

msg("Saving named object:", out_obj)
qs2::qs_save(obj, out_obj)

audit <- data.frame(
  original_label = names(name_map),
  final_label = unname(name_map),
  n_cells = as.integer(table(obj$subtype_k4)[names(name_map)]),
  neftel_anchor = c(
    "NPC2 + cycling axis",
    "OPC + mature oligodendrocyte markers",
    "MES (refined: perivascular program)",
    "MES (refined: inflammatory/myeloid program)"
  ),
  key_markers = c(
    "ASCL1, DLX5, CD24, LMO3 + KIAA0101/NUF2 (cycling)",
    "OPALIN, PLP1, MOBP, MAG, MOG, NKX6-2",
    "ACTA2, TAGLN, COL4A1, FRZB, CAVIN1",
    "C1QC, C3, FCGR3A, CYBB, IL1B, CCL3/4"
  ),
  literature_ref = c(
    "Neftel 2019 (NPC2 + cycling axis)",
    "Neftel 2019 (OPC); Couturier 2020",
    "Richards 2021 (MES2-PVN); Garofano 2021 (MTC)",
    "Hara 2021 (MES myeloid-mimicry); Greenwald 2024"
  ),
  stringsAsFactors = FALSE
)

write.csv(audit, file.path(tab_dir, "06_subtype_naming_audit.csv"), row.names = FALSE)

cat("\n== Naming audit ==\n")
print(audit)
cat("\nDone.\n")
cat("Final object  :", out_obj, "\n")
cat("Main figure   :", file.path(fig_dir, "06_final_main_UMAP_named.pdf"), "\n")
cat("No-label fig  :", file.path(fig_dir, "06_final_main_UMAP_named_nolabel.pdf"), "\n")
cat("Audit table   :", file.path(tab_dir, "06_subtype_naming_audit.csv"), "\n")
