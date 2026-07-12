#!/usr/bin/env Rscript
# =============================================================================
# R9 | Stage 3 preflight: IVY-GAP-like and hypoxia signature set audit
# Purpose:
#   Lock signature sources before scoring.
#   No spatial scoring, no redundancy test, no proximity analysis.
#
# Locked sets:
#   IVY-like: CT / CTmvp / CTpan / LE from local IVY_gap_signatures.xlsx/csv.
#   Hypoxia main: MSigDB BUFFA_HYPOXIA_METAGENE via msigdbr.
#   Hypoxia sensitivity: MSigDB HALLMARK_HYPOXIA via msigdbr.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(msigdbr)
  library(Seurat)
})

base_dir <- getwd()
st_path <- "<DATA_ROOT>/项目/分型/分型代码/0.对象/5.ST_merge.rds"
ivy_xlsx <- "<DATA_ROOT>/项目/分型/分型代码/8.122多组学空间分析（map注释）/IVY_gap_signatures.xlsx"
ivy_csv <- file.path(base_dir, "tables/C3_niche_preflight/IVY_gap_signatures.csv")
out_dir <- file.path(base_dir, "tables/R9_stage3_signature_sets")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- spatial gene universe -------------------------------------------------
st <- readRDS(st_path)
spatial_genes <- rownames(st[["Spatial"]])
sct_genes <- rownames(st[["SCT"]])

## ---- IVY source audit ------------------------------------------------------
ivy_sheets <- readxl::excel_sheets(ivy_xlsx)
stopifnot("ciberGeneSignature_rerun" %in% ivy_sheets)
ivy_x <- as.data.table(readxl::read_excel(ivy_xlsx, sheet = "ciberGeneSignature_rerun"))
ivy_c <- fread(ivy_csv)
stopifnot(identical(names(ivy_x), names(ivy_c)))
stopifnot(all(c("NAME", "Assigned") %in% names(ivy_x)))

ivy_x[, NAME := as.character(NAME)]
ivy_x[, Assigned := as.character(Assigned)]
ivy_sets <- ivy_x[, .(gene = unique(NAME)), by = .(signature = Assigned)]
ivy_sets[, source := "local_IVY_gap_signatures.xlsx:ciberGeneSignature_rerun"]
ivy_sets[, role := fifelse(signature == "CTpan", "PAN-like contextual layer",
                    fifelse(signature == "CTmvp", "MVP-like contextual layer",
                    fifelse(signature == "CT", "CT-like contextual layer",
                            "LE-like contextual layer")))]

## ---- Hypoxia source audit --------------------------------------------------
msig <- as.data.table(msigdbr(species = "Homo sapiens"))
hypoxia_names <- c("BUFFA_HYPOXIA_METAGENE", "HALLMARK_HYPOXIA")
missing <- setdiff(hypoxia_names, unique(msig$gs_name))
if (length(missing)) stop("MSigDB gene set not found: ", paste(missing, collapse = ", "))

hyp <- msig[gs_name %in% hypoxia_names, .(
  gene = unique(gene_symbol),
  gs_collection = first(gs_collection),
  gs_subcollection = first(gs_subcollection),
  gs_description = first(gs_description),
  gs_pmid = first(gs_pmid),
  db_version = first(db_version),
  gs_url = first(gs_url)
), by = .(signature = gs_name)]
hyp[, source := "msigdbr::msigdbr(species='Homo sapiens') / MSigDB"]
hyp[, role := fifelse(signature == "BUFFA_HYPOXIA_METAGENE",
                      "primary hypoxia signature",
                      "hypoxia sensitivity signature")]

all_sets <- rbindlist(list(
  ivy_sets[, .(signature, gene, source, role,
               gs_collection = NA_character_, gs_subcollection = NA_character_,
               gs_description = NA_character_, gs_pmid = NA_character_,
               db_version = NA_character_, gs_url = NA_character_)],
  hyp[, .(signature, gene, source, role,
          gs_collection, gs_subcollection, gs_description, gs_pmid,
          db_version, gs_url)]
), use.names = TRUE, fill = TRUE)

overlap <- all_sets[, .(
  n_genes = uniqueN(gene),
  n_in_spatial = uniqueN(gene[gene %in% spatial_genes]),
  n_in_sct = uniqueN(gene[gene %in% sct_genes]),
  frac_in_spatial = uniqueN(gene[gene %in% spatial_genes]) / uniqueN(gene),
  frac_in_sct = uniqueN(gene[gene %in% sct_genes]) / uniqueN(gene),
  missing_spatial = paste(setdiff(unique(gene), spatial_genes), collapse = ";")
), by = .(signature, source, role, gs_collection, gs_subcollection, gs_pmid, db_version, gs_url)]
setorder(overlap, signature)

provenance <- unique(all_sets[, .(
  signature, source, role, gs_collection, gs_subcollection,
  gs_description, gs_pmid, db_version, gs_url
)])
setorder(provenance, signature)

fwrite(all_sets, file.path(out_dir, "R9_stage3_locked_signature_genes_IVY_hypoxia.csv"))
fwrite(overlap, file.path(out_dir, "R9_stage3_locked_signature_overlap.csv"))
fwrite(provenance, file.path(out_dir, "R9_stage3_locked_signature_provenance.csv"))

cat("== Locked IVY-like signatures:\n")
print(overlap[signature %in% c("CT", "CTmvp", "CTpan", "LE"),
              .(signature, n_genes, n_in_spatial, frac_in_spatial)])

cat("\n== Locked hypoxia signatures:\n")
print(overlap[signature %in% hypoxia_names,
              .(signature, n_genes, n_in_spatial, frac_in_spatial,
                gs_collection, gs_subcollection, gs_pmid, db_version)])

cat("\n== Provenance:\n")
print(provenance[, .(signature, source, role, gs_collection, gs_subcollection,
                     gs_pmid, db_version, gs_url)])

cat("\n[STOP signature audit] Signature sets are locked, but no spatial scoring has been run.\n")
