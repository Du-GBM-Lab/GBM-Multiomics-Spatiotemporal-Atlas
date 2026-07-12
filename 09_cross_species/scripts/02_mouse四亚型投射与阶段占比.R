#!/usr/bin/env Rscript
# =============================================================================
# 07_发育时间_TF_ATAC验证 / 01_mouse_developmental
# Phase 1 - Script 02: 人四亚型 signature -> 小鼠恶性细胞 跨物种投射 + 置信度评估
# -----------------------------------------------------------------------------
# 这是 cross-species representation/scoring, 不是 trajectory。
# 下游表述统一用:
# "cross-species subtype representation across tumorigenesis stages"。
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(babelgene)
  library(UCell)
})

set.seed(20260524)

CFG <- list(
  top_n = 80,
  min_pct_delta = 0.10,
  min_genes_per_sig = 10,
  low_conf_quantile = 0.10
)

module_dir <- normalizePath("07_发育时间_TF_ATAC验证", winslash = "/", mustWork = TRUE)
project_root <- dirname(module_dir)
project_dir <- file.path(module_dir, "01_mouse_developmental")
data_dir <- file.path(project_dir, "data")
tables_dir <- file.path(module_dir, "tables")
outputs_dir <- file.path(project_dir, "outputs")
figures_dir <- file.path(project_dir, "figures")
src_data_dir <- file.path(figures_dir, "source_data")
logs_dir <- file.path(project_dir, "logs")

for (d in c(tables_dir, outputs_dir, figures_dir, src_data_dir, logs_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

mouse_obj_path <- file.path(data_dir, "GSE278511小鼠恶性模型_精简.qs2")
human_marker_path <- file.path(
  project_root,
  "05_恶性细胞分亚群与Neftel对照/tables/10e_one_vs_rest_sparse_avg_log2FC_rank_combined.csv"
)
internal_label_path <- file.path(
  project_root,
  "05_恶性细胞分亚群与Neftel对照/tables/11_subtype_internal_naming_locked.csv"
)

stopifnot(file.exists(mouse_obj_path))
stopifnot(file.exists(human_marker_path))
stopifnot(file.exists(internal_label_path))

# 注意：当前项目框架里 MES-I = MES-Immune-interacting-like。
# 既往锁定 CSV 里可能写作 MES-Antigen-presenting；这是更强、更具体的 claim。
# 这里统一映射到短名 MES-I，并在 QC 日志中显式提醒后续全文统一。
short_label_map <- c(
  "Proliferative-NPC" = "NPC-P",
  "Proliferative-NPC subtype" = "NPC-P",
  "OPC-Myelination" = "OPC-M",
  "OPC-Myelination subtype" = "OPC-M",
  "Vascular-niche MES" = "MES-V",
  "Vascular-niche MES subtype" = "MES-V",
  "MES-Antigen-presenting" = "MES-I",
  "MES-Antigen-presenting subtype" = "MES-I",
  "MES-Immune-interacting" = "MES-I",
  "MES-Immune-interacting subtype" = "MES-I"
)

subtype_levels <- c("NPC-P", "OPC-M", "MES-V", "MES-I")
stage_levels <- c("Preneoplastic", "Early lesion", "Mid lesion", "Endpoint")

qc <- character(0)
log_line <- function(...) {
  line <- paste0(...)
  qc <<- c(qc, line)
  message(line)
}

seurat_v5 <- utils::packageVersion("SeuratObject") >= "5.0.0"
get_counts <- function(obj, assay = "RNA") {
  if (seurat_v5) {
    SeuratObject::LayerData(obj, assay = assay, layer = "counts")
  } else {
    Seurat::GetAssayData(obj, assay = assay, slot = "counts")
  }
}

log_line("== Phase 1 / Script 02 :: mouse cross-species subtype projection ==")
log_line("[time] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

mouse_obj <- qs2::qs_read(mouse_obj_path)
log_line("[mouse] cells = ", ncol(mouse_obj), " ; features = ", nrow(mouse_obj))

if (!"stage_inferred" %in% colnames(mouse_obj@meta.data)) {
  stop("Mouse object missing 'stage_inferred'.")
}

stage_tab <- table(mouse_obj$stage_inferred, useNA = "ifany")
log_line("[stage] ", paste(names(stage_tab), as.integer(stage_tab), sep = "=", collapse = "; "))

malig_fields <- intersect(
  c("malignant", "is_malignant", "tumor_normal", "compartment", "cell_class"),
  colnames(mouse_obj@meta.data)
)
if (length(malig_fields) > 0) {
  for (f in malig_fields) {
    tab <- table(mouse_obj@meta.data[[f]], useNA = "ifany")
    log_line("[malignancy check] field '", f, "': ", paste(names(tab), tab, sep = "=", collapse = "; "))
  }
  log_line("[malignancy check] 请确认上述字段表明对象已限定为恶性细胞；若混入非恶性细胞，投射结果对那部分细胞无意义。")
} else {
  log_line("[malignancy check] 未发现显式恶性标记字段；按对象命名假定已预先过滤为恶性细胞。")
}

internal_labels <- read_csv(internal_label_path, show_col_types = FALSE)
need_cols <- c("subtype_k4", "subtype_label_final")
if (!all(need_cols %in% colnames(internal_labels))) {
  stop("internal naming CSV missing columns: ", paste(setdiff(need_cols, colnames(internal_labels)), collapse = ", "))
}

internal_labels <- internal_labels %>%
  mutate(short_label = unname(short_label_map[subtype_label_final])) %>%
  select(subtype_k4, subtype_label_final, short_label)

if (any(is.na(internal_labels$short_label))) {
  stop(
    "internal naming CSV contains subtype_label_final not covered by short_label_map: ",
    paste(unique(internal_labels$subtype_label_final[is.na(internal_labels$short_label)]), collapse = ", ")
  )
}
if (!setequal(internal_labels$short_label, subtype_levels)) {
  stop("internal naming CSV does not resolve to exactly NPC-P / OPC-M / MES-V / MES-I.")
}

write_csv(internal_labels, file.path(tables_dir, "mouse_phase1_subtype_label_audit.csv"))
log_line("[audit] four subtypes confirmed from locked CSV: ", paste(sort(unique(internal_labels$short_label)), collapse = ", "))
if (any(grepl("Antigen-presenting", internal_labels$subtype_label_final))) {
  log_line(
    "[audit][WARN] MES-I 在锁定 CSV 中名为 'MES-Antigen-presenting'；框架文档写的是 'MES-Immune-interacting'。",
    "这是 strong vs broad claim 的差别，后续正文需锁定一个名字。"
  )
}

marker_req <- c("gene", "subtype_k4", "avg_log2FC", "pct.1", "pct.2")
marker_tbl <- read_csv(human_marker_path, show_col_types = FALSE)
if (!all(marker_req %in% colnames(marker_tbl))) {
  stop("human marker CSV missing columns: ", paste(setdiff(marker_req, colnames(marker_tbl)), collapse = ", "))
}

marker_tbl <- marker_tbl %>%
  left_join(internal_labels, by = "subtype_k4") %>%
  mutate(pct_delta = pct.1 - pct.2) %>%
  filter(!is.na(short_label), avg_log2FC > 0, pct_delta >= CFG$min_pct_delta) %>%
  group_by(short_label) %>%
  arrange(desc(avg_log2FC), desc(pct_delta), .by_group = TRUE) %>%
  slice_head(n = CFG$top_n) %>%
  ungroup()

orth <- babelgene::orthologs(
  genes = unique(marker_tbl$gene),
  species = "mouse",
  human = TRUE
) %>%
  as_tibble() %>%
  select(human_symbol = human_symbol, mouse_symbol = symbol) %>%
  distinct()

mapped <- marker_tbl %>%
  left_join(orth, by = c("gene" = "human_symbol")) %>%
  mutate(mouse_present = !is.na(mouse_symbol) & mouse_symbol %in% rownames(mouse_obj))

write_csv(mapped, file.path(src_data_dir, "mouse_phase1_human_markers_mouse_orthologs.csv"))

signature_tbl <- mapped %>%
  filter(mouse_present) %>%
  distinct(short_label, mouse_symbol)

coverage <- mapped %>%
  group_by(short_label) %>%
  summarise(
    n_human_markers = n_distinct(gene),
    n_mouse_ortholog_used = n_distinct(mouse_symbol[mouse_present]),
    .groups = "drop"
  ) %>%
  mutate(coverage_fraction = n_mouse_ortholog_used / n_human_markers)

write_csv(coverage, file.path(tables_dir, "mouse_phase1_signature_ortholog_coverage.csv"))
for (i in seq_len(nrow(coverage))) {
  log_line(
    "[signature] ", coverage$short_label[i], ": ",
    coverage$n_mouse_ortholog_used[i], "/", coverage$n_human_markers[i],
    " genes mapped (coverage = ", sprintf("%.2f", coverage$coverage_fraction[i]), ")"
  )
}
if (min(coverage$coverage_fraction) < 0.5) {
  log_line("[signature][WARN] 某亚型 ortholog 覆盖率 < 0.5；z 标准化可部分抵消基线偏差，但解读需注意不对称。")
}

signatures <- split(signature_tbl$mouse_symbol, signature_tbl$short_label)
signatures <- signatures[subtype_levels]
too_few <- names(signatures)[vapply(signatures, length, integer(1)) < CFG$min_genes_per_sig]
if (length(too_few) > 0) {
  stop("Too few mouse genes for signatures: ", paste(too_few, collapse = ", "))
}

sanity_markers <- signature_tbl %>%
  left_join(
    mapped %>% select(short_label, mouse_symbol, avg_log2FC) %>% distinct(),
    by = c("short_label", "mouse_symbol")
  ) %>%
  group_by(short_label) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  select(short_label, mouse_symbol)
write_csv(sanity_markers, file.path(tables_dir, "mouse_phase1_sanity_markers.csv"))

DefaultAssay(mouse_obj) <- "RNA"
counts <- get_counts(mouse_obj, assay = "RNA")
ucell_mat <- UCell::ScoreSignatures_UCell(
  matrix = counts,
  features = signatures,
  ncores = 1
)
colnames(ucell_mat) <- sub("_UCell$", "", colnames(ucell_mat))
ucell_mat <- ucell_mat[, subtype_levels, drop = FALSE]

z_mat <- scale(ucell_mat)
top_idx <- max.col(z_mat, ties.method = "first")
z_sorted <- t(apply(z_mat, 1, sort, decreasing = TRUE))
margin_z <- z_sorted[, 1] - z_sorted[, 2]

margin_cut <- as.numeric(quantile(margin_z, CFG$low_conf_quantile))
conf <- ifelse(margin_z >= margin_cut, "high", "low")

mouse_obj$assigned_subtype <- factor(subtype_levels[top_idx], levels = subtype_levels)
mouse_obj$assigned_subtype_score <- z_sorted[, 1]
mouse_obj$assigned_subtype_margin <- margin_z
mouse_obj$assignment_confidence <- factor(conf, levels = c("high", "low"))
mouse_obj$assigned_subtype_hc <- factor(
  ifelse(conf == "high", as.character(mouse_obj$assigned_subtype), NA),
  levels = subtype_levels
)

for (s in subtype_levels) {
  mouse_obj@meta.data[[paste0("UCell_", s)]] <- ucell_mat[, s]
}

log_line("[assign] UCell argmax after per-signature z-normalisation. margin cutoff q", CFG$low_conf_quantile, " = ", sprintf("%.3f", margin_cut))
log_line("[assign] high-confidence cells: ", sum(conf == "high"), " / ", length(conf))
for (s in subtype_levels) {
  log_line("[assign]   ", s, ": all = ", sum(mouse_obj$assigned_subtype == s), " ; high-conf = ", sum(mouse_obj$assigned_subtype_hc == s, na.rm = TRUE))
}

agreement <- NA_real_
sec_ok <- FALSE
try({
  mouse_obj <- Seurat::AddModuleScore(
    mouse_obj,
    features = signatures,
    name = "MS_",
    assay = "RNA"
  )
  ms_cols <- paste0("MS_", seq_along(signatures))
  ms_mat <- as.matrix(mouse_obj@meta.data[, ms_cols, drop = FALSE])
  colnames(ms_mat) <- subtype_levels
  ms_assign <- subtype_levels[max.col(scale(ms_mat), ties.method = "first")]
  agreement <- mean(ms_assign == as.character(mouse_obj$assigned_subtype))
  mouse_obj$assigned_subtype_addmodule <- factor(ms_assign, levels = subtype_levels)
  sec_ok <- TRUE
}, silent = TRUE)

if (sec_ok) {
  log_line("[robustness] UCell vs AddModuleScore agreement = ", sprintf("%.3f", agreement))
  if (agreement < 0.8) {
    log_line("[robustness][WARN] 两种打分法一致性 < 0.8；亚型归类对方法敏感，正文需弱化并在补充中并列两种结果。")
  }
} else {
  log_line("[robustness][WARN] AddModuleScore 第二方法未能运行；仅有 UCell 单一方法结果。")
}

md <- mouse_obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  mutate(
    stage_inferred = factor(stage_inferred, levels = stage_levels),
    assigned_subtype = factor(assigned_subtype, levels = subtype_levels)
  )

per_cell <- md %>%
  select(
    cell,
    any_of(c("sample", "sample_group", "seurat_clusters")),
    stage_inferred,
    assigned_subtype,
    assigned_subtype_hc,
    assignment_confidence,
    assigned_subtype_score,
    assigned_subtype_margin,
    all_of(paste0("UCell_", subtype_levels)),
    any_of("assigned_subtype_addmodule")
  )
write_csv(per_cell, file.path(src_data_dir, "mouse_phase1_per_cell_assigned_subtype.csv"))

counts_full <- md %>%
  count(stage_inferred, assigned_subtype, name = "n_cells", .drop = FALSE)
write_csv(counts_full, file.path(src_data_dir, "mouse_phase1_stage_subtype_counts.csv"))

prop_all <- md %>%
  filter(!is.na(stage_inferred), !is.na(assigned_subtype)) %>%
  count(stage_inferred, assigned_subtype, name = "n_cells", .drop = FALSE) %>%
  group_by(stage_inferred) %>%
  mutate(
    stage_total = sum(n_cells),
    proportion = n_cells / stage_total
  ) %>%
  ungroup()
write_csv(prop_all, file.path(src_data_dir, "mouse_phase1_stage_subtype_proportions.csv"))

prop_hc <- md %>%
  filter(!is.na(stage_inferred), !is.na(assigned_subtype_hc)) %>%
  count(stage_inferred, assigned_subtype_hc, name = "n_cells", .drop = FALSE) %>%
  rename(assigned_subtype = assigned_subtype_hc) %>%
  group_by(stage_inferred) %>%
  mutate(
    stage_total = sum(n_cells),
    proportion = n_cells / stage_total
  ) %>%
  ungroup()
write_csv(prop_hc, file.path(src_data_dir, "mouse_phase1_stage_subtype_proportions_highconf.csv"))

log_line("[tables] per-cell / counts / proportions(all+highconf) 已写出。")
log_line("[note] Preneoplastic 阶段仅约 455 个恶性细胞；该阶段比例噪声大，出图务必标注 n。")

qs2::qs_save(
  mouse_obj,
  file.path(outputs_dir, "GSE278511小鼠恶性模型_四亚型投射.qs2")
)

writeLines(qc, file.path(logs_dir, "02_小鼠四亚型投射_QC日志.txt"))
writeLines(
  capture.output(sessionInfo()),
  file.path(logs_dir, "02_小鼠四亚型投射_sessionInfo.txt")
)

log_line("== done ==")
