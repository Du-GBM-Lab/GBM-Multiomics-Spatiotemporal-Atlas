#!/usr/bin/env Rscript
# =============================================================================
# R9 close-out | robustness heatmap + master summary
#
# Nature:
#   Pure close-out summary. This script only reads locked per-slice source data,
#   builds a reproducibility heatmap and writes the R9 master summary. It does
#   not rerun statistical tests, does not create new claims, and does not change
#   final placement decisions.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

base_dir <- getwd()
out_dir <- file.path(base_dir, "tables/R9_robustness_summary")
fig_dir <- file.path(base_dir, "figures/R9_robustness_summary")
doc_dir <- file.path(base_dir, "docs")
archive_root <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/10_R9_空间转录组"
archive_fig_dir <- file.path(archive_root, "补充图/R9_robustness_summary")
archive_src_dir <- file.path(archive_root, "源数据/补图_R9_robustness_summary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(doc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(archive_fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(archive_src_dir, showWarnings = FALSE, recursive = TRUE)

paths <- list(
  anchor = file.path(base_dir, "tables/C3_4_local_niche/R9_C3_4_v3B_MESV_primary_crossK_forest_source.csv"),
  a1 = file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A1_program_colocalization_toroidal/A1_program_toroidal_perslice.csv"),
  a1_q1 = file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A1_program_colocalization_toroidal/slice_check/q1_stratification_audit/A1_Q1_question3_subset_neuron_preservation.csv"),
  a1_vascular_total = file.path(base_dir, "tables/R9_batch1_unbiased_landscape/A1_program_colocalization_toroidal/slice_check/q1_stratification_audit/A1_Q1_slice_vascular_total_18_values.csv"),
  territory = file.path(base_dir, "tables/R9_batch2_landscape_gradient/territory_redo/subtype_dual_method_consistency_perslice.csv"),
  territory_neuron = file.path(base_dir, "tables/R9_batch2_landscape_gradient/territory_redo/neuron_gate_combined.csv"),
  stop58 = file.path(base_dir, "docs/R9_MESV_vascular_anchor_restore_diagnostic_STOP.md")
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing)) stop("Missing locked source(s): ", paste(missing, collapse = ", "))

anchor <- fread(paths$anchor)
a1 <- fread(paths$a1)
a1_q1 <- fread(paths$a1_q1)
a1_vascular_total <- fread(paths$a1_vascular_total)
territory <- fread(paths$territory)
territory_neuron <- fread(paths$territory_neuron)

slices <- sort(unique(anchor[row_type == "slice", slice]))
slice_meta <- data.table(slice = slices)
slice_meta[, IDH_status := fifelse(grepl("IDHMutant", slice), "IDH-Mutant", "IDH-WT")]

add_col <- function(dt, col_id, col_label, status_col, source, role = "signal") {
  x <- copy(dt)
  x[, `:=`(column_id = col_id, column_label = col_label, source = source, role = role)]
  setnames(x, status_col, "status")
  x[, .(slice, column_id, column_label, status, source, role)]
}

heat <- list()

for (kk in c(6L, 12L)) {
  v <- anchor[row_type == "slice" & response == "MES-V" & predictor == "vascular" & top_rule == "top10" & k == kk,
    .(slice, status = fifelse(emp_sig == TRUE, "pass", "fail"))
  ]
  heat[[paste0("A_k", kk)]] <- add_col(v, paste0("A_MESV_vascular_crossK_k", kk), paste0("A MES-V x vascular cross-K k=", kk), "status", paths$anchor)

  n <- anchor[row_type == "slice" & response == "MES-V" & predictor == "neuron_control" & top_rule == "top10" & k == kk,
    .(slice, status = fifelse(emp_sig == FALSE, "clean", "fail"))
  ]
  heat[[paste0("A_neuron_k", kk)]] <- add_col(n, paste0("A_neuron_control_clean_k", kk), paste0("A neuron clean k=", kk), "status", paths$anchor, "neuron_control")
}

a1_main <- a1[pair == "MESV_vascular" & role == "primary", .(slice, status = fifelse(pass_positive_fdr == TRUE, "pass", "fail"))]
heat[["A1_frozen"]] <- add_col(a1_main, "A1_program_toroidal_frozen", "A1 program toroidal frozen", "status", paths$a1)

a1_neuron <- a1[pair == "neuron_vascular" & role == "neuron_control", .(slice, status = fifelse(pass_positive_fdr == FALSE, "clean", "fail"))]
heat[["A1_neuron"]] <- add_col(a1_neuron, "A1_neuron_clean", "A1 neuron clean", "status", paths$a1, "neuron_control")

q1_slices <- unlist(strsplit(a1_q1$subset_slices[1], ";", fixed = TRUE), use.names = FALSE)
q1_slices <- trimws(q1_slices)
a1_q1_status <- merge(slice_meta[, .(slice)], a1_main, by = "slice", all.x = TRUE)
a1_q1_status[, status := fifelse(slice %in% q1_slices, status, "not_captured")]
heat[["A1_Q1"]] <- add_col(a1_q1_status, "A1_program_Q1_vascular_sufficient", "A1 vascular-sufficient subset", "status", paths$a1_q1)

for (lab in c("OPC-M", "NPC-P", "MES-I", "MES-V")) {
  t <- territory[label == lab, .(slice, status = fifelse(both_pass == TRUE, "pass", "fail"))]
  heat[[paste0("territory_", lab)]] <- add_col(t, paste0("territory_", gsub("[^A-Za-z0-9]+", "_", lab)), paste0("Territory ", lab, " both-pass"), "status", paths$territory)
}

tn <- territory_neuron[, .(slice, status = fifelse(eligible == TRUE, fifelse(join_pass == TRUE & ripley_pass == TRUE, "clean", "fail"), "not_captured"))]
heat[["territory_neuron"]] <- add_col(tn, "territory_neuron_gate", "Territory neuron gate", "status", paths$territory_neuron, "neuron_control")

heatmap_dt <- rbindlist(heat, use.names = TRUE)
heatmap_dt <- merge(slice_meta, heatmap_dt, by = "slice", all.y = TRUE)

col_order <- c(
  "A_MESV_vascular_crossK_k6",
  "A_MESV_vascular_crossK_k12",
  "A_neuron_control_clean_k6",
  "A_neuron_control_clean_k12",
  "A1_program_toroidal_frozen",
  "A1_program_Q1_vascular_sufficient",
  "A1_neuron_clean",
  "territory_OPC_M",
  "territory_NPC_P",
  "territory_MES_I",
  "territory_MES_V",
  "territory_neuron_gate"
)
heatmap_dt[, column_id := factor(column_id, levels = col_order)]
heatmap_dt[, slice := factor(slice, levels = rev(slices))]

count_check <- heatmap_dt[, .(
  pass = sum(status == "pass", na.rm = TRUE),
  clean = sum(status == "clean", na.rm = TRUE),
  fail = sum(status == "fail", na.rm = TRUE),
  not_captured = sum(status == "not_captured", na.rm = TRUE),
  denominator = .N,
  source = first(source)
), by = .(column_id, column_label, role)]
count_check[, reported_count := fifelse(role == "neuron_control", paste0(clean, "/", denominator - not_captured, " clean"), paste0(pass, "/", denominator - not_captured, " pass"))]

expected <- data.table(
  column_id = factor(col_order, levels = col_order),
  expected_count = c("15/18 pass", "15/18 pass", "17/18 clean", "16/18 clean", "14/18 pass", "12/14 pass", "18/18 clean", "17/18 pass", "16/18 pass", "16/18 pass", "14/18 pass", "10/10 clean")
)
count_check <- merge(count_check, expected, by = "column_id", all.x = TRUE)
count_check[, count_matches_locked := reported_count == expected_count]

fwrite(heatmap_dt, file.path(out_dir, "R9_robustness_heatmap_source.csv"))
fwrite(count_check, file.path(out_dir, "R9_robustness_count_check.csv"))
fwrite(data.table(source_name = names(paths), path = unlist(paths)), file.path(out_dir, "R9_robustness_sources.csv"))

plot_dt <- copy(heatmap_dt)
plot_dt[, status_plot := factor(status, levels = c("pass", "clean", "fail", "not_captured"))]
label_map <- unique(plot_dt[, .(column_id, column_label)])
label_vec <- setNames(label_map$column_label, label_map$column_id)

p <- ggplot(plot_dt, aes(x = column_id, y = slice, fill = status_plot)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = fifelse(status == "pass", "P", fifelse(status == "clean", "C", fifelse(status == "fail", "F", "NA")))), size = 2.2, color = "black") +
  scale_x_discrete(labels = label_vec, position = "top") +
  scale_fill_manual(
    values = c(pass = "#2E7D32", clean = "#4C78A8", fail = "#C62828", not_captured = "#D9D9D9"),
    breaks = c("pass", "clean", "fail", "not_captured"),
    labels = c("Pass", "Neuron clean", "Fail", "Not captured / N.A."),
    drop = FALSE
  ) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5),
    axis.text.y = element_text(size = 7),
    legend.position = "bottom"
  )

ggsave(file.path(fig_dir, "R9_robustness_heatmap.pdf"), p, width = 9.2, height = 6.2)
ggsave(file.path(fig_dir, "R9_robustness_heatmap.png"), p, width = 9.2, height = 6.2, dpi = 300)
file.copy(file.path(fig_dir, "R9_robustness_heatmap.pdf"), archive_fig_dir, overwrite = TRUE)
file.copy(file.path(fig_dir, "R9_robustness_heatmap.png"), archive_fig_dir, overwrite = TRUE)
file.copy(list.files(out_dir, full.names = TRUE), archive_src_dir, overwrite = TRUE)
script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1]
script_path <- if (length(script_arg) && !is.na(script_arg)) sub("^--file=", "", script_arg) else "scripts/59_R9_robustness_heatmap_and_master_summary.R"
file.copy(normalizePath(script_path, winslash = "/", mustWork = TRUE), archive_src_dir, overwrite = TRUE)

master_lines <- c(
  "# R9 空间转录组 — Master Summary（单一信源）",
  "",
  "## 顶部强制 caveat（写作时绝不可丢）",
  "",
  "- A1：正文必须同时写 `vascular-sufficient 14片中 12/14` 和 `全18片 frozen primary 14/18`；禁写 18/18，禁写“排除弱片后全过”。`#UKF304_T_ST` 与 `#UKF313_T_ST` 仍在 14 片子集内且未 pass。",
  "- 领地结构：正文主语为“四恶性亚型整体”斑块化，非 MES-V 专属；四个数全列：OPC-M 17/18、NPC-P 16/18、MES-I 16/18、MES-V 14/18。MES-V 14/18 仅作一致性趋势，禁单独抬为正文专属主证据；禁演化、过渡、时间读法。",
  "- border-core：最终 internal-audit / 弃用图件。若内部提及，必须标明仅 11 张切片捕获 low-to-high malignant-density axis；不得示作 18 片，不进正文/补充图。",
  "- 所有“差一点”结论：per-slice 实测数与落位同时摆出，禁把未过线数字写成过线。",
  "",
  "## 措辞天花板",
  "",
  "- 允许：spatially co-enriched / proximal / spatial preference / consistent with / spatially consistent with inferred communication / static self-clustering / patchiness / along the malignant-density axis。",
  "- 禁止：contact / recruitment / interaction（因果义）/ drives / gates / proves / causality / confirmed loop / infiltration front（动态）/ invasion / evolution / trajectory-direction / physical or single-cell colocalization / TAM recruitment / hypoxia drives MES。",
  "",
  "## 正文层（R9 空间佐证；关系承重唯一在档 A）",
  "",
  "- 档 A — MES-V × vascular niche proximity：正文主图。v3-B cross-K + forest，保自相关 null，neuron 干净，15/18（k=6/k=12），median p=0.001，count 守恒。58_ 独立复核稳健，逐片 delta max 0.00938。R9 relation-level 唯一主证据。",
  "- A1 — MES-V/vascular program 空间共定位：正文，分层敏感性，fig.SX。toroidal 检验，frozen 14/18 + vascular-sufficient 14 片中 12/14。caveat 见顶部。",
  "- 领地结构 — 四亚型整体空间斑块化：正文。join-count 解析 null + univariate Ripley 双法，neuron gate 10/10，count 72/72 守恒。caveat 见顶部。",
  "- IVY ROI 四亚型占比：正文总览。",
  "",
  "## 补充层",
  "",
  "- A2/D02 无监督域：MES-V+vascular 共富集；vascular 12/18，6 片三分层均未救活；正文一句引用。",
  "- ③ 双轴血管 × 缺氧：MES-V 偏血管轴而非缺氧轴，作为锁定结论 #2 的二维印证。",
  "- B1 一句 negative-result 回应：COMMOT 未提供特异空间支持，正面回应原稿单细胞 CellChat 套用。",
  "- 复现热图 / robustness 总览：本次新增，作为主要空间结论跨片复现的汇总图，不是新分析。",
  "",
  "## Internal-Audit（备查，不进图位）",
  "",
  "- A3 squidpy：默认 null 为禁用 label permutation，z 不采信；myeloid 自富集锁伪影。",
  "- A4 cellular neighborhood：CN05/CN07 分离，诚实弱。",
  "- B1 全库扫：超时残缺，锁 audit。",
  "- 领地 51_：坏 null 版，count 不守恒，已被 55_ 覆盖。",
  "- ② 距血管梯度：非单调，CV 2.7，退。",
  "- border-core：malignant/non-malignant 近恒等式无意义，亚型层 1/4 clean 无特异，整体舍弃；留 `边界带无亚型特异结构` 备查。",
  "- 57_ MES-V niche 邻居画像：组合闸过严，vascular 锚未复现是闸的问题，不是档 A 不稳。",
  "- 58_ 档 A 保护性诊断：重点备查，作为档 A 稳健性的 reviewer 答辩材料。还原原条件 15/18、locked source 15/18、逐片 delta <0.01；掉线链条归因组合闸。",
  "",
  "## 结案 / 不做",
  "",
  "- B1 空间通讯：舍弃。neuron 2/18 失守，MES-V 9/18、vascular 5/18，信号不特异；补充一句 negative。",
  "- B2 轨迹：不做。数据无时间方向，铁律禁证演化方向 + drives；亚型空间格局需求由领地结构覆盖。",
  "- FOSL1/PLAUR 斑块内富集：不做。踩“按 score 筛 spot”禁令，且有循环论证风险，会反噬领地结构。",
  "",
  "## 贯穿守住的核心资产（response-to-reviewer 重点）",
  "",
  "- myeloid/TAM 空间 null 全程未被复活：A3（自富集伪影）、B1（COMMOT 表达驱动伪影）、③（缺氧轴/坏死区伪影风险）、57_（4/18 FDR 未过 gate）四次独立考验均守住。这是回应审稿人 TAM 诉求的核心资产。",
  "- relation-level 判定唯一出口 = v3-B：cross-K + 保自相关 null + neuron 裁判 + count 守恒 + per-slice n=18。描述性指标（Jaccard/overlap/梯度曲线）永不配 null 充当关系检验。",
  "- 因果归 R5/R6/R7；R9 全程是 association/context 级空间佐证，reviewer #2 “共现写成 causality”诉求已全部降级。",
  "- 每个“差一点”诚实标注：A1 14/18、领地 MES-V 14/18 均未冒充过线；A1 靠合规分层敏感性正当进正文，非下调门槛。",
  "",
  "## 可复现三关确认",
  "",
  "- 方法可复现：58_ 逐片复现档 A，delta <0.01。",
  "- 禁用找不回旧对象：58_ 验证当前代码/当前 RCTD weights 可复现 locked source。",
  "- 参数、版本、来源均写入 source data / STOP。"
)
writeLines(master_lines, file.path(doc_dir, "R9_master_summary.md"), useBytes = TRUE)

stop_lines <- c(
  "# R9 close-out robustness heatmap + master summary STOP",
  "",
  "性质声明：本轮只汇总已锁 per-slice source data，不重算任何检验、不新增关系主张、不改变落位。",
  "",
  "## Robustness Heatmap",
  paste0("- output dir：`", out_dir, "`"),
  paste0("- figure：`", file.path(fig_dir, "R9_robustness_heatmap.pdf"), "`"),
  paste0("- archive figure：`", file.path(archive_fig_dir, "R9_robustness_heatmap.pdf"), "`"),
  "- count consistency check:",
  paste(capture.output(print(count_check[, .(column_label, reported_count, expected_count, count_matches_locked)])), collapse = "\n"),
  "",
  "## Master Summary",
  paste0("- master summary：`", file.path(doc_dir, "R9_master_summary.md"), "`"),
  "",
  "## Guardrail",
  "- 不下新判决；热图是已有结论的复现可视化，不是新证据。"
)
writeLines(stop_lines, file.path(doc_dir, "R9_closeout_robustness_master_STOP.md"), useBytes = TRUE)

message("Wrote robustness heatmap and master summary.")
