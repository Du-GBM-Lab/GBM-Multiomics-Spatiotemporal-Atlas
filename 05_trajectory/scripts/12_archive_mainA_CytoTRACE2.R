options(stringsAsFactors = FALSE)

root <- normalizePath("<DATA_ROOT>/项目/分型/修稿杠生信/重新分析", winslash = "/", mustWork = TRUE)
archive <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/04_拟时序分析"

dir.create(file.path(archive, "正文图"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(archive, "补充图"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(archive, "源数据", "正图A_NPCP根_CytoTRACE2"), recursive = TRUE, showWarnings = FALSE)

fig_src_dir <- file.path(root, "06_恶性细胞拟时序/figures/final_panels/main")
src_data_dir <- file.path(root, "06_恶性细胞拟时序/figures/final_panels/source_data/main")

fig_map <- data.frame(
  category = "main",
  type = "figure",
  source = file.path(fig_src_dir, c("Main_A_root_CytoTRACE2.pdf", "Main_A_root_CytoTRACE2.png")),
  destination = file.path(archive, "正文图", c("正图A_NPCP根_CytoTRACE2.pdf", "正图A_NPCP根_CytoTRACE2.png")),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(fig_map))) {
  file.copy(fig_map$source[i], fig_map$destination[i], overwrite = TRUE)
}

source_files <- c("Main_A_CytoTRACE2_by_subtype.csv", "Main_A_root_metric_summary.csv")
sd_map <- data.frame(
  category = "main",
  type = "source_data",
  source = file.path(src_data_dir, source_files),
  destination = file.path(archive, "源数据", "正图A_NPCP根_CytoTRACE2", source_files),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(sd_map))) {
  file.copy(sd_map$source[i], sd_map$destination[i], overwrite = TRUE)
}

cyto <- read.csv(file.path(src_data_dir, "Main_A_CytoTRACE2_by_subtype.csv"), check.names = FALSE)
cyto_summary <- do.call(rbind, lapply(split(cyto$CytoTRACE2_score, cyto$subtype_short), function(x) {
  data.frame(
    n_cells = length(x),
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    q25 = as.numeric(quantile(x, 0.25, na.rm = TRUE)),
    q75 = as.numeric(quantile(x, 0.75, na.rm = TRUE)),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
}))
cyto_summary <- cbind(subtype_short = rownames(cyto_summary), cyto_summary)
rownames(cyto_summary) <- NULL
cyto_summary <- cyto_summary[match(c("NPC-P", "OPC-M", "MES-V", "MES-I"), cyto_summary$subtype_short), ]

summary_dest <- file.path(archive, "源数据", "正图A_NPCP根_CytoTRACE2", "正图A_CytoTRACE2统计摘要.csv")
write.csv(cyto_summary, summary_dest, row.names = FALSE, fileEncoding = "UTF-8")

stats_src <- file.path(root, "06_恶性细胞拟时序/tables/拟时序_root候选统计检验.csv")
stats <- read.csv(stats_src, check.names = FALSE)
stats <- stats[stats$table %in% c("CytoTRACE2_KW", "CytoTRACE2_Dunn_vs_NPC_P"), ]
stats_dest <- file.path(archive, "源数据", "正图A_NPCP根_CytoTRACE2", "正图A_CytoTRACE2统计检验.csv")
write.csv(stats, stats_dest, row.names = FALSE, fileEncoding = "UTF-8")

now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
fig_manifest <- transform(fig_map, timestamp = now)
sd_manifest <- rbind(
  transform(sd_map, timestamp = now),
  data.frame(
    category = "main",
    type = "source_data",
    source = c(summary_dest, stats_dest),
    destination = c(summary_dest, stats_dest),
    timestamp = now,
    stringsAsFactors = FALSE
  )
)

write.csv(fig_manifest, file.path(archive, "源数据", "图文件复制清单.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(sd_manifest, file.path(archive, "源数据", "源数据复制清单.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("Archived Main A to:", archive, "\n")
print(cyto_summary)
print(stats[, intersect(c("table", "metric", "test", "p_value", "comparison", "z", "p_adj", "star"), colnames(stats))])
