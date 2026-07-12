options(stringsAsFactors = FALSE)

root <- normalizePath("<DATA_ROOT>/项目/分型/修稿杠生信/重新分析", winslash = "/", mustWork = TRUE)
archive <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/04_拟时序分析"

dir.create(file.path(archive, "正文图"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(archive, "源数据", "正图D_MES终点程序"), recursive = TRUE, showWarnings = FALSE)

fig_src <- file.path(root, "06_恶性细胞拟时序/figures/final_panels/main/Main_E_MES_endpoint_comparison.pdf")
fig_dest <- file.path(archive, "正文图", "正图D_MES终点程序.pdf")
file.copy(fig_src, fig_dest, overwrite = TRUE)

src_data_dir <- file.path(root, "06_恶性细胞拟时序/figures/final_panels/source_data/main")
fig_source_dir <- file.path(root, "06_恶性细胞拟时序/figures/source_data")
table_dir <- file.path(root, "06_恶性细胞拟时序/tables")

source_files <- c(
  file.path(src_data_dir, "Main_E_MES_endpoint_boxplot.csv"),
  file.path(fig_source_dir, "06_panel_D_endpoint_boxplot.csv"),
  file.path(table_dir, "terminal_endpoint_comparison.csv"),
  file.path(table_dir, "terminal_endpoint_composition.csv"),
  file.path(table_dir, "terminal_marker_validation.csv"),
  file.path(table_dir, "terminal_marker_validation_marker_summary.csv"),
  file.path(table_dir, "terminal_marker_validation_reference_marker_summary.csv"),
  file.path(table_dir, "paga_connectivity_matrix.csv")
)
source_files <- source_files[file.exists(source_files)]

dest_files <- file.path(archive, "源数据", "正图D_MES终点程序", basename(source_files))
for (i in seq_along(source_files)) {
  file.copy(source_files[i], dest_files[i], overwrite = TRUE)
}

endpoint <- read.csv(file.path(table_dir, "terminal_endpoint_comparison.csv"), check.names = FALSE)
endpoint_wide <- endpoint |>
  subset(feature %in% c("MP02_UCell", "MP04_UCell", "RGS5", "ACTA2", "TAGLN", "HLA-DRA", "HLA-DPB1", "CD74")) |>
  subset(endpoint_group == "MES-I terminal") |>
  transform(
    display_feature = feature_label,
    MESI_mean = mean_MESI,
    MESV_mean = mean_MESV,
    delta_MESV_minus_MESI = delta_MESV_minus_MESI,
    BH_p = p_adj_BH
  ) |>
  subset(select = c("feature", "display_feature", "feature_class", "MESI_mean", "MESV_mean", "delta_MESV_minus_MESI", "BH_p"))
endpoint_wide <- endpoint_wide[order(endpoint_wide$feature_class, endpoint_wide$feature), ]
endpoint_summary_dest <- file.path(archive, "源数据", "正图D_MES终点程序", "正图D_endpoint程序差异摘要.csv")
write.csv(endpoint_wide, endpoint_summary_dest, row.names = FALSE, fileEncoding = "UTF-8")

conn <- read.csv(file.path(table_dir, "paga_connectivity_matrix.csv"), check.names = FALSE)
conn_pairs <- subset(conn, from_cluster != to_cluster)
if ("connectivity" %in% names(conn_pairs)) {
  conn_pairs <- conn_pairs[order(-conn_pairs$connectivity), ]
} else if ("observed_expected_ratio" %in% names(conn_pairs)) {
  conn_pairs <- conn_pairs[order(-conn_pairs$observed_expected_ratio), ]
}
conn_summary_dest <- file.path(archive, "源数据", "正图D_MES终点程序", "正图D_MES相邻性依据.csv")
write.csv(conn_pairs, conn_summary_dest, row.names = FALSE, fileEncoding = "UTF-8")

now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
fig_manifest_path <- file.path(archive, "源数据", "图文件复制清单.csv")
src_manifest_path <- file.path(archive, "源数据", "源数据复制清单.csv")

read_manifest <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(
      category = character(),
      type = character(),
      source = character(),
      destination = character(),
      timestamp = character(),
      stringsAsFactors = FALSE
    ))
  }
  tryCatch(
    read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM"),
    error = function(e) data.frame(
      category = character(),
      type = character(),
      source = character(),
      destination = character(),
      timestamp = character(),
      stringsAsFactors = FALSE
    )
  )
}

fig_manifest_old <- read_manifest(fig_manifest_path)
src_manifest_old <- read_manifest(src_manifest_path)

fig_manifest_new <- data.frame(
  category = "main",
  type = "figure",
  source = fig_src,
  destination = fig_dest,
  timestamp = now,
  stringsAsFactors = FALSE
)
src_manifest_new <- data.frame(
  category = "main",
  type = "source_data",
  source = c(source_files, endpoint_summary_dest, conn_summary_dest),
  destination = c(dest_files, endpoint_summary_dest, conn_summary_dest),
  timestamp = now,
  stringsAsFactors = FALSE
)

fig_manifest <- rbind(fig_manifest_old, fig_manifest_new)
src_manifest <- rbind(src_manifest_old, src_manifest_new)
fig_manifest <- fig_manifest[!duplicated(fig_manifest$destination, fromLast = TRUE), ]
src_manifest <- src_manifest[!duplicated(src_manifest$destination, fromLast = TRUE), ]

write.csv(fig_manifest, fig_manifest_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(src_manifest, src_manifest_path, row.names = FALSE, fileEncoding = "UTF-8")

cat("Archived MES endpoint main figure to:", fig_dest, "\n")
cat("Endpoint program summary:\n")
print(endpoint_wide, row.names = FALSE)
cat("Top connectivity pairs:\n")
print(head(conn_pairs, 8), row.names = FALSE)
