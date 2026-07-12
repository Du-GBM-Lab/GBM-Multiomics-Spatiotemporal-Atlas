options(stringsAsFactors = FALSE)

root <- normalizePath("<DATA_ROOT>/项目/分型/修稿杠生信/重新分析", winslash = "/", mustWork = TRUE)
archive <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/04_拟时序分析"

dir.create(file.path(archive, "补充图"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(archive, "源数据", "补图S_root稳健性"), recursive = TRUE, showWarnings = FALSE)

fig_src <- file.path(root, "06_恶性细胞拟时序/figures/final_panels/main/Main_D_topology_robustness.pdf")
fig_dest <- file.path(archive, "补充图", "补图S_root稳健性.pdf")
file.copy(fig_src, fig_dest, overwrite = TRUE)

src_data_dir <- file.path(root, "06_恶性细胞拟时序/figures/final_panels/source_data/main")
table_dir <- file.path(root, "06_恶性细胞拟时序/tables")

source_files <- c(
  file.path(src_data_dir, "Main_D_multiroot_curves.csv"),
  file.path(src_data_dir, "Main_D_multiroot_lineage_meta.csv"),
  file.path(table_dir, "multiroot_lineage_summary.csv"),
  file.path(table_dir, "multiroot_pseudotime_correlation.csv"),
  file.path(table_dir, "L3_robustness_summary.csv")
)
dest_files <- file.path(archive, "源数据", "补图S_root稳健性", basename(source_files))
for (i in seq_along(source_files)) file.copy(source_files[i], dest_files[i], overwrite = TRUE)

lineage_summary <- read.csv(file.path(table_dir, "multiroot_lineage_summary.csv"), check.names = FALSE)
lineage_meta <- read.csv(file.path(src_data_dir, "Main_D_multiroot_lineage_meta.csv"), check.names = FALSE)
correlation <- read.csv(file.path(table_dir, "multiroot_pseudotime_correlation.csv"), check.names = FALSE)
l3 <- read.csv(file.path(table_dir, "L3_robustness_summary.csv"), check.names = FALSE)

lineage_summary_dest <- file.path(archive, "源数据", "补图S_root稳健性", "补图S_multiroot摘要.csv")
write.csv(lineage_summary, lineage_summary_dest, row.names = FALSE, fileEncoding = "UTF-8")

lineage_meta_dest <- file.path(archive, "源数据", "补图S_root稳健性", "补图S_multiroot路径明细.csv")
write.csv(lineage_meta, lineage_meta_dest, row.names = FALSE, fileEncoding = "UTF-8")

cor_dest <- file.path(archive, "源数据", "补图S_root稳健性", "补图S_root_relative相关性.csv")
write.csv(correlation, cor_dest, row.names = FALSE, fileEncoding = "UTF-8")

l3_dest <- file.path(archive, "源数据", "补图S_root稳健性", "补图S_L3稳健性.csv")
write.csv(l3, l3_dest, row.names = FALSE, fileEncoding = "UTF-8")

now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
fig_manifest_path <- file.path(archive, "源数据", "图文件复制清单.csv")
src_manifest_path <- file.path(archive, "源数据", "源数据复制清单.csv")

read_manifest <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(category = character(), type = character(), source = character(), destination = character(), timestamp = character(), stringsAsFactors = FALSE))
  }
  tryCatch(
    read.csv(path, check.names = FALSE, fileEncoding = "UTF-8-BOM"),
    error = function(e) data.frame(category = character(), type = character(), source = character(), destination = character(), timestamp = character(), stringsAsFactors = FALSE)
  )
}

fig_manifest_old <- read_manifest(fig_manifest_path)
src_manifest_old <- read_manifest(src_manifest_path)

fig_manifest_new <- data.frame(category = "supplement", type = "figure", source = fig_src, destination = fig_dest, timestamp = now, stringsAsFactors = FALSE)
src_manifest_new <- data.frame(
  category = "supplement",
  type = "source_data",
  source = c(source_files, lineage_summary_dest, lineage_meta_dest, cor_dest, l3_dest),
  destination = c(dest_files, lineage_summary_dest, lineage_meta_dest, cor_dest, l3_dest),
  timestamp = now,
  stringsAsFactors = FALSE
)

fig_manifest <- rbind(fig_manifest_old, fig_manifest_new)
src_manifest <- rbind(src_manifest_old, src_manifest_new)
fig_manifest <- fig_manifest[!duplicated(fig_manifest$destination, fromLast = TRUE), ]
src_manifest <- src_manifest[!duplicated(src_manifest$destination, fromLast = TRUE), ]

write.csv(fig_manifest, fig_manifest_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(src_manifest, src_manifest_path, row.names = FALSE, fileEncoding = "UTF-8")

cat("Archived root-sensitivity supplement to:", fig_dest, "\n")
cat("Multiroot lineage summary:\n")
print(lineage_summary)
cat("Lineage path details:\n")
print(lineage_meta)
cat("Root-relative pseudotime correlations:\n")
print(correlation)
cat("L3 robustness:\n")
print(l3)
