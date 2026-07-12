options(stringsAsFactors = FALSE)

root <- normalizePath("<DATA_ROOT>/项目/分型/修稿杠生信/重新分析", winslash = "/", mustWork = TRUE)
archive <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/04_拟时序分析"

dir.create(file.path(archive, "补充图"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(archive, "源数据", "补图S1A_细胞周期"), recursive = TRUE, showWarnings = FALSE)

fig_src <- file.path(root, "06_恶性细胞拟时序/figures/final_panels/supplement/S1A_cell_cycle_phase.pdf")
fig_dest <- file.path(archive, "补充图", "补图S1A_细胞周期.pdf")
file.copy(fig_src, fig_dest, overwrite = TRUE)

src_data_dir <- file.path(root, "06_恶性细胞拟时序/figures/final_panels/source_data/supplement")
source_files <- file.path(src_data_dir, "01_cell_cycle_phase_composition.csv")
dest_files <- file.path(archive, "源数据", "补图S1A_细胞周期", basename(source_files))
file.copy(source_files, dest_files, overwrite = TRUE)

phase <- read.csv(source_files, check.names = FALSE)
phase_summary <- phase
phase_summary$cell_fraction <- ave(phase_summary$n_cells, phase_summary$subtype_short, FUN = function(x) x / sum(x))
phase_summary_dest <- file.path(archive, "源数据", "补图S1A_细胞周期", "补图S1A_细胞周期比例摘要.csv")
write.csv(phase_summary, phase_summary_dest, row.names = FALSE, fileEncoding = "UTF-8")

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
  source = c(source_files, phase_summary_dest),
  destination = c(dest_files, phase_summary_dest),
  timestamp = now,
  stringsAsFactors = FALSE
)

fig_manifest <- rbind(fig_manifest_old, fig_manifest_new)
src_manifest <- rbind(src_manifest_old, src_manifest_new)
fig_manifest <- fig_manifest[!duplicated(fig_manifest$destination, fromLast = TRUE), ]
src_manifest <- src_manifest[!duplicated(src_manifest$destination, fromLast = TRUE), ]

write.csv(fig_manifest, fig_manifest_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(src_manifest, src_manifest_path, row.names = FALSE, fileEncoding = "UTF-8")

cat("Archived S1A cell-cycle supplement to:", fig_dest, "\n")
print(phase_summary[, c("subtype_short", "Phase", "n_cells", "cell_fraction")], row.names = FALSE)
