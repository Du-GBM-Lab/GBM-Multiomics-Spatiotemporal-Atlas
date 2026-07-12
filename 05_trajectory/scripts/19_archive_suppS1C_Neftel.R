options(stringsAsFactors = FALSE)

root <- normalizePath("<DATA_ROOT>/项目/分型/修稿杠生信/重新分析", winslash = "/", mustWork = TRUE)
archive <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/04_拟时序分析"

dir.create(file.path(archive, "补充图"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(archive, "源数据", "补图S1C_Neftel状态"), recursive = TRUE, showWarnings = FALSE)

fig_src <- file.path(root, "06_恶性细胞拟时序/figures/final_panels/supplement/S1C_Neftel_AMS.pdf")
fig_dest <- file.path(archive, "补充图", "补图S1C_Neftel状态.pdf")
file.copy(fig_src, fig_dest, overwrite = TRUE)

src_data_dir <- file.path(root, "06_恶性细胞拟时序/figures/final_panels/source_data/supplement")
source_file <- file.path(src_data_dir, "01_Neftel_6state_AMS.csv")
dest_file <- file.path(archive, "源数据", "补图S1C_Neftel状态", basename(source_file))
file.copy(source_file, dest_file, overwrite = TRUE)

ams <- read.csv(source_file, check.names = FALSE)
ams_summary <- do.call(rbind, lapply(split(ams, paste(ams$module, ams$subtype_short, sep = "|")), function(x) {
  data.frame(
    module = x$module[1],
    subtype_short = x$subtype_short[1],
    n_cells = nrow(x),
    mean = mean(x$score, na.rm = TRUE),
    median = median(x$score, na.rm = TRUE),
    pct_pos = mean(x$score > 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
ams_summary <- ams_summary[order(ams_summary$module, ams_summary$subtype_short), ]
summary_dest <- file.path(archive, "源数据", "补图S1C_Neftel状态", "补图S1C_Neftel_AMS摘要.csv")
write.csv(ams_summary, summary_dest, row.names = FALSE, fileEncoding = "UTF-8")

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
  source = c(source_file, summary_dest),
  destination = c(dest_file, summary_dest),
  timestamp = now,
  stringsAsFactors = FALSE
)

fig_manifest <- rbind(fig_manifest_old, fig_manifest_new)
src_manifest <- rbind(src_manifest_old, src_manifest_new)
fig_manifest <- fig_manifest[!duplicated(fig_manifest$destination, fromLast = TRUE), ]
src_manifest <- src_manifest[!duplicated(src_manifest$destination, fromLast = TRUE), ]

write.csv(fig_manifest, fig_manifest_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(src_manifest, src_manifest_path, row.names = FALSE, fileEncoding = "UTF-8")

cat("Archived S1C Neftel supplement to:", fig_dest, "\n")
print(ams_summary, row.names = FALSE)
