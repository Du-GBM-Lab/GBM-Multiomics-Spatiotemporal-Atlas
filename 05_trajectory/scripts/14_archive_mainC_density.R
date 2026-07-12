options(stringsAsFactors = FALSE)

root <- normalizePath("<DATA_ROOT>/项目/分型/修稿杠生信/重新分析", winslash = "/", mustWork = TRUE)
archive <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/04_拟时序分析"

dir.create(file.path(archive, "正文图"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(archive, "源数据", "正图C_拟时序密度"), recursive = TRUE, showWarnings = FALSE)

fig_src <- file.path(root, "06_恶性细胞拟时序/figures/final_panels/main/Main_C_lineage_density.pdf")
fig_dest <- file.path(archive, "正文图", "正图C_拟时序密度.pdf")
file.copy(fig_src, fig_dest, overwrite = TRUE)

src_data_dir <- file.path(root, "06_恶性细胞拟时序/figures/final_panels/source_data/main")
source_files <- file.path(src_data_dir, c(
  "Main_C_lineage_pseudotime_density.csv",
  "Main_C_lineage_pseudotime_density_plotted.csv",
  "Main_C_lineage_density_counts.csv"
))
dest_files <- file.path(archive, "源数据", "正图C_拟时序密度", basename(source_files))
for (i in seq_along(source_files)) {
  file.copy(source_files[i], dest_files[i], overwrite = TRUE)
}

d <- read.csv(file.path(src_data_dir, "Main_C_lineage_pseudotime_density.csv"), check.names = FALSE)
summary <- do.call(rbind, lapply(split(d, paste(d$lineage_id, d$subtype_short, sep = "|")), function(x) {
  data.frame(
    lineage_id = x$lineage_id[1],
    lineage_label = x$lineage_label[1],
    subtype_short = x$subtype_short[1],
    n_cells = nrow(x),
    median_pseudotime = median(x$pseudotime, na.rm = TRUE),
    q25_pseudotime = as.numeric(quantile(x$pseudotime, 0.25, na.rm = TRUE)),
    q75_pseudotime = as.numeric(quantile(x$pseudotime, 0.75, na.rm = TRUE)),
    min_pseudotime = min(x$pseudotime, na.rm = TRUE),
    max_pseudotime = max(x$pseudotime, na.rm = TRUE)
  )
}))
summary <- summary[order(summary$lineage_id, summary$median_pseudotime), ]
summary_dest <- file.path(archive, "源数据", "正图C_拟时序密度", "正图C_pseudotime分布摘要.csv")
write.csv(summary, summary_dest, row.names = FALSE, fileEncoding = "UTF-8")

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
  source = c(source_files, summary_dest),
  destination = c(dest_files, summary_dest),
  timestamp = now,
  stringsAsFactors = FALSE
)

fig_manifest <- rbind(fig_manifest_old, fig_manifest_new)
src_manifest <- rbind(src_manifest_old, src_manifest_new)
fig_manifest <- fig_manifest[!duplicated(fig_manifest$destination, fromLast = TRUE), ]
src_manifest <- src_manifest[!duplicated(src_manifest$destination, fromLast = TRUE), ]

write.csv(fig_manifest, fig_manifest_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(src_manifest, src_manifest_path, row.names = FALSE, fileEncoding = "UTF-8")

cat("Archived Main C to:", fig_dest, "\n")
print(summary[, c("lineage_id", "subtype_short", "n_cells", "median_pseudotime", "q25_pseudotime", "q75_pseudotime")], row.names = FALSE)
