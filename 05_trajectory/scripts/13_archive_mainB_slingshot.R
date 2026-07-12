options(stringsAsFactors = FALSE)

root <- normalizePath("<DATA_ROOT>/项目/分型/修稿杠生信/重新分析", winslash = "/", mustWork = TRUE)
archive <- "<DATA_ROOT>/项目/分型/修稿杠生信/图片表格/04_拟时序分析"

dir.create(file.path(archive, "正文图"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(archive, "源数据", "正图B_三分支轨迹"), recursive = TRUE, showWarnings = FALSE)

fig_src <- file.path(root, "06_恶性细胞拟时序/figures/final_panels/main/Main_B_slingshot_curves.pdf")
fig_dest <- file.path(archive, "正文图", "正图B_三分支轨迹.pdf")
file.copy(fig_src, fig_dest, overwrite = TRUE)

src_data_dir <- file.path(root, "06_恶性细胞拟时序/figures/final_panels/source_data/main")
table_dir <- file.path(root, "06_恶性细胞拟时序/tables")

source_files <- c(
  file.path(src_data_dir, "Main_B_umap_cells.csv"),
  file.path(src_data_dir, "Main_B_slingshot_curves.csv"),
  file.path(src_data_dir, "Main_B_curve_start.csv"),
  file.path(src_data_dir, "Main_B_curve_endpoints.csv"),
  file.path(table_dir, "slingshot_branch_points.csv"),
  file.path(table_dir, "lineage_subtype_crosstab.csv"),
  file.path(table_dir, "slingshot_lineage_summary.csv")
)

dest_files <- file.path(archive, "源数据", "正图B_三分支轨迹", basename(source_files))
for (i in seq_along(source_files)) {
  file.copy(source_files[i], dest_files[i], overwrite = TRUE)
}

crosstab <- read.csv(file.path(table_dir, "lineage_subtype_crosstab.csv"), check.names = FALSE)
lineage_total <- aggregate(n_cells ~ lineage_id, data = crosstab, sum)
names(lineage_total)[2] <- "lineage_n_cells"
crosstab <- merge(crosstab, lineage_total, by = "lineage_id")
crosstab$lineage_percent <- 100 * crosstab$n_cells / crosstab$lineage_n_cells
crosstab$lineage_label <- c(lineage1 = "L1: NPC-P to MES-I", lineage2 = "L2: NPC-P to MES-V", lineage3 = "L3: NPC-P to OPC-M")[crosstab$lineage_id]
crosstab <- crosstab[order(crosstab$lineage_id, -crosstab$n_cells), ]

summary_dest <- file.path(archive, "源数据", "正图B_三分支轨迹", "正图B_lineage组成摘要.csv")
write.csv(crosstab, summary_dest, row.names = FALSE, fileEncoding = "UTF-8")

branch <- read.csv(file.path(table_dir, "slingshot_branch_points.csv"), check.names = FALSE)
branch_dest <- file.path(archive, "源数据", "正图B_三分支轨迹", "正图B_branch_point摘要.csv")
write.csv(branch, branch_dest, row.names = FALSE, fileEncoding = "UTF-8")

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
  out <- tryCatch(
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
  out
}

fig_manifest_old <- read_manifest(fig_manifest_path)
src_manifest_old <- read_manifest(src_manifest_path)

main_a_pdf <- file.path(archive, "正文图", "正图A_NPCP根_CytoTRACE2.pdf")
if (file.exists(main_a_pdf) && !main_a_pdf %in% fig_manifest_old$destination) {
  fig_manifest_old <- rbind(
    fig_manifest_old,
    data.frame(
      category = "main",
      type = "figure",
      source = file.path(root, "06_恶性细胞拟时序/figures/final_panels/main/Main_A_root_CytoTRACE2.pdf"),
      destination = main_a_pdf,
      timestamp = now,
      stringsAsFactors = FALSE
    )
  )
}

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
  source = c(source_files, summary_dest, branch_dest),
  destination = c(dest_files, summary_dest, branch_dest),
  timestamp = now,
  stringsAsFactors = FALSE
)

fig_manifest <- rbind(fig_manifest_old, fig_manifest_new)
src_manifest <- rbind(src_manifest_old, src_manifest_new)
fig_manifest <- fig_manifest[!duplicated(fig_manifest$destination, fromLast = TRUE), ]
src_manifest <- src_manifest[!duplicated(src_manifest$destination, fromLast = TRUE), ]

write.csv(fig_manifest, fig_manifest_path, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(src_manifest, src_manifest_path, row.names = FALSE, fileEncoding = "UTF-8")

cat("Archived Main B to:", fig_dest, "\n")
cat("Lineage totals:\n")
print(lineage_total)
cat("Lineage composition:\n")
print(crosstab[, c("lineage_id", "lineage_label", "subtype_short", "n_cells", "lineage_n_cells", "lineage_percent")])
cat("Branch point:\n")
print(branch)
