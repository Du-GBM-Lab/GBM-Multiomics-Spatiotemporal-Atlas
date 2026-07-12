suppressPackageStartupMessages({
  library(qs2)
  library(data.table)
  library(ggplot2)
})

root <- "/home/data/t010639/projects/GBM_R9_spatial_RCTD"
big_path <- file.path(root, "outputs/R9_A2_RCTD_weights_allslices_long.qs2")
out_dir <- file.path(root, "outputs/R9_C1_spatial_localization")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(root, "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

big <- as.data.table(qs2::qs_read(big_path))
needed <- c(
  "spot_id", "slice", "image", "x", "y",
  "Subtype1", "Subtype2", "Subtype3", "Subtype4",
  "Macrophages", "Microglial", "Monocytes",
  "Endothelial", "Mural cells"
)
stopifnot(all(needed %in% colnames(big)))

big[, `MES-V` := Subtype3]
big[, `MES-I` := Subtype4]
big[, `MES-lineage` := Subtype3 + Subtype4]
big[, Myeloid := Macrophages + Microglial + Monocytes]
big[, Vascular := Endothelial + `Mural cells`]

features <- c("MES-V", "MES-I", "MES-lineage", "Myeloid", "Vascular")
source_cols <- c("spot_id", "slice", "image", "x", "y", features)
source_dt <- big[, ..source_cols]
data.table::fwrite(source_dt, file.path(tab_dir, "R9_C1_spatial_localization_source.csv"))
qs2::qs_save(source_dt, file.path(out_dir, "R9_C1_spatial_localization_source.qs2"))

summary_dt <- big[, c(
  .(n_spots = .N),
  lapply(.SD, mean),
  setNames(lapply(.SD, max), paste0("max_", features))
), by = slice, .SDcols = features]
setorder(summary_dt, -`MES-V`)
data.table::fwrite(summary_dt, file.path(tab_dir, "R9_C1_perslice_feature_summary.csv"))

selected_slices <- c(
  "#UKF313_T_ST",  # high MES-V
  "#UKF304_T_ST",  # regular-size MES-I candidate
  "#UKF251_T_ST",  # regular-size MES-I candidate
  "#UKF265_T_ST"   # high MES-I but small slice caveat
)
selected_slices <- selected_slices[selected_slices %in% unique(big$slice)]

plot_one <- function(dt, feature, out_prefix) {
  max_val <- max(dt[[feature]], na.rm = TRUE)
  p <- ggplot(dt, aes(x = x, y = y, color = .data[[feature]])) +
    geom_point(size = 0.85, alpha = 0.95) +
    scale_y_reverse() +
    coord_fixed() +
    scale_color_gradientn(
      colours = c("#eeeeee", "#fee08b", "#f46d43", "#7f0000"),
      limits = c(0, max_val),
      name = feature
    ) +
    labs(title = paste0(unique(dt$slice), " | ", feature), x = NULL, y = NULL) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      legend.position = "right",
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  ggsave(paste0(out_prefix, ".pdf"), p, width = 4.2, height = 3.8, bg = "white")
  ggsave(paste0(out_prefix, ".png"), p, width = 4.2, height = 3.8, dpi = 300, bg = "white")
}

for (sl in selected_slices) {
  dt <- big[slice == sl]
  for (ft in features) {
    prefix <- file.path(
      fig_dir,
      paste0(gsub("[^A-Za-z0-9]+", "_", sl), "_", gsub("[^A-Za-z0-9]+", "_", ft))
    )
    plot_one(dt, ft, prefix)
  }
}

cat("== C1 source spots:", nrow(source_dt), "\n")
cat("== selected slices:\n")
print(selected_slices)
cat("== per-slice feature summary ordered by MES-V:\n")
print(summary_dt)
cat("== figures written:", length(list.files(fig_dir, pattern = "\\.(pdf|png)$")), "\n")
cat("\n[STOP C1] Spatial localization preview generated. Review figure candidates before PAN/niche or PLAU-PLAUR/FOSL1 co-enrichment.\n")
