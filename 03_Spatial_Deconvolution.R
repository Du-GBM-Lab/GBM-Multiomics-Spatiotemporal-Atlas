################################################################################
# Script: 03_Spatial_Deconvolution.R
# Description: Deconvoluting Visium data using scRNA-seq reference via CARD.
#              Mapping subtype proportions to spatial locations (Fig S2E).
# Dependencies: CARD, Seurat, ggplot2, Matrix
################################################################################

library(CARD)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(patchwork)

# --- 1. Load Data ---
# Single-cell reference (Processed in Step 03)
sc_data <- readRDS("data/processed/GBM_malignant_processed.rds")
# Spatial data (Merged object containing all slices)
ST_merge <- readRDS("data/raw/ST_merge.rds")

# Ensure memory is sufficient for large matrix operations
options(future.globals.maxSize = 40000 * 1024^2)

# --- 2. Prepare Reference Data (Run Once) ---
message("Preparing Single-cell Reference for CARD...")

# Intersect genes between SC and ST
common_genes <- intersect(rownames(ST_merge), rownames(sc_data))
message(paste("Number of common genes:", length(common_genes)))

# Prepare SC Count Matrix and Metadata
sc_counts <- GetAssayData(sc_data, assay = "RNA", slot = "counts")[common_genes, ]
# Use the annotation column defined in previous steps (e.g., 'subtypes_corl' or 'ident')
sc_meta <- data.frame(
  cellType = sc_data@meta.data$anno_ident, # Or your specific subtype column
  sample = "sc_reference",
  row.names = colnames(sc_counts)
)

# --- 3. Define Deconvolution Function ---
# Encapsulates the logic for coordinate extraction and CARD execution
run_card_per_sample <- function(sample_id, full_st_obj, sc_count_mat, sc_meta_df, gene_set) {
  
  message(paste(">>> Processing Sample:", sample_id))
  
  # A. Subset ST Object
  sample_data <- subset(full_st_obj, orig.ident == sample_id)
  sp_counts <- GetAssayData(sample_data, assay = "Spatial", slot = "counts")[gene_set, ]
  
  # B. Extract Coordinates (Robust Handling)
  image_name <- names(sample_data@images)[1]
  sp_location_df <- NULL
  
  # Try standard Seurat coordinates
  try({
    coords <- GetTissueCoordinates(sample_data, image = image_name)
    sp_location_df <- data.frame(x = coords[,1], y = coords[,2])
    rownames(sp_location_df) <- colnames(sp_counts)
  }, silent = TRUE)
  
  # Fallback: Manual extraction or Grid generation (From your original code)
  if (is.null(sp_location_df)) {
    message("   Using fallback coordinate extraction...")
    n_spots <- ncol(sp_counts)
    grid_size <- ceiling(sqrt(n_spots))
    sp_location_df <- data.frame(
      x = rep(1:grid_size, each = grid_size)[1:n_spots],
      y = rep(1:grid_size, times = grid_size)[1:n_spots]
    )
    rownames(sp_location_df) <- colnames(sp_counts)
  }
  
  # C. Run CARD
  card_obj <- createCARDObject(
    sc_count = as.matrix(sc_count_mat),
    sc_meta = sc_meta_df,
    spatial_count = as.matrix(sp_counts),
    spatial_location = sp_location_df,
    ct.varname = "cellType",
    ct.select = unique(sc_meta_df$cellType),
    sample.varname = "sample"
  )
  
  card_res <- CARD_deconvolution(card_obj)
  return(card_res@Proportion_CARD)
}

# --- 4. Execute Deconvolution Loop ---
sample_list <- unique(ST_merge$orig.ident)
# List to store results temporarily
deconv_results <- list()

for (s in sample_list) {
  # Skip samples with too few spots if necessary
  tryCatch({
    props <- run_card_per_sample(s, ST_merge, sc_counts, sc_meta, common_genes)
    deconv_results[[s]] <- props
  }, error = function(e) {
    message(paste("Error in sample", s, ":", e$message))
  })
  gc() # Clean up memory
}

# --- 5. Merge Results Back to Seurat Object ---
message("Merging CARD proportions into Spatial Metadata...")

for (s in names(deconv_results)) {
  props <- deconv_results[[s]]
  cells <- rownames(props)
  
  # Normalize names if needed and add to metadata
  for (ctype in colnames(props)) {
    # Create metadata column name: e.g., "CARD_Mesen-sub"
    col_name <- paste0("CARD_", ctype)
    ST_merge@meta.data[cells, col_name] <- props[, ctype]
  }
}

# Save the intermediate object with deconvolution results
saveRDS(ST_merge, "data/processed/ST_merge_with_CARD.rds")

# --- 6. Visualization (Figure S2E Style) ---
# Generate Spatial Feature Plots for the 4 subtypes
message("Generating Spatial Feature Plots...")

subtypes_to_plot <- c("CARD_Astro-sub", "CARD_Prolif-sub", "CARD_Mesen-sub", "CARD_Oligo-sub")
pdf("figures/FigS2E_Spatial_Subtypes.pdf", width = 12, height = 10)

for (s in sample_list) {
  obj_sub <- subset(ST_merge, orig.ident == s)
  
  # Check if metadata exists
  if (all(subtypes_to_plot %in% colnames(obj_sub@meta.data))) {
    p_list <- list()
    for (feat in subtypes_to_plot) {
      p <- SpatialFeaturePlot(obj_sub, features = feat, pt.size.factor = 1.6, image.alpha = 0.5) +
        scale_fill_viridis_c(option = "magma") +
        ggtitle(paste(s, "-", feat)) +
        theme(legend.position = "right")
      p_list[[feat]] <- p
    }
    print(wrap_plots(p_list, ncol = 2))
  }
}
dev.off()

message("Script 06 Completed Successfully.")