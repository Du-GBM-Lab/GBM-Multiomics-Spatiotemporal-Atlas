################################################################################
# Script: 12_Spatial_CellChat_Analysis.R
# Description: Spatially-resolved ligand-receptor analysis.
#              Focus on Global Network and PLAU-PLAUR Specificity.
# Dependencies: Seurat, CellChat, dplyr, ggplot2, patchwork, ComplexHeatmap
################################################################################

library(Seurat)
library(CellChat)
library(tidyverse)
library(ComplexHeatmap)
library(patchwork)

# Set parallel processing
future::plan("multisession", workers = 4)

# --- 1. Data Loading & Preparation ---
# Load the merged spatial object with CARD predictions
st_merged <- readRDS("data/processed/ST_merge_with_CARD.rds")
st_samples <- SplitObject(st_merged, split.by = "orig.ident")

# --- 2. Batch Spatial CellChat ---
cellchat_list <- list()

for (sample_name in names(st_samples)) {
  tryCatch({
    message(paste("Running Spatial CellChat:", sample_name))
    obj <- st_samples[[sample_name]]
    
    # Define Major Cell Type from CARD proportions
    card_cols <- grep("CARD_", colnames(obj@meta.data), value = T)
    props <- obj@meta.data[, card_cols]
    major_types <- colnames(props)[max.col(props, ties.method = "first")]
    obj$major_celltype <- gsub("CARD_", "", major_types)
    Idents(obj) <- "major_celltype"
    
    # Prepare Input
    data.input <- GetAssayData(obj, slot = "data", assay = "SCT")
    meta <- obj@meta.data
    meta$labels <- Idents(obj)
    
    # Extract Spatial Coordinates (Crucial for Spatial Mode)
    # Using the standard image coordinates
    img_name <- names(obj@images)[1]
    spatial.locs <- GetTissueCoordinates(obj, image = img_name)[, 1:2]
    colnames(spatial.locs) <- c("x", "y")
    
    # Scale factors (Standard Visium)
    scale.factors <- list(spot.diameter = 65, spot = 65)
    
    # Create Object with Spatial Mode
    cc <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                         datatype = "spatial", 
                         coordinates = as.matrix(spatial.locs),
                         scale.factors = scale.factors)
    
    cc@DB <- CellChatDB.human
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    
    # Compute Probability with Distance Constraint
    cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1, 
                            distance.use = TRUE, interaction.range = 250) # 250 micron radius
    cc <- filterCommunication(cc, min.cells = 10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    
    cellchat_list[[sample_name]] <- cc
    
  }, error = function(e) { message(paste("Skip:", sample_name)) })
}

# Merge for group analysis
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))

# --- 3. Visualization: Global Network (Fig 7A, 7D) ---
# Aggregated Heatmap (Fig 7A)
pdf("figures/Fig7A_Interaction_Strength.pdf", width = 8, height = 7)
netVisual_heatmap(cellchat_merged, measure = "weight", color.heatmap = "Reds")
dev.off()

# Bubble Plot: Mesen-sub -> Macrophage Mechanisms (Fig 7D)
# Focusing on SPP1, FN1, etc.
pdf("figures/Fig7D_Mesen_to_Macro_Mechanisms.pdf", width = 10, height = 8)
netVisual_bubble(cellchat_merged, sources.use = "Mesen-sub", targets.use = "Macrophages", 
                 remove.isolate = FALSE)
dev.off()

# --- 4. PLAU-PLAUR Specificity Analysis (Fig 7F, 7G) ---
message("Analyzing PLAU-PLAUR Axis...")

# Extract high-sensitivity interactions for PLAU
target_pair <- "PLAU - PLAUR"
plau_senders <- data.frame()

for (n in names(cellchat_list)) {
  # Relaxed threshold to capture all potential sources
  tmp <- subsetCommunication(cellchat_list[[n]], thresh = 1) 
  if(!is.null(tmp)) {
    tmp_plau <- tmp %>% filter(interaction_name_2 == target_pair)
    if(nrow(tmp_plau) > 0) {
      tmp_plau$Sample <- n
      plau_senders <- rbind(plau_senders, tmp_plau)
    }
  }
}

# Boxplot of PLAU Senders (Fig 7F)
# Corrected logic: Showing Macrophages can also be senders
p_sender <- ggplot(plau_senders, aes(x = source, y = prob, fill = source)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
  labs(title = "PLAU Ligand Source Probability", y = "Interaction Probability") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/Fig7F_PLAU_Senders.pdf", p_sender, width = 6, height = 5)

# Interaction Matrix (Fig 7G)
# Visualizing who talks to whom via PLAU
pdf("figures/Fig7G_PLAU_Matrix.pdf", width = 6, height = 5)
netVisual_circle(cellchat_merged@net$weight, signaling = "PLAU", vertex.weight = 10, arrow.width = 1)
dev.off()

saveRDS(cellchat_merged, "data/processed/Spatial_CellChat_Final.rds")