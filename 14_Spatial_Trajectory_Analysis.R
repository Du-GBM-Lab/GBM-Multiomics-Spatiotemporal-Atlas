################################################################################
# Script: 13_Spatial_Trajectory_Analysis.R
# Description: Spatial Monocle3 analysis.
#              1. Auto-Root identification (Astro-score).
#              2. Spatial Trajectory Mapping (Fig 7I-K).
#              3. Gene Expression Dynamics (Fig 7L-M).
# Dependencies: Seurat, monocle3, ggplot2, viridis, ComplexHeatmap
################################################################################

library(Seurat)
library(monocle3)
library(tidyverse)
library(viridis)
library(ComplexHeatmap)
library(circlize)

# --- 1. Load Single Sample (Example Workflow) ---
# In practice, loop this over the samples (e.g., #UKF313, #UKF242)
st_obj <- readRDS("data/processed/ST_sample_UKF313.rds")

# --- 2. Tumor Cell Isolation & Subtyping ---
message("Isolating Tumor Cells...")

# Define signatures for scoring
sigs <- list(
  Astro = c("GPM6A", "GAP43", "IGFBP2"),
  Mesen = c("PLAU", "PLAUR", "VIM", "MMP9"),
  Oligo = c("MAG", "MBP", "PLP1"),
  Prolif = c("TOP2A", "MKI67")
)
st_obj <- AddModuleScore(st_obj, features = sigs, name = "Score_")
colnames(st_obj@meta.data)[grep("Score_", colnames(st_obj@meta.data))] <- c("Score_Astro", "Score_Mesen", "Score_Oligo", "Score_Prolif")

# Filter Tumor Cells (Malignant Score threshold)
st_obj$malignant_score <- rowSums(st_obj@meta.data[, c("Score_Astro", "Score_Mesen", "Score_Oligo", "Score_Prolif")])
tumor_cells <- WhichCells(st_obj, expression = malignant_score > 0.4)

# Create Tumor-Only Object (Avoids V5 subset bugs)
tumor_counts <- GetAssayData(st_obj, assay = "Spatial", layer = "counts")[, tumor_cells]
tumor_obj <- CreateSeuratObject(counts = tumor_counts, meta.data = st_obj@meta.data[tumor_cells, ])

# Unsupervised Clustering
tumor_obj <- NormalizeData(tumor_obj) %>% ScaleData() %>% RunPCA() %>% FindNeighbors() %>% FindClusters(resolution = 0.4)

# --- 3. Monocle 3 Trajectory Construction ---
message("Constructing Spatial Trajectory...")

# Convert to CDS
cds <- new_cell_data_set(
  GetAssayData(tumor_obj, layer = "counts"),
  cell_metadata = tumor_obj@meta.data,
  gene_metadata = data.frame(gene_short_name = rownames(tumor_obj), row.names = rownames(tumor_obj))
)

# Learn Graph
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE)

# --- 4. Auto-Root Identification (Fig 7H) ---
# Identify the cluster with highest Astro-sub score
cluster_scores <- tumor_obj@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(Mean_Astro = mean(Score_Astro))

root_cluster <- as.character(cluster_scores$seurat_clusters[which.max(cluster_scores$Mean_Astro)])
message(paste("Auto-detected Root Cluster:", root_cluster))

# Helper to find root node
get_root_node <- function(cds, target_cluster, cluster_col="seurat_clusters"){
  cells <- rownames(colData(cds))[colData(cds)[[cluster_col]] == target_cluster]
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex)
  names(which.max(table(closest_vertex[cells, ])))
}

root_node <- get_root_node(cds, root_cluster)
cds <- order_cells(cds, root_pr_nodes = root_node)

# --- 5. Spatial Visualization (Fig 7I-K) ---
# Map pseudotime back to Seurat for spatial plotting
tumor_obj$Pseudotime <- pseudotime(cds)
st_obj$Pseudotime_Plot <- NA
st_obj$Pseudotime_Plot[colnames(tumor_obj)] <- tumor_obj$Pseudotime

# Plot: Spatial Pseudotime (Fig 7J)
p_time <- SpatialFeaturePlot(st_obj, features = "Pseudotime_Plot", pt.size.factor = 2.5, stroke = 0) +
  scale_fill_viridis_c(option = "plasma", na.value = "transparent") +
  ggtitle("Spatial Trajectory (Pseudotime)")

ggsave("figures/Fig7J_Spatial_Pseudotime.pdf", p_time, width = 8, height = 6)

# --- 6. Gene Dynamics (Fig 7L, 7M) ---
message("Plotting Gene Dynamics...")

target_genes <- c("GPM6A", "SOX2", "TOP2A", "VIM", "FN1", "SPP1", "PLAU", "PLAUR", "VEGFA", "CA9")
valid_genes <- intersect(target_genes, rownames(cds))
cds_subset <- cds[valid_genes, ]

# Waterfall Heatmap (Fig 7L)
# (Simplified matrix preparation for heatmap)
pt_matrix <- normalized_counts(cds_subset)
pt_order <- order(pseudotime(cds_subset))
pt_matrix <- log1p(pt_matrix[, pt_order])
pt_matrix <- t(scale(t(as.matrix(pt_matrix))))

pdf("figures/Fig7L_Waterfall_Heatmap.pdf", width = 8, height = 10)
Heatmap(pt_matrix, cluster_columns = FALSE, show_column_names = FALSE, 
        name = "Expression", col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
dev.off()

# Trend Plots (Fig 7M)
p_trend <- plot_genes_in_pseudotime(cds_subset, color_cells_by = "pseudotime", 
                                    min_expr = 0.1, ncol = 5) +
  scale_color_viridis_c(option = "plasma")

ggsave("figures/Fig7M_Gene_Trends.pdf", p_trend, width = 15, height = 8)

message("Spatial Trajectory Analysis Complete.")