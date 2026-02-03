################################################################################
# Script: 04_Trajectory_Analysis.R
# Description: Pseudotime trajectory inference using Monocle3.
#              Identification of root state and gene expression kinetics.
# Dependencies: Seurat, monocle3, tidyverse, patchwork
################################################################################

library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)

# --- 1. Load Data ---
# Load the processed Seurat object from Step 03
obj <- readRDS("data/processed/GBM_malignant_processed.rds")

# --- 2. Initialize Monocle3 CDS ---
message("Converting Seurat object to Monocle3 CDS...")

# Extract data
data <- GetAssayData(obj, assay = 'RNA', slot = 'counts')
cell_metadata <- obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# Create CDS
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# Transfer UMAP coordinates from Seurat to CDS
# This ensures the trajectory overlays perfectly with Seurat's UMAP
umap_coords <- Embeddings(obj, reduction = "umap")
reducedDims(cds)$UMAP <- umap_coords

# Transfer Clustering Info
# Note: Using 'subtypes_corl' (renamed to Astro/Mesen/etc. in preprocessing)
clusters <- obj$subtypes_corl 
names(clusters) <- colnames(obj)
cds@clusters$UMAP$clusters <- clusters
cds@clusters$UMAP$partitions <- clusters # Use subtypes as partitions

# --- 3. Learn Trajectory Graph ---
message("Learning Principal Graph...")
# Key parameters from your code
cds <- learn_graph(cds, use_partition = FALSE, close_loop = TRUE)

# --- 4. Order Cells (Pseudotime) ---
message("Ordering cells...")

# Function to automatically find the root node in Astro-sub
get_earliest_principal_node <- function(cds, root_subtype = "Astro-sub"){
  cell_ids <- which(colData(cds)[, "subtypes_corl"] == root_subtype)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # Find the node with the most Astro-sub cells
  root_node <- names(which.max(table(closest_vertex[cell_ids,])))
  return(root_node)
}

root_node <- get_earliest_principal_node(cds, root_subtype = "Astro-sub")
cds <- order_cells(cds, root_pr_nodes = root_node)

# --- 5. Visualization (Figure 3A, 3B) ---
# Trajectory by Subtype
p_traj <- plot_cells(cds,
                     color_cells_by = "subtypes_corl",
                     label_groups_by_cluster = FALSE,
                     label_leaves = FALSE,
                     label_branch_points = TRUE,
                     graph_label_size = 3) +
  ggtitle("Trajectory by Subtype")

# Trajectory by Pseudotime
p_pseudo <- plot_cells(cds,
                       color_cells_by = "pseudotime",
                       label_cell_groups = FALSE,
                       label_leaves = FALSE,
                       label_branch_points = FALSE) + 
  scale_color_viridis_c(option = "plasma") +
  ggtitle("Pseudotime")

ggsave("figures/Fig3A_Trajectory.pdf", p_traj, width = 8, height = 6)
ggsave("figures/Fig3B_Pseudotime.pdf", p_pseudo, width = 8, height = 6)

# --- 6. Gene Expression Kinetics (Figure 3C) ---
message("Plotting Gene Kinetics...")

# Define markers (Ensure they match dataset gene symbols)
genes_to_plot <- c("GAP43", "GPM6A", "IGFBP2",   # Early (Astro)
                   "COL1A1", "ISLR",             # Mid (Mesen)
                   "TOP2A", "MKI67",             # Late (Prolif)
                   "MAG", "MBP")                 # Late (Oligo)

# Filter for available genes
genes_avail <- intersect(genes_to_plot, rowData(cds)$gene_short_name)
cds_subset <- cds[genes_avail, ]

p_genes <- plot_genes_in_pseudotime(cds_subset,
                                    color_cells_by = "subtypes_corl",
                                    min_expr = 0.5,
                                    ncol = 3)

ggsave("figures/Fig3C_Gene_Kinetics.pdf", p_genes, width = 12, height = 8)

saveRDS(cds, "data/processed/monocle3_cds_final.rds")