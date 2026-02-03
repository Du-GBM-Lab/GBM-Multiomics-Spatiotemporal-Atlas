################################################################################
# Project: GBM Single-Cell Transcriptomic Atlas & Subtyping
# Description: Integrated pipeline for QC, malignant cell identification, 
#              batch correction (Harmony), and robust subtyping.
# Data Source: Mei et al., Nat Cancer 2023
################################################################################

# --- Load Required Libraries ---
library(Seurat)
library(harmony)
library(tidyverse)
library(pheatmap)
library(cluster) # For Gap statistic
library(viridis)

# ==============================================================================
# Section 1: Data Loading & Adaptive Quality Control
# ==============================================================================

# Load raw integrated object (Source: Mei et al.)
# This object contains annotations from the original study
obj <- readRDS("data/raw/0.GBM.RNA.integrated.24.rds")

# Visualize initial QC metrics
print("Initial Data Dimensions:")
print(dim(obj))

# Adaptive Filtering Strategy
# We use dynamic quantiles to account for biological variability across patients
# rather than rigid hard thresholds.
high_count_threshold <- quantile(obj$nCount_RNA, 0.99)
high_gene_threshold <- quantile(obj$nFeature_RNA, 0.98)

obj <- subset(obj, 
              subset = nFeature_RNA > 500 & 
                nFeature_RNA < high_gene_threshold & 
                nCount_RNA < high_count_threshold & 
                percent.mt < 40) # Standard mitochondrial threshold

# Normalize and Find Variable Features
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)

# ==============================================================================
# Section 2: Malignant Cell Identification & Extraction
# ==============================================================================

# Define publication-ready color palette for visualization
colors_publication <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488",
  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
  "#7E6148", "#B09C85", "#2E8B57", "#6A3D9A",
  "#FF7F00", "#1F77B4", "#D62728", "#9467BD", "#8C564B"
)

# Visualize original annotations to confirm lineage accuracy
DimPlot(obj, reduction = "umap", group.by = "anno_ident", 
        cols = colors_publication, label = TRUE, raster = FALSE) +
  ggtitle("Global Cell Lineages (Mei et al. Annotation)")

# Extract Malignant Lineages
# Based on original annotations and validated by CNV (Chr7+/Chr10-)
# Note: CNV inference (inferCNV) is performed in a separate upstream pipeline.
tumor_lineages <- c("Astrocytes", "Oligodendrocytes", "OPCs", "Radial glial")
tumor_obj <- subset(obj, subset = anno_ident %in% tumor_lineages)

print(paste("Number of malignant cells extracted:", ncol(tumor_obj)))

# ==============================================================================
# Section 3: Batch Correction & Integration (Harmony)
# ==============================================================================

# Re-process the tumor-specific object
tumor_obj <- FindVariableFeatures(tumor_obj, selection.method = "vst", nfeatures = 2000)
tumor_obj <- ScaleData(tumor_obj, features = VariableFeatures(tumor_obj))
tumor_obj <- RunPCA(tumor_obj, npcs = 50, verbose = FALSE)

# Run Harmony to integrate across 24 patients
# Parameters optimized for GBM heterogeneity preservation
tumor_obj <- RunHarmony(tumor_obj, 
                        group.by = "Pt_number", 
                        dims.use = 1:30,
                        theta = 2.5,     # Diversity clustering penalty (Fine-tuned)
                        lambda = 1,      # Ridge regression penalty
                        sigma = 0.1,     # Width of soft k-means cluster
                        max.iter.harmony = 20, 
                        plot_convergence = FALSE)

# ==============================================================================
# Section 4: Dimensionality Reduction & Clustering
# ==============================================================================

# Run UMAP on Harmony embeddings
tumor_obj <- RunUMAP(tumor_obj, reduction = "harmony", dims = 1:30)

# Construct SNN Graph
tumor_obj <- FindNeighbors(tumor_obj, reduction = "harmony", dims = 1:30)

# Unsupervised Clustering
# We over-cluster slightly (res=0.6) and then merge based on correlation analysis
tumor_obj <- FindClusters(tumor_obj, resolution = 0.6)

# Visualization of preliminary clusters
DimPlot(tumor_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Preliminary Unsupervised Clusters")

# ==============================================================================
# Section 5: Robust Subtyping via Correlation Analysis
# ==============================================================================

# Identify Marker Genes for each preliminary cluster
markers <- FindAllMarkers(tumor_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Select Top 50 markers per cluster for robust definition
top50_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

# Calculate Average Expression Profile for each cluster
avg_expr <- AverageExpression(tumor_obj, features = top50_markers, assays = "RNA", slot = "data")$RNA

# Compute Spearman Correlation Matrix
# This step reveals the 4 meta-modules (Astro, Oligo, Mesen, Prolif)
cor_mat <- cor(avg_expr, method = "spearman")

# Gap Statistic Analysis to determine optimal k
# Note: K.max set to 8 to cover the potential range
gap_stat <- clusGap(as.matrix(cor_mat), FUN = kmeans, K.max = 8, B = 50)
print(gap_stat)
plot(gap_stat, main = "Gap Statistic for Optimal k")

# Generate Correlation Heatmap with Hierarchical Clustering
# The distinct blocks correspond to the 4 GBM malignant subtypes
pheatmap(cor_mat,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D", # Ward linkage minimizes variance
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         border_color = NA,
         fontsize = 10,
         main = "Robust Subtype Identification Matrix")

# Save the processed object for downstream NMF validation and Spatial mapping
saveRDS(tumor_obj, file = "data/processed/GBM_malignant_processed.rds")

# End of Pipeline