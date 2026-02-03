################################################################################
# Script: 06_scATAC_Motif_and_Footprinting.R
# Description: scATAC-seq Label Transfer, ChromVAR Motif Analysis, 
#              and Footprinting (Fig 2L, 3A, S3D, S3E).
# Dependencies: Signac, Seurat, JASPAR2020, TFBSTools, BSgenome.Hsapiens.UCSC.hg38,
#               chromVAR, pheatmap, ggplot2
################################################################################

library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(chromVAR)
library(pheatmap)
library(dplyr)

# --- 1. Label Transfer from scRNA-seq (Fig 2L) ---
# Load data (Assuming pre-processed objects)
sc_data <- readRDS("data/processed/GBM_malignant_processed.rds")
atac_data <- readRDS("data/raw/scATAC_GBM.rds")

DefaultAssay(atac_data) <- "RNA" # Use gene activity assay for integration

message("Transferring labels...")
transfer.anchors <- FindTransferAnchors(
  reference = sc_data,
  query = atac_data,
  reference.assay = "RNA",
  query.assay = "RNA",
  reduction = "cca"
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = sc_data@meta.data[["subtypes_corl"]], # Use our robust subtypes
  weight.reduction = atac_data[["lsi"]],
  dims = 2:50
)

atac_data <- AddMetaData(atac_data, metadata = predicted.labels)

# Visualization
p_dim <- DimPlot(atac_data, group.by = "predicted.id", label = FALSE) + 
  ggtitle("scATAC-seq with Predicted Labels")
ggsave("figures/Fig2L_ATAC_UMAP.pdf", p_dim, width = 8, height = 6)

# --- 2. ChromVAR Motif Analysis (Fig 3A, S3D) ---
message("Running ChromVAR...")

DefaultAssay(atac_data) <- "ATAC"
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)

atac_data <- AddMotifs(atac_data, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
atac_data <- RunChromVAR(atac_data, genome = BSgenome.Hsapiens.UCSC.hg38)

# Differential Motif Activity
# Focus on Astro-sub specific motifs as per manuscript
diff_motifs_astro <- FindMarkers(
  object = atac_data,
  ident.1 = "Astro-sub",
  only.pos = TRUE,
  test.use = "wilcox",
  mean.fxn = rowMeans,
  fc.name = "avg_diff",
  assay = 'chromvar'
)

# Prepare Heatmap Data (Diagonal Alignment)
# Using the specific motif list derived in your analysis
focus_motifs <- c("NFIC", "NFIB", "ASCL1", "Tcf12", "Hic1",  # Astro
                  "JUNB", "FOSL1", "JUND", "FOS::JUND",      # Mesen
                  "SOX8", "SOX9", "SOX15", "SOX13", "SOX2",  # Oligo
                  "LMX1B", "MNX1", "ARGFX", "VAX1")          # Prolif

# Helper function to get ID from Name
motif_ids <- sapply(focus_motifs, function(n) {
  names(which(sapply(pfm, name) == n))[1]
})
motif_ids <- unlist(motif_ids)

# Generate Heatmap
p_heatmap <- DoHeatmap(
  object = atac_data,
  features = motif_ids,
  group.by = "predicted.id",
  assay = 'chromvar',
  slot = 'data'
) + scale_fill_gradient2(low = "blue", mid = "white", high = "red")

ggsave("figures/Fig3A_Motif_Activity_Heatmap.pdf", p_heatmap, width = 10, height = 8)

# --- 3. Footprinting Analysis (Fig S3E) ---
message("Generating Footprints...")

# Perform footprinting on key motifs
atac_data <- Footprint(
  object = atac_data,
  motif.name = focus_motifs,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# Plotting Loop
pdf("figures/FigS3E_Footprints.pdf", width = 8, height = 6)
for (motif in focus_motifs) {
  print(PlotFootprint(atac_data, features = motif, group.by = "predicted.id"))
}
dev.off()

saveRDS(atac_data, "data/processed/scATAC_analyzed.rds")


################################################################################
# Script: MultiOmics_CoveragePlots.R
# Description: Linked visualization of Chromatin Accessibility and Gene Expression
#              (Figure S4A).
# Dependencies: Signac, Seurat, ggplot2, patchwork
################################################################################

library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)

# --- 1. Load Data ---
atac_data <- readRDS("data/processed/scATAC_analyzed.rds")
# Ensure RNA assay is present or linked for ExpressionPlot
# Assuming 'atac_data' has an 'RNA' assay from integration or loaded separately

# --- 2. Define Targets ---
genes_of_interest <- c(
  "NFIC", "NFIB", "ASCL1", "Tcf12", "Hic1",
  "JUNB", "FOSL1", "JUND", "FOS::JUND",
  "SOX8", "SOX9", "SOX15", "SOX13", "SOX2",
  "LMX1B", "mix-a", "MNX1", "ARGFX", "VAX1"
)

# --- 3. Generate Combined Plots ---
output_dir <- "figures/multiomics_plots"
dir.create(output_dir, showWarnings = FALSE)

tile_idents <- c("Astro-sub", "Mesen-sub", "Oligo-sub", "Prolif-sub")

for (gene in genes_of_interest) {
  tryCatch({
    message(paste("Plotting:", gene))
    
    # 1. Coverage Plot (ATAC Peaks)
    cov_plot <- CoveragePlot(
      object = atac_data,
      region = gene,
      annotation = TRUE,
      peaks = TRUE,
      links = FALSE
    )
    
    # 2. Tile Plot (Fragment Density)
    tile_plot <- TilePlot(
      object = atac_data,
      region = gene,
      idents = tile_idents,
      tile.size = 200
    )
    
    # 3. Expression Plot (RNA Levels)
    # Handle gene name differences (e.g., 'FOS::JUND' -> 'FOS')
    gene_clean <- unlist(strsplit(gene, "::"))[1] 
    
    expr_plot <- ExpressionPlot(
      object = atac_data,
      features = gene_clean,
      assay = "RNA"
    )
    
    # 4. Combine
    combined_plot <- (cov_plot | expr_plot) / tile_plot + 
      plot_layout(heights = c(3, 2), guides = 'collect')
    
    ggsave(file.path(output_dir, paste0(gene, "_multiomics.png")), 
           combined_plot, width = 12, height = 10)
    
  }, error = function(e) {
    message(paste("Error plotting", gene, ":", e$message))
  })
}

message("Multi-omics visualization completed.")