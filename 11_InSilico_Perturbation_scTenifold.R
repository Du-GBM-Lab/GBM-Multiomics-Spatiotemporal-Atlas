################################################################################
# Script: 11_InSilico_Perturbation_scTenifold.R
# Description: Virtual Knockout of PLAUR in Mesen-sub using scTenifoldKnk.
#              1. Data preparation (Subsetting & Downsampling).
#              2. Network construction & Virtual KO.
#              3. Manifold Alignment & Z-score calculation.
#              4. Functional Enrichment of perturbed genes (Fig 5D, 5E).
# Dependencies: Seurat, scTenifoldKnk, clusterProfiler, ggplot2, ggrepel
################################################################################

library(Seurat)
library(scTenifoldKnk)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(dplyr)

# --- 1. Data Preparation ---
message("Preparing data for virtual knockout...")

# Load processed Seurat object
obj <- readRDS("data/processed/GBM_malignant_processed.rds")

# Ensure subtype names are correct
# Assuming 'subtypes_corl' contains "Subtype_3" or "Mesen-sub"
# Here we standardize to "Mesen-sub" for clarity
if("Subtype_3" %in% obj$subtypes_corl) {
  target_cells <- WhichCells(obj, expression = subtypes_corl == "Subtype_3")
} else {
  target_cells <- WhichCells(obj, expression = subtypes_corl == "Mesen-sub")
}

mes_obj <- subset(obj, cells = target_cells)

# Downsample to 2000 cells for efficiency (Standard practice for scTenifold)
set.seed(42)
if(ncol(mes_obj) > 2000){
  mes_obj <- subset(mes_obj, cells = sample(colnames(mes_obj), 2000))
}

# Select Features: HVGs + Target Gene
mes_obj <- FindVariableFeatures(mes_obj, nfeatures = 2000)
hvg_genes <- VariableFeatures(mes_obj)
target_gene <- "PLAUR"

if(!target_gene %in% hvg_genes){
  final_genes <- unique(c(target_gene, hvg_genes))
} else {
  final_genes <- hvg_genes
}

# Extract Count Matrix
# scTenifoldKnk expects raw counts
counts_matrix <- GetAssayData(mes_obj, assay = "RNA", slot = "counts")[final_genes, ]
counts_matrix <- as.matrix(counts_matrix)

# --- 2. Run Virtual Knockout ---
message(paste("Running scTenifoldKnk for:", target_gene))

# Note: This step is computationally intensive.
# nc_nNet = 50 for robust results in final publication
res_ko <- scTenifoldKnk(
  countMatrix = counts_matrix, 
  gKO = target_gene, 
  nc_nNet = 50, 
  nCores = 8 # Adjust based on available cores
)

# Save intermediate results
saveRDS(res_ko, "data/processed/scTenifold_PLAUR_KO.rds")

# --- 3. Result Processing ---
# Extract differential regulation table
dr_table <- res_ko$diffRegulation
dr_table <- dr_table[order(dr_table$Z, decreasing = TRUE), ]

# Save full list
write.csv(dr_table, "results/tables/PLAUR_KO_Perturbation_List.csv", row.names = FALSE)

# --- 4. Visualization: Ranked Z-Score Plot (Figure 5E) ---
message("Plotting Ranked Gene Disruption...")

dr_table$rank <- 1:nrow(dr_table)
# Label top genes for visualization
dr_table$label <- ifelse(dr_table$rank <= 10 | dr_table$gene == target_gene, dr_table$gene, "")

p_rank <- ggplot(dr_table, aes(x = rank, y = Z)) +
  geom_point(aes(color = Z > 1.2), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("grey70", "#D55E00")) +
  geom_text_repel(aes(label = label), max.overlaps = 50, 
                  fontface = "bold", box.padding = 0.5, segment.color = "grey50") +
  geom_hline(yintercept = 1.2, linetype = "dashed", color = "blue") +
  theme_classic(base_size = 14) +
  labs(title = paste("Virtual KO Effect:", target_gene, "in Mesen-sub"),
       x = "Gene Rank", y = "Disruption Score (Z-score)") +
  theme(legend.position = "none")

ggsave("figures/Fig5E_Ranked_ZScore.pdf", p_rank, width = 8, height = 6)

# --- 5. Functional Enrichment (Figure 5F) ---
message("Running Functional Enrichment on Perturbed Genes...")

# Select significantly perturbed genes
sig_genes <- dr_table$gene[dr_table$Z > 1.2]

# GO Enrichment
ego <- enrichGO(
  gene          = sig_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP", 
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# Visualization
p_enrich <- dotplot(ego, showCategory=15) + 
  ggtitle(paste("GO Enrichment (Perturbed Genes Z>1.2)")) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("figures/Fig5F_Perturbation_Enrichment.pdf", p_enrich, width = 9, height = 7)

message("Virtual Knockout Analysis Complete.")