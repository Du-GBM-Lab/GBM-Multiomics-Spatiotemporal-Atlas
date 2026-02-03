################################################################################
# Script: 10_hdWGCNA_PLAUR_Network.R
# Description: High-dimensional WGCNA to identify Mesen-sub specific modules
#              and validate PLAUR as a hub gene (Figure 5 & S5).
# Dependencies: Seurat, hdWGCNA, WGCNA, tidyverse, clusterProfiler, ggplot2
################################################################################

library(Seurat)
library(tidyverse)
library(WGCNA)
library(hdWGCNA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)

# --- 1. Load Data & Preprocessing ---
obj <- readRDS("data/processed/GBM_malignant_processed.rds")

# Ensure proper naming for visualization
# (Hidden detail: This mapping ensures consistency with manuscript figures)
obj@meta.data$subtypes_corl <- recode(
  obj@meta.data$subtypes_corl,
  "Subtype_1" = "Astro-sub", "Subtype_2" = "Prolif-sub",
  "Subtype_3" = "Mesen-sub", "Subtype_4" = "Oligo-sub"
)

# --- 2. Setup hdWGCNA ---
message("Initializing hdWGCNA...")

obj <- SetupForWGCNA(
  obj,
  gene_select = "fraction", 
  fraction = 0.1, # Adjusted parameter for robust gene selection
  wgcna_name = "GBM_Malignant_Network"
)

# Construct Metacells to reduce sparsity
# Parameter Tweaking: k=30 used here for higher robustness than default
obj <- MetacellsByGroups(
  obj,
  group.by = "subtypes_corl",
  k = 30, 
  max_shared = 15,
  ident.group = "subtypes_corl"
)

obj <- NormalizeMetacells(obj)

# --- 3. Network Construction ---
message("Constructing Co-expression Network...")

# Set expression data
obj <- SetDatExpr(
  obj,
  group_name = unique(obj$subtypes_corl),
  group.by = "subtypes_corl",
  assay = "RNA", slot = "data"
)

# Soft Power Selection
# We simulate the selection process here. 
# In practice, we selected power based on Scale Free Topology Model Fit > 0.8
# obj <- TestSoftPowers(obj, networkType = "signed")
# power_table <- GetPowerTable(obj)
selected_power <- 12 # Hardcoded based on optimization

# Construct Network
obj <- ConstructNetwork(
  obj,
  soft_power = selected_power,
  setDatExpr = FALSE,
  tom_name = "Malignant_TOM",
  overwrite_tom = TRUE
)

# Dendrogram Visualization (Figure 5A)
pdf("figures/Fig5A_WGCNA_Dendrogram.pdf", width = 8, height = 6)
PlotDendrogram(obj, main = "Malignant Co-expression Modules")
dev.off()

# --- 4. Module Analysis & kME Calculation ---
message("Calculating Module Eigengenes and Connectivity...")

# Scale data for PCA
obj <- ScaleData(obj, features = VariableFeatures(obj))

# Calculate MEs
obj <- ModuleEigengenes(obj, group.by.vars = "subtypes_corl")

# Calculate kME (Module Connectivity)
obj <- ModuleConnectivity(
  obj, 
  group.by = "subtypes_corl", 
  group_name = unique(obj$subtypes_corl)
)

# --- 5. Identify PLAUR Module & Hub Status (Figure S5E) ---
modules <- GetModules(obj)
target_gene <- "PLAUR"

# Find which module contains PLAUR
target_module <- modules %>% filter(gene_name == target_gene) %>% pull(module)

if(length(target_module) > 0) {
  message(paste(target_gene, "belongs to module:", target_module))
  
  # Extract kME for the target module
  kme_col <- paste0("kME_", target_module)
  module_genes <- modules %>% filter(module == target_module) %>% arrange(desc(get(kme_col)))
  
  # Visualization: kME Barplot (Waterfall plot)
  # Highlight PLAUR position
  module_genes$Rank <- 1:nrow(module_genes)
  module_genes$IsTarget <- ifelse(module_genes$gene_name == target_gene, "Yes", "No")
  
  p_kme <- ggplot(module_genes, aes(x = Rank, y = get(kme_col), fill = IsTarget)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = c("Yes" = "red", "No" = "grey80")) +
    labs(title = paste("Hub Gene Analysis:", target_module, "Module"),
         y = "kME (Connectivity)", x = "Gene Rank") +
    theme_classic() + theme(legend.position = "none")
  
  ggsave("figures/FigS5E_kME_Waterfall.pdf", p_kme, width = 6, height = 4)
  
} else {
  warning(paste(target_gene, "not found in any module."))
}

# --- 6. Functional Enrichment (Figure 5C) ---
message("Running Enrichment Analysis...")

if(length(target_module) > 0) {
  genes_use <- modules %>% filter(module == target_module) %>% pull(gene_name)
  genes_entrez <- bitr(genes_use, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  
  # GO Enrichment
  go_res <- enrichGO(
    gene = genes_entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05
  )
  
  # Visualization
  p_go <- dotplot(go_res, showCategory = 10) + 
    ggtitle(paste("GO Enrichment:", target_module, "Module"))
  
  ggsave("figures/Fig5C_Module_GO.pdf", p_go, width = 7, height = 6)
}

# --- 7. Subtype Specificity of Module (Figure 5B/S5G) ---
message("Visualizing Module Specificity...")

if(length(target_module) > 0) {
  # Get Module Eigengene values
  me_col <- paste0("ME", target_module)
  plot_data <- obj@meta.data %>% select(subtypes_corl, all_of(me_col))
  colnames(plot_data)[2] <- "ME_Score"
  
  # Violin Plot with Boxplot overlay
  p_vln <- ggplot(plot_data, aes(x = subtypes_corl, y = ME_Score, fill = subtypes_corl)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, fill = "white") +
    scale_fill_brewer(palette = "Set1") +
    labs(title = paste(target_module, "Module Expression"), y = "Eigengene Expression") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("figures/Fig5B_Module_Specificity.pdf", p_vln, width = 5, height = 4)
}

saveRDS(obj, "data/processed/GBM_hdWGCNA_Result.rds")