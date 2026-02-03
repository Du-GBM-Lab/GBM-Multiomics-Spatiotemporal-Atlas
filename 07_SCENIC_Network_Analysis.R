#Script: 07_SCENIC_Network_Analysis.py
#Description: Generates Figures 3B, 3C, 3D, and 3E.
#Dependencies: omicverse, scanpy, pandas, networkx, matplotlib, seaborn, gseapy


import os
import pandas as pd
import numpy as np
import scanpy as sc
import omicverse as ov
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from omicverse.external.pyscenic.rss import regulon_specificity_scores

# ================= Configuration =================
# Define subtype order and color palette
SUBTYPE_ORDER = ['Astro-sub', 'Prolif-sub', 'Oligo-sub', 'Mesen-sub']
COLORS = {
  'Mesen-sub': '#E63946', 
  'Prolif-sub': '#F4A261', 
  'Oligo-sub': '#457B9D', 
  'Astro-sub': '#1D3557'
}

# Define core driver TFs for each subtype (Used for Fig 3D & 3E)
CORE_DRIVERS = {
  'Mesen-sub': 'FOSL1', 
  'Astro-sub': 'ASCL1', 
  'Oligo-sub': 'SOX10', 
  'Prolif-sub': 'E2F1'
}

# ================= 1. Preprocessing & SCENIC Inference =================
print("Loading data and initializing SCENIC...")
adata = sc.read_h5ad("gbm_assembled.h5ad")

# Filter genes: min_cells=50, remove MT genes
sc.pp.filter_genes(adata, min_cells=50)
adata = adata[:, [g for g in adata.var_names if not g.startswith('MT-')]]

# Select Top 10,000 HVGs for robust network inference
sc.pp.highly_variable_genes(adata, n_top_genes=10000, flavor='seurat_v3', subset=True)

# Initialize OmicVerse SCENIC object
ov_obj = ov.single.SCENIC(
  adata=adata, 
  db_glob="./scenic_db/*feather", 
  motif_path="./scenic_db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl", 
  n_jobs=20
)

# Run GRN and Regulon prediction
edgelist = ov_obj.cal_grn(layer='counts')
ov_obj.cal_regulons()
auc_mtx = ov_obj.auc_mtx

# Calculate RSS (Regulon Specificity Score)
rss = regulon_specificity_scores(auc_mtx, adata.obs['cell_type'])

# ================= 2. Figure 3B: Regulon Activity Heatmap =================
print("Generating Figure 3B...")
# Subsample for balanced visualization (2000 cells total)
ad_subset = sc.pp.subsample(adata, n_obs=2000, copy=True)
ad_subset.obs['cell_type'] = ad_subset.obs['cell_type'].astype('category').cat.reorder_categories(SUBTYPE_ORDER)

# Select Top 5 specific regulons per subtype
top_regs = []
for ct in SUBTYPE_ORDER:
  top_regs.extend(rss[ct].sort_values(ascending=False).head(5).index.tolist())

sc.pl.heatmap(
  ad_subset, var_names=top_regs, groupby='cell_type', 
  cmap='RdYlBu_r', standard_scale='var', 
  categories_order=SUBTYPE_ORDER, figsize=(10, 6), save='_Fig3B.pdf'
)

# ================= 3. Figure 3C: RSS Specificity Bubble Plot =================
print("Generating Figure 3C...")
plot_data = []
for ct in SUBTYPE_ORDER:
  # Top 10 TFs per subtype
  top_10 = rss[ct].sort_values(ascending=False).head(10)
for tf, score in top_10.items():
  plot_data.append({'TF': tf, 'Subtype': ct, 'RSS': score})
df_rss = pd.DataFrame(plot_data)

plt.figure(figsize=(6, 8))
sns.scatterplot(data=df_rss, x='Subtype', y='TF', size='RSS', hue='Subtype', 
                palette=COLORS, sizes=(50, 400))
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.tight_layout()
plt.savefig("Fig3C_RSS_Bubble.pdf")

# ================= 4. Figure 3D: GO Enrichment for Core Drivers =================
print("Generating Figure 3D...")
go_res = []
for ct, tf in CORE_DRIVERS.items():
  # Get top 300 targets for the core driver
  targets = edgelist[edgelist['TF'] == tf].sort_values('importance', ascending=False).head(300)['target'].tolist()
# Run Enrichr
enr = gp.enrichr(gene_list=targets, gene_sets='GO_Biological_Process_2023', organism='Human').results
enr['TF'] = tf
go_res.append(enr.head(3)) # Keep top 3 terms

df_go = pd.concat(go_res)
plt.figure(figsize=(8, 6))
sns.barplot(data=df_go, y='Term', x='Combined Score', hue='TF', dodge=False, palette='viridis')
plt.savefig("Fig3D_GO_Enrichment.pdf", bbox_inches='tight')

# ================= 5. Figure 3E: Force-Directed Network =================
print("Generating Figure 3E...")
G = nx.Graph()
focus_genes = ['PLAUR', 'CD44', 'VIM', 'COL1A1']

for ct, tf in CORE_DRIVERS.items():
  tf_edges = edgelist[edgelist['TF'] == tf]
# Select Top 15 targets + any focused genes
targets = set(tf_edges.nlargest(15, 'importance')['target']).union(
  set(tf_edges[tf_edges['target'].isin(focus_genes)]['target']))

# Add nodes and edges
G.add_node(tf, type='TF', color=COLORS[ct], size=2000)
for t in targets:
  G.add_node(t, type='Target', color=COLORS[ct], size=300)
G.add_edge(tf, t, weight=0.5)

# Highlight PLAUR specifically
if 'PLAUR' in G.nodes:
  G.nodes['PLAUR']['color'] = '#FFD700' # Gold
G.nodes['PLAUR']['size'] = 600

pos = nx.spring_layout(G, k=0.5, seed=42)
node_colors = [nx.get_node_attributes(G, 'color')[n] for n in G.nodes()]
node_sizes = [nx.get_node_attributes(G, 'size')[n] for n in G.nodes()]

plt.figure(figsize=(10, 10))
nx.draw(G, pos, with_labels=True, node_color=node_colors, node_size=node_sizes, 
        edge_color='grey', alpha=0.8, font_size=8)
plt.savefig("Fig3E_Network.pdf")

# ================= Export Data for R =================
# Save FOSL1 targets for Spatial Analysis
fosl1_targets = edgelist[edgelist['TF']=='FOSL1'].sort_values('importance', ascending=False).head(50)['target'].tolist()
pd.Series(fosl1_targets).to_csv("FOSL1_targets_for_R.csv", index=False)
print("Done. FOSL1 targets exported.")


# Spatial_Analysis.R
# Description: Generates Figures 3F, 3G, and 3H.
# Dependencies: Seurat, CARD, ggplot2, dplyr, patchwork

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(CARD)

# ================= 1. Load Data & Define Signatures =================
# Load merged Seurat object (18 slides)
load("GBM_Spatial_Merged.RData") 

# Define Genes for PAN Niche (IVY-GAP signature)
pan_genes <- c("CA9", "VEGFA", "SLC2A1", "PGK1", "LDHA", "IGFBP2")

# Define Core Regulons for each subtype (For Fig 3G logic)
# FOSL1 targets are imported from Python SCENIC output
fosl1_targets <- read.csv("FOSL1_targets_for_R.csv")$x
regulon_lists <- list(
  FOSL1_Activity = fosl1_targets,                          # Mesen-sub Driver
  ASCL1_Activity = c("ASCL1", "DLL3", "SOX4", "TCF4"),     # Astro-sub Driver
  SOX10_Activity = c("SOX10", "MBP", "PLP1", "MAG"),       # Oligo-sub Driver
  E2F1_Activity  = c("E2F1", "MKI67", "TOP2A", "PCNA")     # Prolif-sub Driver
)

# ================= 2. Score Calculation & Deconvolution =================
# Calculate PAN Score
st_obj <- AddModuleScore(st_obj, features = list(IVY_PAN = pan_genes), name = "IVY_PAN_Score")

# Calculate Regulon Activities
for(reg_name in names(regulon_lists)) {
  st_obj <- AddModuleScore(st_obj, features = list(regulon_lists[[reg_name]]), name = reg_name)
}

# Clean column names (Seurat adds '1' suffix)
colnames(st_obj@meta.data) <- gsub("IVY_PAN_Score1", "IVY_PAN_Score", colnames(st_obj@meta.data))
colnames(st_obj@meta.data) <- gsub("_Activity1", "_Activity", colnames(st_obj@meta.data))

# Run CARD Deconvolution (Pseudocode, assuming run previously)
# Output columns expected: "prop_Mesen_sub", "prop_Astro_sub", etc.
# st_obj <- RunCARD(...) 

# ================= 3. Figure 3F: Spatial Congruency Maps =================
# Select a representative sample showing PAN regions
sample_id <- "Sample_UKF248"
subset_obj <- subset(st_obj, subset = sample_id == sample_id)

p1 <- SpatialFeaturePlot(subset_obj, features = "prop_Mesen_sub", pt.size.factor = 1.6) + 
  scale_fill_viridis_c(option = "magma", name="Mesen Prop") + ggtitle("Mesen-sub Proportion")
p2 <- SpatialFeaturePlot(subset_obj, features = "FOSL1_Activity", pt.size.factor = 1.6) + 
  scale_fill_viridis_c(option = "magma", name="FOSL1 Act") + ggtitle("FOSL1 Regulon Activity")
p3 <- SpatialFeaturePlot(subset_obj, features = "IVY_PAN_Score", pt.size.factor = 1.6) + 
  scale_fill_viridis_c(option = "cividis", name="PAN Score") + ggtitle("IVY-GAP PAN Niche")

ggsave("Fig3F_Spatial_Congruency.pdf", p1 | p2 | p3, width = 15, height = 5)

# ================= 4. Figure 3G: Spatial Correlation Boxplot =================
# Logic: Calculate correlation of PAN Score vs. Cell Prop AND PAN Score vs. Specific Driver
# Note: FOSL1 is paired with Mesen, ASCL1 with Astro, etc.

sample_ids <- unique(st_obj$sample_id)
plot_df <- data.frame()

# Mapping: Subtype Name -> c(Cell Proportion Column, Specific Regulon Column)
mapping <- list(
  "Mesen"  = c("prop_Mesen_sub",  "FOSL1_Activity"),
  "Astro"  = c("prop_Astro_sub",  "ASCL1_Activity"),
  "Oligo"  = c("prop_Oligo_sub",  "SOX10_Activity"),
  "Prolif" = c("prop_Prolif_sub", "E2F1_Activity")
)

for(s in sample_ids) {
  sub_meta <- st_obj@meta.data[st_obj$sample_id == s, ]
  
  for(type in names(mapping)) {
    cols <- mapping[[type]]
    prop_col <- cols[1]
    reg_col <- cols[2]
    
    # 1. Correlation: PAN vs Cell Proportion
    r_prop <- cor(sub_meta[[prop_col]], sub_meta$IVY_PAN_Score, method="pearson")
    # 2. Correlation: PAN vs Specific Regulon Activity
    r_reg  <- cor(sub_meta[[reg_col]], sub_meta$IVY_PAN_Score, method="pearson")
    
    plot_df <- rbind(plot_df, data.frame(
      Sample = s, Subtype = type, 
      Metric = "Cell Proportion", Correlation = r_prop
    ))
    plot_df <- rbind(plot_df, data.frame(
      Sample = s, Subtype = type, 
      Metric = "Regulon Activity", Correlation = r_reg
    ))
  }
}

# Plot Boxplot
p_box <- ggplot(plot_df, aes(x = Subtype, y = Correlation, fill = Metric)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha=0.5, size=1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  scale_fill_manual(values = c("Cell Proportion"="#455A64", "Regulon Activity"="#E64A19")) +
  theme_classic() +
  labs(y = "Spatial Similarity to PAN Niche (Pearson r)", x = "")

ggsave("Fig3G_Correlation_Boxplot.pdf", p_box, width = 8, height = 5)

# ================= 5. Figure 3H: Global Scatter Density =================
# Global correlation between FOSL1 Activity and PAN Score across all spots
df_scatter <- st_obj@meta.data[, c("FOSL1_Activity", "IVY_PAN_Score")]
r_val <- round(cor(df_scatter$FOSL1_Activity, df_scatter$IVY_PAN_Score), 3)

p_scatter <- ggplot(df_scatter, aes(x = IVY_PAN_Score, y = FOSL1_Activity)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", color="white") +
  scale_fill_viridis_c(option = "magma") +
  geom_smooth(method = "lm", color = "yellow", linetype="dashed") +
  annotate("text", x = min(df_scatter$IVY_PAN_Score), y = max(df_scatter$FOSL1_Activity), 
           label = paste0("R = ", r_val), hjust=0, vjust=1, size=6, color="black") +
  theme_bw() + labs(title = "Global Spatial Synchronization")

ggsave("Fig3H_Scatter_Density.pdf", p_scatter, width = 6, height = 6)