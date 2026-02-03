################################################################################
# Script: 05_Cross_Species_Validation.R
# Description: Analysis of Mouse GBM dataset (Susek et al.).
#              Cross-species mapping, subtype annotation, and evolutionary analysis.
# Dependencies: Seurat, ggplot2, pheatmap, dplyr, RColorBrewer
################################################################################

library(Seurat)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# --- 1. Load Mouse Data ---
# Note: Assuming 'mouse_obj' is the raw integrated object
mouse_obj <- readRDS("data/raw/mouse_gbm_integrated.rds")

# --- 2. Subtype Annotation via Homologous Markers ---
message("Annotating Mouse Subtypes...")

# Define Homologous Markers (Mouse Symbols)
# Ordered: Astro -> Prolif -> Mesen -> Oligo
mouse_markers <- c(
  "Gap43", "Apoe", "Igfbp2", "Cxcl14", "Pdpn",      # Astro-sub
  "Top2a", "Mki67", "Nuf2", "Kifc1",                # Prolif-sub
  "Col1a1", "Col1a2", "Col8a1", "Pcolce",           # Mesen-sub
  "Mag", "Mog", "Mbp", "Plp1"                       # Oligo-sub
)

# Check for missing genes
mouse_markers <- mouse_markers[mouse_markers %in% rownames(mouse_obj)]

# Calculate Average Expression per Cluster (0-38)
# Using 'seurat_clusters' from original analysis
avg_expr <- AverageExpression(mouse_obj, features = mouse_markers, group.by = "seurat_clusters")$RNA
avg_expr_scaled <- t(scale(t(avg_expr)))

# Define Mapping (Manual mapping based on heatmap patterns - from your code)
# In a real pipeline, this map comes from interpreting the heatmap
cluster_to_subtype <- c(
  "0" = "Astro-sub", "1" = "Oligo-sub", "2" = "Prolif-sub", "3" = "Oligo-sub",
  "4" = "Oligo-sub", "5" = "Oligo-sub", "6" = "Mesen-sub", "7" = "Oligo-sub",
  "8" = "Astro-sub", "9" = "Mesen-sub", "10" = "Mesen-sub", "11" = "Prolif-sub",
  "12" = "Oligo-sub", "13" = "Prolif-sub", "14" = "Mesen-sub", "15" = "Astro-sub",
  "16" = "Oligo-sub", "17" = "Oligo-sub", "18" = "Astro-sub", "19" = "Prolif-sub",
  "20" = "Oligo-sub", "21" = "Prolif-sub", "22" = "Mesen-sub", "23" = "Mesen-sub",
  "24" = "Prolif-sub", "25" = "Mesen-sub", "26" = "Oligo-sub", "27" = "Prolif-sub",
  "28" = "Prolif-sub", "29" = "Oligo-sub", "30" = "Mesen-sub", "31" = "Prolif-sub",
  "32" = "Mesen-sub", "33" = "Mesen-sub", "34" = "Astro-sub", "35" = "Oligo-sub",
  "36" = "Prolif-sub", "37" = "Mesen-sub", "38" = "Prolif-sub"
)

mouse_obj$subtypes <- cluster_to_subtype[as.character(mouse_obj$seurat_clusters)]

# --- 3. Visualization: Heatmap (Figure S3C) ---
# Create annotation row for heatmap
anno_row <- data.frame(Subtype = factor(rep(c("Astro-sub", "Prolif-sub", "Mesen-sub", "Oligo-sub"), 
                                            times = c(5, 4, 4, 4)))) # Counts match 'mouse_markers' above roughly
rownames(anno_row) <- mouse_markers

pheatmap(avg_expr_scaled,
         cluster_rows = FALSE, cluster_cols = TRUE,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         annotation_row = anno_row,
         main = "Cross-species Marker Expression",
         filename = "figures/FigS3C_Mouse_Heatmap.pdf", width = 10, height = 8)

# --- 4. Visualization: UMAP (Figure 3D) ---
cols_subtypes <- c("Astro-sub" = "#1F77B4", "Prolif-sub" = "#FF7F0E", 
                   "Mesen-sub" = "#2CA02C", "Oligo-sub" = "#D62728")

p_umap <- DimPlot(mouse_obj, reduction = "umap", group.by = "subtypes", 
                  cols = cols_subtypes, pt.size = 0.5) +
  ggtitle("Mouse GBM Subtypes")

ggsave("figures/Fig3D_Mouse_UMAP.pdf", p_umap, width = 8, height = 6)

# --- 5. Evolutionary Analysis (Figure 3E) ---
message("Analyzing Subtype Evolution...")

# Extract metadata
meta <- mouse_obj@meta.data %>% select(sample, subtypes)

# Assign Disease Stage based on Sample Name
meta <- meta %>%
  mutate(stage = case_when(
    grepl("Preneoplastic", sample) ~ "Preneoplastic",
    grepl("Early_lesion", sample) ~ "Early lesion",
    grepl("Mid_lesion", sample) ~ "Mid lesion",
    grepl("Endpoint", sample) ~ "Endpoint"
  ))

# Calculate Proportions
props <- meta %>%
  group_by(stage, subtypes) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count))

# Factorize Stage for ordering
props$stage <- factor(props$stage, levels = c("Preneoplastic", "Early lesion", "Mid lesion", "Endpoint"))

# Stacked Area Plot
p_evol <- ggplot(props, aes(x = stage, y = proportion, fill = subtypes, group = subtypes)) +
  geom_area(alpha = 0.9) +
  scale_fill_manual(values = cols_subtypes) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Subtype Evolution", y = "Proportion", x = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/Fig3E_Evolution_Plot.pdf", p_evol, width = 8, height = 5)

saveRDS(mouse_obj, "data/processed/mouse_gbm_annotated.rds")