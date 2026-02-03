################################################################################
# Script: 02_Functional_Enrichment_and_Scoring.R
# Description: 1. GO/KEGG Enrichment for interpretation.
#              2. Single-cell Hallmark scoring using irGSEA (Multi-method).
################################################################################

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(irGSEA)
library(Seurat)
library(RColorBrewer)

# --- Load Data ---
tumor_obj <- readRDS("data/processed/GBM_malignant_processed.rds")
markers_df <- read_csv("results/tables/markers_subtypes_corl.csv")

# ==============================================================================
# Part 1: Traditional Enrichment (GO/KEGG)
# ==============================================================================
# --- 2. Load Marker Data ---
# This input comes from the subtyping step (Script 03/04)
# Ensure the file path is correct relative to the repo root
reginal_subtype_markers <- read_csv("results/tables/markers_subtypes_corl.csv")

# --- 3. Define Enrichment Helper Function ---
# Encapsulates the logic for filtering, ID conversion, and enrichment
perform_enrichment <- function(markers_df, subtype) {
  
  message(paste0("Processing enrichment for: ", subtype))
  
  # Filter robust markers: strict logFC > 0.5 for cleaner pathways
  subtype_markers <- markers_df %>%
    filter(cluster == subtype,
           p_val_adj < 0.05,
           abs(avg_log2FC) > 0.5)
  
  if(nrow(subtype_markers) == 0) return(NULL)
  
  # Convert Symbol to Entrez ID (Required for KEGG/GO)
  gene_ids <- bitr(subtype_markers$gene,
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db)
  
  # KEGG Pathway Enrichment
  kegg_res <- enrichKEGG(
    gene          = gene_ids$ENTREZID,
    organism      = "hsa",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.25
  ) 
  
  # Set readable gene names for KEGG results
  if(!is.null(kegg_res)) {
    kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  }
  
  # GO Biological Process Enrichment
  go_res <- enrichGO(
    gene          = gene_ids$ENTREZID,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.25,
    readable      = TRUE
  )
  
  return(list(kegg = kegg_res, go = go_res))
}

# --- 4. Execute Analysis for All Subtypes ---
subtypes <- unique(reginal_subtype_markers$cluster)
enr_list <- lapply(subtypes, function(st) perform_enrichment(reginal_subtype_markers, st))
names(enr_list) <- subtypes

# --- 5. Aggregate Results for Visualization ---
# Extract and combine KEGG results
all_kegg <- bind_rows(lapply(names(enr_list), function(st) {
  res <- enr_list[[st]]$kegg
  if (!is.null(res) && nrow(res@result) > 0) {
    res@result %>% mutate(Subtype = st)
  }
}))

# Extract and combine GO results
all_go <- bind_rows(lapply(names(enr_list), function(st) {
  res <- enr_list[[st]]$go
  if (!is.null(res) && nrow(res@result) > 0) {
    res@result %>% mutate(Subtype = st)
  }
}))

# --- 6. Select Top Terms for Plotting ---
top_n <- 5

# Process KEGG Data
top_kegg <- all_kegg %>%
  group_by(Subtype) %>%
  arrange(p.adjust) %>%
  slice_head(n = top_n) %>%
  ungroup() %>%
  mutate(
    negLogP = -log10(p.adjust),
    # Wrap long text for better plotting
    Description = stringr::str_wrap(Description, width = 40),
    Description = factor(Description, levels = unique(Description))
  )

# Process GO Data
top_go <- all_go %>%
  group_by(Subtype) %>%
  arrange(p.adjust) %>%
  slice_head(n = top_n) %>%
  ungroup() %>%
  mutate(
    negLogP = -log10(p.adjust),
    Description = stringr::str_wrap(Description, width = 40),
    Description = factor(Description, levels = unique(Description))
  )

# --- 7. Visualization: KEGG Barplot (Fig 2C Right) ---
p_kegg <- ggplot(top_kegg, aes(x = Description, y = negLogP, fill = Subtype)) +
  geom_col(width = 0.7, color = "black", size = 0.2) +
  facet_grid(Subtype ~ ., scales = "free_y", space = "free_y") +
  coord_flip() +
  scale_fill_manual(values = brewer.pal(max(3, length(subtypes)), "Set1")) +
  labs(title = "KEGG Pathway Enrichment",
       x = NULL, y = expression(-log[10](adj.P))) +
  theme_bw() +
  theme(
    strip.text.y       = element_text(angle = 0, size = 10, face = "bold"),
    axis.text.y        = element_text(size = 10, color = "black"),
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    legend.position    = "none",
    panel.grid.major.y = element_blank()
  )

ggsave("figures/Fig2C_KEGG_Barplot.pdf", p_kegg, width = 7, height = 8)

# --- 8. Visualization: GO Barplot (Fig 2C Left) ---
p_go <- ggplot(top_go, aes(x = Description, y = negLogP, fill = Subtype)) +
  geom_col(width = 0.7, color = "black", size = 0.2) +
  facet_grid(Subtype ~ ., scales = "free_y", space = "free_y") +
  coord_flip() +
  scale_fill_manual(values = brewer.pal(max(3, length(subtypes)), "Set2")) +
  labs(title = "GO Biological Process Enrichment",
       x = NULL, y = expression(-log[10](adj.P))) +
  theme_bw() +
  theme(
    strip.text.y       = element_text(angle = 0, size = 10, face = "bold"),
    axis.text.y        = element_text(size = 10, color = "black"),
    plot.title         = element_text(hjust = 0.5, face = "bold"),
    legend.position    = "none",
    panel.grid.major.y = element_blank()
  )

ggsave("figures/Fig2C_GO_Barplot.pdf", p_go, width = 7, height = 8)

# ==============================================================================
# Part 2: Single-cell Hallmark Scoring (irGSEA)
# ==============================================================================

message("Running irGSEA scoring...")

# 1. Prepare Gene Sets (MSigDB Hallmark)
# Note: irGSEA has built-in MSigDB support, but we define custom list for control
# or use msigdb = T directly.
# We perform scoring on a subset if dataset is too large (>50k cells)
if(ncol(tumor_obj) > 30000) {
  set.seed(123)
  obj_scoring <- subset(tumor_obj, cells = sample(colnames(tumor_obj), 20000))
} else {
  obj_scoring <- tumor_obj
}

# 2. Run Multiple Scoring Methods
# We use 'UCell', 'ssGSEA', 'AUCell' for robust consensus
obj_scoring <- irGSEA.score(
  object = obj_scoring,
  assay = "RNA",
  slot = "data",
  msigdb = TRUE,
  species = "Homo sapiens",
  category = "H",  # Hallmark
  method = c("UCell", "ssGSEA", "AUCell"),
  kcdf = 'Gaussian' 
)

# 3. Integrate Scores (Robust Rank Aggregation)
# This calculates the RRA p-value for each pathway across methods
result_irGSEA <- irGSEA.integrate(
  object = obj_scoring,
  group.by = "ident", # Use the subtypes defined in Step 03
  method = c("UCell", "ssGSEA", "AUCell")
)

# 4. Visualization (Bubble Plot - Figure 2E style)
# Showing RRA score and P-value
p_bubble <- irGSEA.bubble(
  object = result_irGSEA,
  method = "RRA",
  direction = "up",
  top = 10
) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/Fig2E_Hallmark_Bubble.pdf", p_bubble, width = 8, height = 10)

# 5. Visualization (Heatmap - Figure 2D style)
p_heatmap <- irGSEA.heatmap(
  object = result_irGSEA,
  method = "RRA",
  top = 20,
  show_row_names = TRUE
)

ggsave("figures/Fig2D_Hallmark_Heatmap.pdf", p_heatmap, width = 10, height = 8)

saveRDS(result_irGSEA, "data/processed/irGSEA_results.rds")