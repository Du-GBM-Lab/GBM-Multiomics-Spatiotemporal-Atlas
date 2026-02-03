# Script: 08_Clinical_Validation.R
# Description: Generates Figure 4A (Heatmap) and 4C (Survival) using CGGA cohort.
# Dependencies: ComplexHeatmap, survival, survminer, dplyr, readr, circlize

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)
library(survival)
library(survminer)

# ================= 1. Data Loading & Preprocessing =================
# Paths
data_path <- "./data/hgg_data.rds"       
marker_path <- "./data/markers_subtypes_corl.csv"

cat("Loading data...\n")
gbm_data <- readRDS(data_path)
markers_data <- read_csv(marker_path)

# Helper: Prepare Matrix
prepare_matrix <- function(expr_df) {
  mat <- as.matrix(expr_df[, -1]) 
  rownames(mat) <- expr_df$Gene_Name
  return(mat)
}

expr_matrix <- prepare_matrix(gbm_data$expression)

# Filter Markers: Intersection with bulk data + Top 50 per subtype
available_genes <- intersect(markers_data$gene, rownames(expr_matrix))
markers_filtered <- markers_data %>%
  filter(gene %in% available_genes) %>%
  group_by(cluster) %>%
  top_n(50, wt = abs(avg_log2FC)) %>%
  ungroup()

# Subset Expression Matrix
selected_genes <- markers_filtered$gene
expr_subset <- expr_matrix[selected_genes, ]
# Z-score normalization for heatmap
scaled_expr <- t(scale(t(expr_subset)))

# ================= 2. Hierarchical Clustering =================
# Distance: 1 - Pearson Correlation
sample_dist <- as.dist(1 - cor(scaled_expr))
# Linkage: Ward.D2
sample_clust <- hclust(sample_dist, method = "ward.D2")

# Define 4 Clusters
clusters <- cutree(sample_clust, k = 4)
clinical_data <- gbm_data$clinical

# Match IDs and assign clusters
common_ids <- intersect(names(clusters), clinical_data$CGGA_ID)
clinical_data <- clinical_data[match(common_ids, clinical_data$CGGA_ID), ]
clinical_data$Cluster <- paste0("C", clusters[common_ids])

# ================= 3. Figure 4A: Annotated Heatmap =================
# Define Colors
col_clusters <- c("C1"="#FF7F00", "C2"="#FFD700", "C3"="#8A2BE2", "C4"="#00CED1")
col_subtypes <- c("Subtype_1"="#FFB6C1", "Subtype_2"="#ADD8E6", "Subtype_3"="#90EE90", "Subtype_4"="#DDA0DD")

# Row Annotation (Subtypes)
# Map gene list to subtype order
gene_split <- markers_filtered$cluster[match(rownames(scaled_expr), markers_filtered$gene)]
row_ha <- rowAnnotation(
  Subtype = gene_split,
  col = list(Subtype = col_subtypes),
  show_annotation_name = FALSE
)

# Column Annotation (Clinical Features)
col_ha <- HeatmapAnnotation(
  Cluster = clinical_data$Cluster,
  Grade = clinical_data$Grade,
  IDH = clinical_data$IDH_mutation_status,
  p19q = clinical_data$X1p19q_codeletion_status,
  MGMT = clinical_data$MGMTp_methylation_status,
  col = list(
    Cluster = col_clusters,
    Grade = c("WHO III"="#FF8C00", "WHO IV"="#B22222"),
    IDH = c("Mutant"="#FF9999", "Wildtype"="#66B2FF"),
    p19q = c("Codel"="#BA55D3", "Non-codel"="#FFA07A"),
    MGMT = c("methylated"="#6495ED", "un-methylated"="#CD5C5C")
  )
)

pdf("results/Fig4A_Heatmap.pdf", width = 12, height = 10)
Heatmap(scaled_expr[, common_ids],
        name = "Z-score",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        cluster_rows = FALSE, # Genes are pre-ordered by subtype
        cluster_columns = sample_clust,
        top_annotation = col_ha,
        left_annotation = row_ha,
        show_row_names = FALSE, show_column_names = FALSE)
dev.off()

# ================= 4. Figure 4C: Survival Analysis =================
cat("Running Survival Analysis...\n")
surv_obj <- Surv(clinical_data$OS, clinical_data$`Censor..alive.0..dead.1.`)
fit <- survfit(surv_obj ~ Cluster, data = clinical_data)

pdf("results/Fig4C_Survival.pdf", width = 8, height = 7)
ggsurvplot(fit, data = clinical_data,
           pval = TRUE, conf.int = TRUE,
           palette = col_clusters,
           risk.table = TRUE,
           title = "Overall Survival by Molecular Cluster",
           ggtheme = theme_classic())
dev.off()

cat("Analysis Complete. Figures 4A and 4C generated.\n")