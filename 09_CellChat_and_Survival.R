# Script: 09_CellChat_and_Survival.R
# Description: Generates Figure 4D-L and Figure S5 (Communication & PLAU Axis).
# Dependencies: CellChat, Seurat, ggplot2, ComplexHeatmap, survival, survminer, patchwork

library(CellChat)
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(patchwork)
library(survival)
library(survminer)
library(dplyr)

# ================= 1. Data Prep & CellChat Object =================
# Load Seurat object (user provided)
# obj <- readRDS("path/to/seurat_object.rds")

# Ensure clean annotations
obj_filtered <- subset(obj, anno_ident != "Ambiguous")
obj_filtered$anno_ident <- factor(obj_filtered$anno_ident)

# Define groups
tumor_subtypes <- c("Astro-sub", "Mesen-sub", "Oligo-sub", "Prolif-sub")
all_idents <- levels(obj_filtered$anno_ident)
tme_cells <- setdiff(all_idents, tumor_subtypes)

# Initialize CellChat
data.input <- GetAssayData(obj_filtered, assay = "RNA", slot = "data")
meta <- data.frame(group = obj_filtered$anno_ident, row.names = colnames(obj_filtered))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

# DB & Processing
cellchat@DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, nboot = 100)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# ================= 2. Figure 4D & 4E: Global Network =================
pdf("results/Fig4D_Network_Circle.pdf", width = 8, height = 8)
netVisual_circle(cellchat@net$weight, weight.scale = TRUE, label.edge=FALSE, 
                 title.name = "Global Interaction Strength")
dev.off()

pdf("results/Fig4E_Network_Heatmap.pdf", width = 8, height = 6)
netVisual_heatmap(cellchat, title.name = "Aggregated Interaction Strength")
dev.off()

# ================= 3. Figure 4F & S5A: Specific LR Chords =================
# Fig 4F: Focused Chord (Mesen-sub Outgoing)
pdf("results/Fig4F_Mesen_Chord.pdf", width = 10, height = 10)
netVisual_chord_gene(cellchat, sources.use = "Mesen-sub", targets.use = tme_cells, 
                     lab.cex = 0.8, title.name = "Signals from Mesen-sub")
dev.off()

# Fig S5A: Global LR Chord (All cells)
pdf("results/FigS5A_Global_Chord.pdf", width = 12, height = 12)
netVisual_chord_gene(cellchat, lab.cex = 0.6, title.name = "Global Ligand-Receptor Landscape")
dev.off()

# ================= 4. Figure 4G: Pathway Bubble Plot =================
pdf("results/Fig4G_Pathway_Bubble.pdf", width = 14, height = 8)
# Focusing on outgoing signals from tumor subtypes
netVisual_bubble(cellchat, sources.use = tumor_subtypes, targets.use = tme_cells, 
                 remove.isolate = TRUE, title.name = "Tumor-TME Interactions")
dev.off()

# ================= 5. Figure 4H & 4I: PLAU Axis Analysis =================
# Fig 4H: PLAU Role Heatmap
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("results/Fig4H_PLAU_Roles.pdf", width = 8, height = 5)
netAnalysis_signalingRole_network(cellchat, signaling = "PLAU", width = 10, height = 5)
dev.off()

# Fig 4I: PLAU-PLAUR Strength Barplot (Mesen -> Macrophage focus)
prob <- cellchat@net$prob
lr_pair <- "PLAU_PLAUR"
target <- "Macrophages"
df_bar <- data.frame(Source = tumor_subtypes, Prob = 0)

for(i in 1:nrow(df_bar)) {
  src <- df_bar$Source[i]
  # Check if interaction exists to avoid error
  if(src %in% dimnames(prob)[[1]] && target %in% dimnames(prob)[[2]] && lr_pair %in% dimnames(prob)[[3]]) {
    df_bar$Prob[i] <- prob[src, target, lr_pair]
  }
}

pdf("results/Fig4I_PLAU_Barplot.pdf", width = 5, height = 6)
ggplot(df_bar, aes(x = reorder(Source, -Prob), y = Prob, fill = Source)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label=sprintf("%.4f", Prob)), vjust=-0.5) +
  theme_classic() + labs(y = "PLAU-PLAUR Prob (to Macrophages)", x="")
dev.off()

# ================= 6. Figure 4K & S5B: Expression Violins =================
# Fig 4K: Tumor Subtypes Only
pdf("results/Fig4K_PLAU_Tumor_Vln.pdf", width = 6, height = 5)
VlnPlot(subset(obj_filtered, anno_ident %in% tumor_subtypes), 
        features = c("PLAU", "PLAUR"), stack = TRUE, flip = TRUE)
dev.off()

# Fig S5B: All Cells (Evidence for Autocrine + Paracrine)
pdf("results/FigS5B_PLAU_All_Vln.pdf", width = 12, height = 6)
VlnPlot(obj_filtered, features = c("PLAU", "PLAUR"), group.by = "anno_ident", 
        stack = TRUE, flip = TRUE, split.by = "anno_ident") # Split view if needed or standard
dev.off()

# ================= 7. Figure 4L: Survival Analysis =================
# Load CGGA data
cgga <- readRDS("path/to/hgg_data.rds")
expr <- cgga$expression
clin <- cgga$clinical

# Extract PLAUR expression
plaur_idx <- which(expr$Gene_Name == "PLAUR")
plaur_vals <- as.numeric(expr[plaur_idx, -1])
names(plaur_vals) <- colnames(expr)[-1]

# Match clinical
clin <- clin %>% filter(CGGA_ID %in% names(plaur_vals))
clin$PLAUR <- plaur_vals[clin$CGGA_ID]
clin$Group <- ifelse(clin$PLAUR > median(clin$PLAUR, na.rm=TRUE), "High", "Low")

# Survival Plot
fit <- survfit(Surv(OS, Censor..alive.0..dead.1.) ~ Group, data = clin)
pdf("results/Fig4L_PLAUR_Survival.pdf", width = 6, height = 6)
ggsurvplot(fit, data = clin, pval = TRUE, conf.int = TRUE, 
           palette = c("#E7B800", "#2E9FDF"),
           title = "PLAUR Survival Analysis")
dev.off()

message("All figures for Section 1.8 generated.")