#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(qs2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(Matrix)
})

module_dir <- "07_细胞通讯"
outputs_dir <- file.path(module_dir, "outputs")
tables_dir <- file.path(module_dir, "tables")
figures_dir <- file.path(module_dir, "figures")
source_data_dir <- file.path(figures_dir, "source_data")
logs_dir <- file.path(module_dir, "logs")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

obj_path <- file.path(outputs_dir, "CellChat输入_四亚型加纯TME.qs2")
cellchat_comm_path <- file.path(tables_dir, "CellChat纯TME_PLAUR相关通讯_按prob排序.csv")
stopifnot(file.exists(obj_path), file.exists(cellchat_comm_path))

group_order <- c(
  "Malignant_NPC-P", "Malignant_OPC-M", "Malignant_MES-V", "Malignant_MES-I",
  "Macrophages", "Microglial", "Monocytes", "cDCs", "pDCs",
  "T cells", "NK cells", "B cells", "Endothelial", "Mural cells"
)
malignant_groups <- group_order[grepl("^Malignant_", group_order)]
focus_groups <- c("Malignant_MES-V", "Malignant_MES-I", "Macrophages", "Microglial", "Monocytes", "cDCs", "Endothelial", "Mural cells")
focus_genes <- c("PLAU", "PLAUR")

obj <- qs2::qs_read(obj_path)
DefaultAssay(obj) <- "RNA"

genes_present <- intersect(focus_genes, rownames(obj))
if (!setequal(genes_present, focus_genes)) {
  stop("Missing genes in object: ", paste(setdiff(focus_genes, genes_present), collapse = ", "))
}

expr_mat <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(obj, assay = "RNA", layer = "counts")
)
expr <- expr_mat[focus_genes, , drop = FALSE]

md <- obj@meta.data %>%
  as_tibble(rownames = "cell") %>%
  transmute(
    cell,
    Pt_number = as.character(Pt_number),
    comm_group = factor(as.character(comm_group), levels = group_order)
  )

expr_long <- lapply(focus_genes, function(g) {
  x <- as.numeric(expr[g, ])
  tibble(
    cell = colnames(obj),
    gene = g,
    expr = x,
    detected = x > 0
  )
}) %>%
  bind_rows() %>%
  left_join(md, by = "cell") %>%
  filter(!is.na(comm_group))

write_csv(expr_long, file.path(source_data_dir, "PLAU_PLAUR_纯TME_per_cell表达.csv"))

group_summary <- expr_long %>%
  group_by(gene, comm_group) %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(Pt_number),
    pct_detected = mean(detected) * 100,
    mean_expr = mean(expr),
    median_expr = median(expr),
    .groups = "drop"
  ) %>%
  arrange(gene, desc(mean_expr), desc(pct_detected))
write_csv(group_summary, file.path(tables_dir, "PLAU_PLAUR_纯TME_comm_group表达摘要.csv"))

patient_summary <- expr_long %>%
  group_by(gene, comm_group, Pt_number) %>%
  summarise(
    n_cells = n(),
    pct_detected = mean(detected) * 100,
    mean_expr = mean(expr),
    .groups = "drop"
  ) %>%
  arrange(gene, comm_group, desc(n_cells))
write_csv(patient_summary, file.path(tables_dir, "PLAU_PLAUR_纯TME_patient_x_comm_group表达审计.csv"))

patient_dominance <- patient_summary %>%
  mutate(detected_weight = pct_detected * n_cells) %>%
  group_by(gene, comm_group) %>%
  summarise(
    n_patients = n_distinct(Pt_number),
    top_patient = Pt_number[which.max(detected_weight)],
    top_patient_detected_weight_fraction = max(detected_weight, na.rm = TRUE) / sum(detected_weight, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(gene, desc(top_patient_detected_weight_fraction))
write_csv(patient_dominance, file.path(tables_dir, "PLAU_PLAUR_纯TME_patient主导性审计.csv"))

mes_tests <- expr_long %>%
  filter(comm_group %in% c("Malignant_MES-V", "Malignant_MES-I")) %>%
  group_by(gene) %>%
  summarise(
    MES_V_mean = mean(expr[comm_group == "Malignant_MES-V"]),
    MES_I_mean = mean(expr[comm_group == "Malignant_MES-I"]),
    MES_V_pct = mean(detected[comm_group == "Malignant_MES-V"]) * 100,
    MES_I_pct = mean(detected[comm_group == "Malignant_MES-I"]) * 100,
    wilcox_p = wilcox.test(expr ~ comm_group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(wilcox_p, method = "BH"))
write_csv(mes_tests, file.path(tables_dir, "PLAU_PLAUR_MESV_vs_MESI表达检验.csv"))

comm_plaur <- read_csv(cellchat_comm_path, show_col_types = FALSE)
expr_lookup <- group_summary %>%
  select(gene, comm_group, mean_expr, pct_detected) %>%
  mutate(comm_group = as.character(comm_group))

direction_support <- comm_plaur %>%
  filter(ligand == "PLAU", receptor == "PLAUR") %>%
  filter(source %in% focus_groups, target %in% focus_groups) %>%
  left_join(
    expr_lookup %>%
      filter(gene == "PLAU") %>%
      rename(source = comm_group, source_PLAU_mean = mean_expr, source_PLAU_pct = pct_detected) %>%
      select(source, source_PLAU_mean, source_PLAU_pct),
    by = "source"
  ) %>%
  left_join(
    expr_lookup %>%
      filter(gene == "PLAUR") %>%
      rename(target = comm_group, target_PLAUR_mean = mean_expr, target_PLAUR_pct = pct_detected) %>%
      select(target, target_PLAUR_mean, target_PLAUR_pct),
    by = "target"
  ) %>%
  arrange(desc(prob), source, target)
write_csv(direction_support, file.path(tables_dir, "PLAU_PLAUR_纯TME方向表达支撑表.csv"))

plot_df <- expr_long %>%
  filter(comm_group %in% group_order)

p_violin <- ggplot(plot_df, aes(x = comm_group, y = expr, fill = comm_group)) +
  geom_violin(scale = "width", linewidth = 0.15, trim = TRUE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.2, fill = "white", alpha = 0.55) +
  facet_wrap(~gene, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = "log-normalized expression") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(
  file.path(figures_dir, "PLAU_PLAUR_纯TME表达violin.pdf"),
  p_violin,
  width = 8.8,
  height = 6.8,
  useDingbats = FALSE
)

dot_df <- group_summary %>%
  mutate(comm_group = factor(as.character(comm_group), levels = group_order))
write_csv(dot_df, file.path(source_data_dir, "PLAU_PLAUR_纯TME_dotplot_source.csv"))

p_dot <- ggplot(dot_df, aes(x = comm_group, y = gene)) +
  geom_point(aes(size = pct_detected, color = mean_expr)) +
  scale_size_continuous(range = c(1.5, 7), name = "% detected") +
  scale_color_gradient(low = "grey90", high = "#B2182B", name = "Mean expr") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )
ggsave(
  file.path(figures_dir, "PLAU_PLAUR_纯TME表达dotplot.pdf"),
  p_dot,
  width = 7.8,
  height = 3.4,
  useDingbats = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(logs_dir, "05_PLAU_PLAUR表达验证_纯TME_sessionInfo.txt")
)

message("PLAU/PLAUR expression validation outputs written.")
