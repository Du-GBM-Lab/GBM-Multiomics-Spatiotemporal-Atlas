options(encoding = "UTF-8")

suppressPackageStartupMessages({
  library(qs2)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

project_root <- "<DATA_ROOT>/项目/分型/修稿杠生信"
object_path <- file.path(
  project_root,
  "重新分析/05_恶性细胞分亚群与Neftel对照/outputs",
  "GBM.malignant.subtyped.neftel_scored.v2.final_labeled.qs2"
)
figure_dir <- file.path(
  project_root,
  "图片表格/05_发育时间_TF_ATAC验证/02_TF_regulon_SCENIC/正文图候选"
)
source_dir <- file.path(
  project_root,
  "图片表格/05_发育时间_TF_ATAC验证/02_TF_regulon_SCENIC/source_data",
  "FOSL1_PLAUR四亚型表达小提琴图"
)
record_dir <- file.path(
  project_root,
  "图片表格/05_发育时间_TF_ATAC验证/02_TF_regulon_SCENIC/记录"
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(record_dir, recursive = TRUE, showWarnings = FALSE)

reader <- if ("qs_read" %in% getNamespaceExports("qs2")) qs2::qs_read else qs2::qread
obj <- reader(object_path)
stopifnot(inherits(obj, "Seurat"))
DefaultAssay(obj) <- "RNA"

subtype_col <- "subtype_label_final"
genes <- c("FOSL1", "PLAUR")
stopifnot(subtype_col %in% colnames(obj@meta.data))
stopifnot(all(genes %in% rownames(obj)))

subtype_map <- c(
  "Proliferative-NPC" = "NPC-P",
  "OPC-Myelination" = "OPC-M",
  "Vascular-niche MES" = "MES-V",
  "MES-Antigen-presenting" = "MES-I"
)
subtype_order <- c("NPC-P", "OPC-M", "MES-V", "MES-I")
pal <- c(
  "NPC-P" = "#0072B5",
  "OPC-M" = "#E18727",
  "MES-V" = "#20854E",
  "MES-I" = "#BC3C29"
)

dat <- tryCatch(
  FetchData(obj, vars = c(genes, subtype_col), layer = "data"),
  error = function(e) FetchData(obj, vars = c(genes, subtype_col))
)
dat$Subtype <- subtype_map[as.character(dat[[subtype_col]])]
dat <- dat[!is.na(dat$Subtype), c("Subtype", genes)]
dat$Subtype <- factor(dat$Subtype, levels = subtype_order)
dat$cell <- rownames(dat)

long <- bind_rows(lapply(genes, function(g) {
  data.frame(
    cell = dat$cell,
    Subtype = dat$Subtype,
    gene = g,
    expression = as.numeric(dat[[g]])
  )
}))

summary_tbl <- long %>%
  group_by(gene, Subtype) %>%
  summarise(
    n_cells = n(),
    mean_expression = mean(expression),
    median_expression = median(expression),
    pct_detected = mean(expression > 0) * 100,
    mean_positive_expression = ifelse(any(expression > 0), mean(expression[expression > 0]), NA_real_),
    q25 = quantile(expression, 0.25),
    q75 = quantile(expression, 0.75),
    .groups = "drop"
  )

write.csv(long, file.path(source_dir, "FOSL1_PLAUR四亚型表达_source.csv"), row.names = FALSE)
write.csv(summary_tbl, file.path(source_dir, "FOSL1_PLAUR四亚型表达_summary.csv"), row.names = FALSE)

pairwise_stats <- function(d, gene_name) {
  x <- d[d$gene == gene_name, ]
  pairs <- combn(subtype_order, 2, simplify = FALSE)
  out <- lapply(pairs, function(p) {
    a <- x[x$Subtype == p[1], "expression"]
    b <- x[x$Subtype == p[2], "expression"]
    detection_tab <- matrix(
      c(sum(a > 0), sum(a <= 0), sum(b > 0), sum(b <= 0)),
      nrow = 2,
      byrow = TRUE
    )
    data.frame(
      gene = gene_name,
      group1 = p[1],
      group2 = p[2],
      wilcox_p = wilcox.test(a, b)$p.value,
      detection_fisher_p = fisher.test(detection_tab)$p.value,
      pct_group1 = mean(a > 0) * 100,
      pct_group2 = mean(b > 0) * 100,
      detection_delta_group1_minus_group2 = (mean(a > 0) - mean(b > 0)) * 100
    )
  })
  out <- bind_rows(out)
  out$wilcox_p_BH <- p.adjust(out$wilcox_p, method = "BH")
  out$detection_fisher_p_BH <- p.adjust(out$detection_fisher_p, method = "BH")
  out
}

pairwise_tbl <- bind_rows(lapply(genes, function(g) pairwise_stats(long, g)))
write.csv(pairwise_tbl, file.path(source_dir, "FOSL1_PLAUR四亚型表达_pairwise_wilcox.csv"), row.names = FALSE)

theme_pub <- function() {
  theme_classic(base_size = 8.5) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 8.5, color = "black"),
      axis.text.x = element_text(size = 8.2, color = "black", angle = 25, hjust = 1),
      axis.text.y = element_text(size = 8.2, color = "black"),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 7.5, hjust = 0.5, color = "grey35"),
      legend.position = "none",
      plot.margin = margin(7, 8, 5, 8)
    )
}

plot_one <- function(gene_name) {
  d <- long[long$gene == gene_name, ]
  ymax <- max(d$expression, na.rm = TRUE)
  ggplot(d, aes(Subtype, expression, fill = Subtype)) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      color = "#303030",
      linewidth = 0.28,
      alpha = 0.88
    ) +
    geom_boxplot(
      width = 0.14,
      outlier.shape = NA,
      fill = "white",
      color = "#202020",
      linewidth = 0.28,
      alpha = 0.85
    ) +
    scale_fill_manual(values = pal) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title = gene_name,
      subtitle = "Malignant scRNA four-subtype expression",
      y = "log-normalized expression"
    ) +
    theme_pub()
}

for (g in genes) {
  p <- plot_one(g)
  ggsave(
    filename = file.path(figure_dir, paste0(g, "四亚型表达小提琴图.pdf")),
    plot = p,
    width = 2.55,
    height = 2.35,
    device = cairo_pdf
  )
}

combined <- ggplot(long, aes(Subtype, expression, fill = Subtype)) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    color = "#303030",
    linewidth = 0.25,
    alpha = 0.88
  ) +
  geom_boxplot(
    width = 0.13,
    outlier.shape = NA,
    fill = "white",
    color = "#202020",
    linewidth = 0.25,
    alpha = 0.85
  ) +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(y = "log-normalized expression", x = NULL) +
  theme_pub() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    panel.spacing.x = unit(0.35, "cm")
  )
ggsave(
  filename = file.path(figure_dir, "FOSL1_PLAUR四亚型表达小提琴图.pdf"),
  plot = combined,
  width = 5.1,
  height = 2.35,
  device = cairo_pdf
)

fosl1_metric <- summary_tbl %>%
  filter(gene == "FOSL1") %>%
  mutate(
    Subtype = factor(Subtype, levels = subtype_order),
    positive_label = sprintf("%.1f%%", pct_detected),
    mean_pos_label = sprintf("%.2f", mean_positive_expression)
  )

fosl1_metric_long <- bind_rows(
  data.frame(
    Subtype = fosl1_metric$Subtype,
    metric = "Detected cells (%)",
    value = fosl1_metric$pct_detected,
    label = fosl1_metric$positive_label
  ),
  data.frame(
    Subtype = fosl1_metric$Subtype,
    metric = "Mean expression in detected cells",
    value = fosl1_metric$mean_positive_expression,
    label = fosl1_metric$mean_pos_label
  )
)

fosl1_metric_plot <- ggplot(fosl1_metric_long, aes(Subtype, value, fill = Subtype)) +
  geom_col(width = 0.66, color = "#303030", linewidth = 0.25, alpha = 0.92) +
  geom_text(aes(label = label), vjust = -0.35, size = 2.4, color = "#202020") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    title = "FOSL1",
    subtitle = "Zero-inflated expression summarized by detection rate and positive-cell expression",
    y = NULL,
    x = NULL
  ) +
  theme_pub() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 8.5, face = "bold"),
    panel.spacing.x = unit(0.38, "cm")
  )

ggsave(
  filename = file.path(figure_dir, "FOSL1四亚型表达_检测率阳性表达.pdf"),
  plot = fosl1_metric_plot,
  width = 5.15,
  height = 2.35,
  device = cairo_pdf
)

ggsave(
  filename = file.path(figure_dir, "FOSL1_PLAUR四亚型表达小提琴图.pdf"),
  plot = combined,
  width = 5.1,
  height = 2.35,
  device = cairo_pdf
)

record <- c(
  "FOSL1 / PLAUR 四亚型表达小提琴图",
  "",
  paste0("对象: ", object_path),
  paste0("subtype 字段: ", subtype_col),
  "表达层: RNA assay, log-normalized data",
  "",
  "输出图:",
  paste0("- ", file.path(figure_dir, "FOSL1四亚型表达小提琴图.pdf")),
  paste0("- ", file.path(figure_dir, "FOSL1四亚型表达_检测率阳性表达.pdf")),
  paste0("- ", file.path(figure_dir, "PLAUR四亚型表达小提琴图.pdf")),
  paste0("- ", file.path(figure_dir, "FOSL1_PLAUR四亚型表达小提琴图.pdf")),
  "",
  "source data:",
  paste0("- ", file.path(source_dir, "FOSL1_PLAUR四亚型表达_source.csv")),
  paste0("- ", file.path(source_dir, "FOSL1_PLAUR四亚型表达_summary.csv")),
  paste0("- ", file.path(source_dir, "FOSL1_PLAUR四亚型表达_pairwise_wilcox.csv")),
  "",
  "summary:",
  paste(capture.output(print(summary_tbl)), collapse = "\n"),
  "",
  "边界:",
  "- 该图是 malignant scRNA expression-level evidence。",
  "- 可用于说明 FOSL1 / PLAUR 在四个恶性亚型中的表达格局。",
  "- 不能单独证明 FOSL1 因果调控 PLAUR；因果由 ChIP-qPCR / luciferase / perturbation 承重。"
)
writeLines(record, file.path(record_dir, "FOSL1_PLAUR四亚型表达小提琴图记录.txt"))
cat(paste(record, collapse = "\n"), "\n")
