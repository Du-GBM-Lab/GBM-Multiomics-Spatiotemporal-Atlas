# 04_inferCNV_免疫参考验证

本目录用于 Genes & Diseases major revision 中 malignant cell identification 的严格 reference 版本。

定位：

1. 使用最严格 lymphoid/pDC-only reference 重新运行 sample-wise inferCNV。
2. 基于 `cnv_burden_z` 和 `cnv_correlation_ref_z` 两轴定义 high-confidence malignant cells。
3. 与 `03_inferCNV_恶性识别` 中 broad-reference inferCNV 结果做 sensitivity / concordance analysis。
4. 输出 reviewer response 需要的 marker validation 图：
   - 原始 cell type annotation 的谱系 / 恶性 / reference marker bubble plot。
   - malignant high-confidence 与 non-malignant/CNV-low 之间的 marker expression validation plot。

脚本顺序：

```text
01_samplewise_infercnv_immune_reference.R
02_summarize_infercnv_immune_reference_calls.R
03_marker_bubble_plots_immune_reference.R
```

当前只写代码，尚未运行正式分析。
