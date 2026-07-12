# GBM Multi-omics Spatiotemporal Atlas — analysis code

Code accompanying the manuscript. Reproduces the malignant four-subtype taxonomy (NPC-P / OPC-M / MES-V / MES-I) and the FOSL1–PLAU–PLAUR analyses.

## Data source

- Raw single-cell / spatial data: GSA-Human **HRA004677**, GSA **CRA011176** (Mei et al.).
- Processed integrated object (four-subtype annotations, metadata, reductions): Zenodo **10.5281/zenodo.21311384**.

## Environment

Install dependencies first: `00_install/00_install_required_packages.R`. R 4.5.x; Python for SCENIC steps. *(TODO before upload: pin package versions / `sessionInfo`.)*

## Pipeline (run in order)

| 阶段 | 目录 | 源（`修稿杠生信\重新分析\`） |
|---|---|---|
| Preprocessing / QC / malignant calling | `01_preprocessing_QC/` | `02_scRNA_QC` + `03_inferCNV_恶性识别` + `04_inferCNV_免疫参考验证` |
| Malignant subtyping (NMF + Neftel) | `02_malignant_subtyping_NMF/` | `05_恶性细胞分亚群与Neftel对照` |
| SCENIC gene regulatory network | `03_SCENIC_GRN/` | `08_发育时间_TF_ATAC验证\02_TF_regulon_SCENIC` |
| RCTD spatial deconvolution | `04_RCTD_spatial_deconvolution/` | `R9_空间转录组\scripts` |
| Trajectory / pseudotime | `05_trajectory/` | `06_恶性细胞拟时序` |
| Cell–cell communication (CellChat) | `06_cellchat/` | `07_细胞通讯` |
| Clinical / survival | `07_clinical_survival/` | `R4_临床生存分析` |
| Co-expression (hdWGCNA) | `08_scWGCNA/` | `R2_scWGCNA` |
| Cross-species projection | `09_cross_species/` | `08_发育时间_TF_ATAC验证\01_mouse_developmental` |

## 收录范围

- 仅**分析脚本**（`.R/.py/.md/.sh/.ipynb`）；已排除内嵌库（`packages/`）、数据产物（`.qs2/.rds/.csv`）、图件、日志。
- **不含**已弃旧版：CARD 去卷积、scTenifoldKnk 虚拟敲除、旧 6 样本 ATAC（不在当前 code availability 内）；也不含内部审计小工具。

## 运行与路径

- 运行顺序、数据放置见 **`00_RUN_ORDER.md`**。
- 脚本内本机绝对路径已统一去本地化为占位符 **`<DATA_ROOT>/`**（R 库路径为 `<R_LIBS>`）；运行前按 `00_RUN_ORDER.md` 设置。
- `.gitignore` 已排除数据/产物/库。

## 上传前 TODO（剩余）

- 版本锁定（`renv.lock` 或每阶段 `sessionInfo`）。
- 选开源许可证并加 `LICENSE`。
- 各阶段目录可补一行入口说明。
