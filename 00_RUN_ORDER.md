# Run order & data setup

## Data setup

Scripts reference input/output data by a placeholder root **`<DATA_ROOT>/`** (the local machine-specific absolute path was removed for release). Before running, do one of:

- Set an environment variable and adapt the scripts to read it, **or**
- Find-and-replace `<DATA_ROOT>/` with the base directory that holds the analysis data on your machine.

Raw data: GSA-Human **HRA004677** / GSA **CRA011176** (Mei et al.). Processed integrated object: Zenodo **10.5281/zenodo.21311384**. Some `<DATA_ROOT>/环境/…` references are R-library locations — point them at your own library or delete the `.libPaths(...)` lines.

Install dependencies first: `00_install/00_install_required_packages.R` (R 4.5.x; Python 3 for SCENIC — set the `PYTHON` env var if auto-detection fails).

## Pipeline order

| # | 目录 | 作用 |
|---|---|---|
| 1 | `01_preprocessing_QC/` | QC（<20% MT + 3-MAD + DoubletFinder）、inferCNV 恶性判定 |
| 2 | `02_malignant_subtyping_NMF/` | de novo NMF + Neftel 对照 → NPC-P/OPC-M/MES-V/MES-I |
| 3 | `03_SCENIC_GRN/` | SCENIC regulon / 基因调控网络 |
| 4 | `04_RCTD_spatial_deconvolution/` | RCTD 空间去卷积 + niche 统计 |
| 5 | `05_trajectory/` | 拟时序 / 定根稳健性 |
| 6 | `06_cellchat/` | 细胞通讯（PLAU–PLAUR 提名） |
| 7 | `07_clinical_survival/` | bulk 队列生存关联 |
| 8 | `08_scWGCNA/` | hdWGCNA 共表达 |
| 9 | `09_cross_species/` | 跨物种投射验证 |

1–3 有先后依赖（对象逐级产出）；4–9 多为并行的下游分析，各自读上游产出对象。

## 仍待补（上传前）

- 版本锁定：`renv.lock` 或每阶段 `sessionInfo()`（当前仅有 `00_install` 依赖清单）。
- 各阶段目录内可加一行 README 说明入口脚本与产出。
- 选一个开源许可证（如 MIT）并加 `LICENSE`。
