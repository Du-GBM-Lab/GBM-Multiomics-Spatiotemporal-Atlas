# R4 临床生存分析

## 定位

R4 是意义层 / rationale：回答 MES-like / mesenchymal program 高表达是否与较差 OS 相关，从而解释后续为什么值得追 MES-V / PLAU-PLAUR / FOSL1-PLAUR 机制线。

## 硬边界

- R4 不提 PLAUR。PLAUR 在 R5 之后才进入叙事。
- R4 不包含 PLAUR 单基因生存分析。
- R4 不包含 FOSL1 / FOSL2。
- R4 只写 associated / prognostic / stratified / consistent with，不写 causes / drives / causal。
- 主终点是 OS；PFS 仅作有数据时的 sensitivity。
- 主效应量是 continuous per-SD Cox；KM 只是可视化。
- optimal cutpoint 只作 sensitivity，不作主结果。

## 当前 v2 设计

1. 主分析回到原稿思路：四个新亚型 gene sets（NPC-P / OPC-M / MES-V / MES-I）投射到 CGGA HGG，ssGSEA 后无监督聚类，形成 bulk clusters 并做四臂 KM/log-rank。
2. R4 主线在 cluster survival + cluster clinical annotation 收口；CellChat / PLAU-PLAUR / PLAUR 单基因生存不进 R4。
3. 旧 C1-C4 注释不能照搬；新 cluster 按实际四亚型 score fingerprint 与 IDH / 1p19q / MGMT / grade 分布重新注释。
4. IDH-wt GBM 内 MES-V continuous per-SD Cox 作为 backstop，用于回应“是否只是 IDH/grade 梯度”的审稿风险。
5. 打分方式：ssGSEA 使用 GSVA 新 API `gsva(ssgseaParam(expr, geneSets))`；z-score mean 作稳健性对照。

## STOP 点

- STOP 0：包版本、signature 提取、禁用基因检查。
- STOP 1：CGGA 文件结构、样本口径、grade/IDH/MGMT/OS 字段、表达量纲与 ID 对齐审计。
- STOP 2：signature score 分布与缺失值审计。
- STOP 3：Cox / KM 结果判读。
