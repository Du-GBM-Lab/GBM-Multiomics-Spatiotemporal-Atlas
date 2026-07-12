# 恶性细胞拟时序章节 figure 规划与 response letter 骨架

更新时间：2026-05-21

## 1. 总体定位

Trajectory analysis 不作为“确定性发育分化”的单独证据，而是作为 subtype framework 之后的 secondary, hypothesis-supporting analysis。

本章回答四个问题：

1. 四个 malignant subtype 是孤立 cluster，还是位于同一 transcriptional landscape 上。
2. NPC-P 是否有足够证据作为最 plastic / root-like state。
3. MES-V 与 MES-I 是不是可以在 shared MES-like backbone 中被 endpoint program 区分。
4. 轨迹分析能否提供 subtype analysis 之外的 switch / driver gene 信息。

安全表述：

- inferred trajectory
- lineage-compatible ordering
- transcriptional continuum
- shared MES-like transcriptional backbone
- endpoint divergence
- niche-associated program engagement

避免表述：

- deterministic lineage commitment
- real-time progression
- MES-I differentiates into MES-V
- MES-V transforms into MES-I

## 2. Main Figure 规划

核心原则：main figure 只放跟命名一致、能直接回答 R1/R2 的 panel。任何可能造成误读或自相矛盾的 panel 降 supplementary 或删除。

| Panel | 内容 | 来源 | 主要回应 |
|---|---|---|---|
| A | CytoTRACE2 by subtype，NPC-P 最高；旁边一行 multi-metric root summary | Phase 0 | R2 root |
| B | malignant-only UMAP + Slingshot 3-lineage curves，NPC-P root，star pattern | Phase 1 | R1 core / R2 |
| C | Lineage-specific pseudotime density，L1/L2/L3 分开 | Phase 1 | R1 continuum |
| D | Cross-method + multi-root robustness 缩略图：NPC-P 在不同 root 下恒为中心 hub | Phase 2A/2B | R2 stability |
| E | MES-I vs MES-V endpoint comparison：MP02/MP04 + RGS5/ACTA2/TAGLN + CD74/HLA-DRA/HLA-DPB1，附 Fig 07 validation caveat | Phase 3 Panel D + Phase 3 validation | R1 contribution / MES separability |
| F | MES-I vs MES-V diffEnd volcano：switch genes，marker 方向标注 | Phase 4 | R1 switch gene contribution |

Main figure 的逻辑重心：

- A + D：回答 R2 的 root selection 和 stability。
- E + F：回答 R1 的 trajectory contribution。
- B + C：把 trajectory topology 和 subtype continuum 可视化。

## 3. Supplementary Figure 规划

| Supp | 内容 | 来源 | 主要回应 |
|---|---|---|---|
| S1 | Phase 0 multi-metric root justification：cell cycle + SOX family + GAP43 + Neftel 6-state，带统计 | Phase 0 | R2 GAP43 / SOX family |
| S2 | Monocle3 cross-validation：principal graph / pseudotime / scatter / Sankey + branch-len diagnostic table | Phase 2A / 03b | R2 method robustness |
| S3 | Root sensitivity 4-panel + correlation heatmap。Caption 必须说明 topology invariant，pseudotime is root-relative | Phase 2B | R2 root stability |
| S4 | CytoTRACE2 vs pseudotime，标题用诚实版，不写 global decrease | Phase 2B Fig 05 | R2 stemness |
| S5 | Marker validation：CNV 100% / doublet 0% / 共表达 / MHC-II vs TAM / vascular vs mural | Phase 3 validation | MES-I/V marker genuineness |
| S6 | MP01-06 + Neftel 6-state along pseudotime，L1/L2 dashed 或 caveat，L3 solid | Phase 3 A/B | Biology overlay detail |
| S7 | L3 driver heatmap + GO enrichment。Caption 说明 GO background = 600-gene tradeSeq universe，GO themes 是 focused screen 的 pathway theme，不定义 subtype identity | Phase 4 | switch gene / pathway theme |

明确删除或不进 figure：

- L1/L2 association saturated bar chart：零信息量，只保留 supplementary source table。
- 未加 caveat 的 cross-run pseudotime correlation heatmap：容易被误读为 trajectory 不稳。

## 4. GO enrichment background lock

当前 Phase 4 enrichment 已修正为：

```text
GO background = 600-gene tradeSeq gene universe
Mapped Entrez universe = 559 genes
```

对应文件：

```text
tables/driver_enrichment_background.csv
tables/driver_enrichment_per_cluster.csv
figures/08_driver_enrichment_dotplot.pdf
```

写作边界：

- 可以引用 GO as pathway theme within focused tradeSeq screen。
- 不要把 GO enrichment 写成 genome-wide unbiased pathway discovery。
- 不要用 GO crossover 反向定义 subtype identity；subtype identity 仍由 Phase 3 endpoint marker / metaprogram divergence 支撑。

S7 caption 可用：

```text
GO enrichment was performed using the 600-gene tradeSeq testing universe as background. Enriched terms are therefore interpreted as pathway themes among trajectory-associated candidate genes, rather than as genome-wide subtype-defining programs.
```

## 5. Reviewer 1 response skeleton — motivation and contribution

Reviewer 1 concern：trajectory / pseudotime analysis motivation 不清，figure 不够 conclusive。

Draft response：

```text
We thank the reviewer for this comment and have substantially reframed and extended the trajectory analysis. We now state its purpose explicitly: rather than serving as standalone evidence for deterministic lineage transitions, the trajectory analysis is a secondary, hypothesis-supporting analysis that tests whether the four malignant subtypes form a reproducible transcriptional continuum and characterizes how proliferation, stemness, and metaprogram engagement vary across that continuum.

Specifically, the revised analysis provides four contributions beyond the static subtype and NMF analyses:

1. It demonstrates that the four subtypes are organized as a connected transcriptional landscape with Proliferative-NPC-like (NPC-P) as the central, most stem-like state, from which three lineage-compatible branches diverge toward OPC-M-like, MES-V-like, and MES-I-like states.

2. It shows that MES-V-like and MES-I-like cells, although transcriptionally adjacent, can be resolved as distinct MES-like termini, supporting their separation as biologically meaningful niche-associated malignant states.

3. It identifies the molecular programs distinguishing the two MES-like termini: a vascular / mural-mimicry program (MP02, RGS5, ACTA2, TAGLN) in MES-V-like cells and an immune-reactive / antigen-presentation-associated program (MP04, HLA-DRA, HLA-DPB1, CD74) in MES-I-like cells.

4. It nominates trajectory-associated switch / driver genes and pathway themes, including endpoint genes distinguishing MES-I-like and MES-V-like termini.

These results are summarized in revised Fig. X and Supplementary Figs. S1-S7.
```

可填数字：

```text
MP02: MES-V terminal mean = 0.552, MES-I terminal mean = 0.406, BH p = 5.62e-134.
MP04: MES-I terminal mean = 0.428, MES-V terminal mean = 0.288, BH p = 1.80e-159.
Strict endpoints: MES-I n = 694, MES-V n = 1,135.
```

注意：

- R1 回复里不要把 GO enrichment 当 subtype contribution 的主证据。
- Contribution 落在 endpoint marker divergence + switch gene nomination。

## 6. Reviewer 2 response skeleton — root justification and stability

Reviewer 2 concern：GAP43 不足以定 root；需要 root 生物学理由、SOX family 对比、stemness hierarchy stability。

### R2.1 GAP43

```text
We agree that GAP43, a neuronal growth-associated gene, is insufficient to define a developmental origin. In our re-analysis, although GAP43 mean expression is nominally highest in NPC-P-like cells, its median expression is close to zero across all four subtypes, confirming that it is too sparsely expressed to serve as a single-marker root criterion. We therefore replaced GAP43-based root selection with a multi-metric framework.
```

### R2.2 Root rationale + SOX family

```text
Root selection is now based on the convergence of four independent indicators, all of which identify NPC-P-like as the most stem-like / least differentiated origin candidate:

(i) CytoTRACE2 stemness score, which is highest in NPC-P-like cells;
(ii) cell cycle phase composition, with NPC-P-like cells showing the highest S/G2M fraction and OPC-M-like cells being predominantly G1;
(iii) SOX family expression across all subtypes (SOX2/SOX4/SOX9/SOX11), with SOX2 and SOX11 enriched in NPC-P-like cells;
(iv) Neftel NPC1/NPC2 AMS scores, which are highest in NPC-P-like cells.
```

可填数字：

```text
Patient n = 24.
High-confidence malignant cells = 28,213.
OPC-M G1 fraction = 98.1%.
```

### R2.3 Stability

```text
We assessed robustness of the inferred hierarchy along four axes.

First, cross-method validation using Monocle3 reproduced the NPC-P -> MES-V (L2) and NPC-P -> OPC-M (L3) trajectories with high concordance (Spearman rho = 0.95 for both). The NPC-P -> MES-I (L1) trajectory showed lower fine-scale concordance (rho = 0.11); a parameter sensitivity analysis (minimal_branch_len = 10/15/20) confirmed that this reflects the transcriptional continuity between MES-V-like and MES-I-like states rather than a simple parameter artifact.

Second, root sensitivity analysis across four root specifications (NPC-P, OPC-M, MES-V, and algorithm-determined) showed that the topology was qualitatively stable: NPC-P-like consistently occupied the central branching position through which lineages passed. Absolute pseudotime values differed between runs, as expected because pseudotime is measured relative to the chosen root.

Third, stemness coupling supported the hierarchy most strongly along the OPC-M lineage, where CytoTRACE2 declined sharply with pseudotime. MES-like lineages showed weaker and non-monotonic stemness coupling, consistent with MES-like states arising through stress / niche-associated state remodeling rather than classical developmental differentiation.
```

可填数字：

```text
Slingshot vs Monocle3:
- L1 rho = 0.112
- L2 rho = 0.948
- L3 rho = 0.949

Branch length diagnostic:
- L3 rho remains 0.83-0.95 across minimal_branch_len settings.

kNN:
- MES-V <-> MES-I obs/exp = 0.38, highest inter-subtype pair.

CytoTRACE2 vs pseudotime:
- L1 rho = -0.205
- L2 rho = -0.095
- L3 rho = -0.749
```

## 7. Limitations paragraph

```text
We note three limitations. First, pseudotime represents an inferred transcriptional ordering and does not imply real-time progression or deterministic lineage commitment. Second, fine-scale ordering within the MES-like compartment is method-sensitive owing to the transcriptional continuity of MES-V-like and MES-I-like states; we therefore base MES-related conclusions on endpoint comparisons rather than continuous ordering. Third, ambient-RNA correction tools were not applied; however, the MES-I-like immune program and MES-V-like mural program are supported by CNV-confirmed malignancy, single-cell co-expression, absence of doublet enrichment, and expression levels below bona fide TAM/mural references, supporting cell-intrinsic rather than wholesale contaminating signal.
```

可填验证数字：

```text
MES-I terminal CNV+ = 100%.
MES-V terminal CNV+ = 100%.
DoubletFinder doublet_pct = 0% for both strict termini.
MES-I terminal: CNV+ and any MHC-II positive = 93.4%.
MES-V terminal: CNV+ and any vascular marker positive = 93.6%.
No SoupX / CellBender / DecontX output was found, so do not claim ambient-correction-proven expression.
```

## 8. Final writing locks

```text
Main claim:
MES-V and MES-I share a MES-like transcriptional backbone but diverge into vascular-like and immune-reactive endpoint programs.

Do not claim:
- deterministic MES-V/MES-I transition;
- functional antigen presentation by MES-I cells;
- genome-wide unbiased GO discovery from Phase 4;
- CD74 or RGS5 as top tradeSeq diffEnd genes.

Can claim:
- Phase 3 strict endpoint comparison strongly supports canonical marker divergence;
- Phase 3 validation supports CNV-confirmed malignant, non-doublet endpoint cells;
- Phase 4 focused tradeSeq screen nominates additional endpoint switch genes and pathway themes.
```
