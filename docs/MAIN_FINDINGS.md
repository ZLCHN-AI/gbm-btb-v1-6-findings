# GBM–BTB Spatial Remodeling | Main Findings (Manuscript V1.6)

**Status:** Working manuscript, V1.6 · **Updated:** 2026-11-01  
**Study type:** Secondary computational analysis of public, de-identified human glioblastoma transcriptomic resources.  
**Central question:** How do vascular and blood–tumor barrier (BTB)-associated transcriptional features vary between histopathologically annotated glioblastoma regions, and which observations remain interpretable after accounting for tissue composition?

> **Scope:** These are **transcriptional associations**, not direct measurements of vascular permeability, endothelial transcytosis, barrier opening, reversible target modulation, or drug delivery.

## 1. Opposing anatomical patterns of two barrier-associated transcriptional modules

In the GBM-Space discovery resource, the infiltrative-tumor (**IT**) to cellular-tumor (**CT**) comparison retained **17 eligible Visium sections from 4 donors**. The two separately computed tissue-level programs moved in opposite directions:

| Program | Gene set used in V1.6 | Adjusted donor-median CT − IT effect | Donor direction |
|---|---|---:|---|
| Junction-associated core | `CLDN5`, `OCLN`, `TJP1` | **−0.25278** | **4/4 lower in CT** |
| Transcytosis-associated core | `PLVAP`, `CAV1`, `CAV2`, `CAVIN1` | **+0.21665** | **4/4 higher in CT** |
| Previously reported transcytosis-minus-junctional contrast | A derived composite, not an additional independent biological program | **+0.21296** | **4/4 higher in CT** |

These data reveal **regionally patterned tissue-level vascular/BTB-associated transcription**. Four concordant donor directions correspond to a **two-sided exact sign-test P = 0.125**; they do not establish conventional statistical significance. The two module effects are not interchangeable with the separately defined composite contrast.

## 2. PLVAP retains a modest regional association after myeloid-abundance adjustment

Gene-level IT-to-CT effects were evaluated under a base tissue-composition model accounting for endothelial, pericyte, and malignant abundances, then an expanded model additionally accounting for **myeloid abundance**.

| Gene | Base-adjusted CT − IT | Expanded (+ myeloid) CT − IT | Donor direction after expanded adjustment |
|---|---:|---:|---|
| **PLVAP** | **+0.161** | **+0.135** | **4/4 positive** |
| CAV1 | +0.343 | +0.185 | 3/4 positive; one approximately zero |
| CAVIN1 | +0.262 | +0.181 | 4/4 positive |

**PLVAP** preserved a small, consistently positive regional signal after the tested adjustments. **CAV1** was more sensitive to composition adjustment. A sensitivity analysis excluding CAV1 attenuated the adjusted transcytosis-associated module effect, although the reported donor-level direction remained positive in 4/4 donors.

**Interpretation:** PLVAP is a comparatively robust **tissue-level correlate**, not a proven endothelial-specific spatial source or a causally validated permeability regulator. The spatial effect is small, with one donor close to zero after adjustment.

## 3. A donor-consistent PLVAP-high / MFSD2A-low endothelial transcriptional state

A separately analyzed **single-nucleus layer within the same GBM-Space atlas** included **2,241 endothelial nuclei from 12 donors** (1,363 canonical and 878 remodeling-associated endothelial nuclei). In annotation-defined remodeling-associated versus canonical endothelium:

| Readout | Donor-level effect | Direction across donors |
|---|---:|---|
| PLVAP detection | **+0.26409** | **12/12 higher** |
| PLVAP pseudobulk expression | **+1.61784 log1pCPM** | **11/12 higher** |
| MFSD2A pseudobulk expression | **−0.71448 log1pCPM** | **10/12 lower** |

These results provide **cell-resolved endothelial-state evidence**. The nuclei and Visium spots originate from the **same parent atlas**, not independent validation cohorts. They were **not linked cell-by-cell or spot-by-spot** to specific IT/CT locations. In particular, the lower single-nucleus MFSD2A expression **was not recapitulated as a lower composition-adjusted spatial effect**.

## 4. Complementary external and localization evidence

- **GSE162631:** Patient-paired **core-versus-periphery endothelial** analyses support related regional endothelial remodeling patterns; core/periphery is **not equivalent** to GBM-Space's strict IT/CT pathology classification.
- **Ivy GAP / GSE107559 (ASTR122):** Histopathology-resolved **bulk RNA-seq** provides tumor-level evidence for anatomical vascular transcription differences; bulk-region expression cannot establish cell-of-origin.
- **Sonpatki et al. Xenium S2/S6:** Restricted in situ transcript-localization context indicates **limited CAV1 endothelial specificity** and heterogeneous PLVAP enrichment between samples. Incomplete labels, panel coverage, and cell-segmentation uncertainties constrain interpretation.

These resources offer **complementary evidence, not four interchangeable replications of one spatial statistic**.

## Biological significance

Together, V1.6 identifies **opposing IT-to-CT tissue-level junction- and transcytosis-associated transcriptional patterns**, a **modest PLVAP regional signal that persists after the tested composition adjustments**, and a **reproducible remodeling-associated endothelial PLVAP-high/MFSD2A-low transcriptional profile**. The integrative contribution is to separate **anatomical tissue-level associations** from **cell-resolved endothelial-state evidence**, while evaluating marker specificity and composition sensitivity.

### Evidence boundaries

1. The principal IT-to-CT spatial comparison involves **4 independent donors**; **P = 0.125** for 4/4 direction concordance by two-sided exact sign test.
2. Visium measures **mixed-cell tissue expression**. Higher transcytosis-associated RNA module scores **do not demonstrate increased endothelial vesicular transport**.
3. PLVAP has a **small adjusted effect**; CAV1 is **composition-sensitive** and is not an endothelial-specific tissue marker.
4. The endothelial snRNA results and tissue spatial results are **compatible but not spatially linked**. MFSD2A findings differ across these analytical levels.
5. This study does **not** establish a causal target, reversible BTB permeability modulation, or improved therapeutic delivery; these require dedicated functional experiments.

## Primary data resources

- [GBM-Space multimodal atlas](https://www.gbmspace.org/) — primary discovery resource (Visium and snRNA layers).
- [GSE162631](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162631) — independent patient-paired endothelial comparisons.
- [Ivy GAP / GSE107559](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107559) — histopathology-resolved bulk regional transcriptomics.
- [Xenium S2/S6 source records](https://doi.org/10.5281/zenodo.17622242) — limited cross-technology localization.

**Reproducibility note:** The repository's source-code and analysis artifacts should be interpreted according to their existing provenance and availability documentation. The present Markdown summary does not assert that the entire analysis is runnable end-to-end from raw data in one command. Final figures and supplementary tables are being reconciled against manuscript V1.6.

---

## 中文摘要｜V1.6 主要科学发现

**研究定位：** 基于公开 GBM-Space 图谱的原创二次计算分析，以病理解剖区域的血管/BTB 相关转录异质性为核心，不将转录评分等同于真实通透性。

**发现一｜IT→CT 的两类屏障相关转录模块呈相反变化。** 在17张符合筛选条件的空间切片、4名独立供体中，连接相关模块（CLDN5/OCLN/TJP1）的校正后供体中位效应为 **−0.25278（4/4下降）**；转胞运输相关模块（PLVAP/CAV1/CAV2/CAVIN1）为 **+0.21665（4/4上升）**。方向一致，但4名供体双侧精确符号检验 **P=0.125**，属于探索性证据。

**发现二｜PLVAP 的小幅空间上升对所检验的细胞组成校正相对稳健。** IT→CT 效应由 **+0.161** 变为额外校正髓系丰度后的 **+0.135（4/4同向）**。相较之下，CAV1 由 **+0.343** 降至 **+0.185（3/4同向）**。这些是组织层面的相关性，不能证明 PLVAP 在所有空间点位均来源于内皮。

**发现三｜存在 PLVAP 较高、MFSD2A 较低的重塑相关内皮细胞状态。** GBM-Space 单核层12供体共2,241个内皮细胞核，PLVAP 检出率差异 **+0.26409（12/12同向）**，PLVAP pseudobulk **+1.61784（11/12）**，MFSD2A pseudobulk **−0.71448（10/12下降）**。单核层的 MFSD2A 下降不能直接推广为 IT→CT 的空间校正后下降。

**研究价值：** 在同一个多模态发现平台上，用空间区域比较、供体级统计、细胞组成敏感性和单核内皮状态分析，区分组织整体的血管相关转录现象与内皮细胞自身的表达状态，为后续寻找可干预的 BTB 调控因素提供依据；**目前没有直接证明可逆性屏障开放或药物递送改善。**