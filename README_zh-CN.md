# 胶质母细胞瘤血脑肿瘤屏障相关转录重塑：V1.6

**独立科研项目展示仓库（工作稿，2026-11-01）**  
[English README](README.md) · [完整科研发现](docs/MAIN_FINDINGS.md)

## 科学问题

GBM 病理学浸润区（IT）与细胞性肿瘤区（CT）是否存在不同的血管/BTB 相关转录程序？在考虑细胞组成后，哪些区域性分子关联依然保持？

## 三项主要发现

1. **IT→CT 两类组织转录模块呈相反方向变化**：连接相关模块 **−0.25278**、转胞运输相关模块 **+0.21665**；17张符合条件的 Visium 切片来自4名供体，两项结果各有 **4/4** 供体同向。
2. **PLVAP 具有较小但方向一致的校正后区域性关联**：基础组成校正后 **+0.161**，进一步校正髓系丰度后 **+0.135（4/4供体同向）**。相较之下，CAV1 的效应由 **+0.343** 降至 **+0.185（3/4）**。
3. **PLVAP-high/MFSD2A-low 重塑相关内皮转录状态**：GBM-Space 单核层12供体、2,241个内皮细胞核；PLVAP 检出率12/12供体升高、pseudobulk 11/12升高、MFSD2A pseudobulk 10/12降低。

**证据边界：** 4/4供体同向的双侧精确符号检验为 **P=0.125**。Visium 为多细胞混合组织点位，不是内皮细胞特异的通透性测量；内皮单核状态与空间点位没有一一对应，也没有证实药物递送或屏障可逆调控。

## 本仓库提供什么？

- `docs/MAIN_FINDINGS.md`：中英文完整发现和限制。
- `results/`：三张由已报告数值整理的**汇总数据表**，**不是原始测序数据或供体级完整输出**。
- `figures/`：Figure 2B 汇总图初稿，并非供体逐点主图。
- `docs/DATA_SOURCES.md`：GBM-Space 与外部数据资源来源及用途。
- `docs/RELEASE_NOTES.md`：后续需要完善的事项。
- `figure_code/`：可重新运行的出图代码。`make_figures.py` 用 `figure_code/data/` 中冻结的供体级和汇总数据表重画稿件的 Fig. 2、Fig. 4 和补充图 S2–S3（`python make_figures.py`，需要 numpy、pandas、matplotlib）。
- `ivygap_ish.py`、`cache_ivygap/`、`results_ivygap/`：Ivy GAP 原位杂交数据的获取与分析脚本、可离线复现的 API 缓存和结果。

原有仓库 `ZLCHN-AI/gbm-btb-spatial-remodeling` 保持不变。本项目作为**独立新仓库**公开 V1.6 工作稿的科研发现；上游分析脚本见该原有仓库。本仓库不包含原始 RNA-seq、Visium 或单核测序数据，不从原始数据重跑上游分析，也不涉及功能性实验；第三方数据须从原始提供方获取。