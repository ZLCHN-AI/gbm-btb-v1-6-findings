#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ivygap_ish.py — Ivy GAP 原位杂交(ISH)按解剖结构的表达量提取与配对比较

用途
----
从 Allen Institute 的 Ivy Glioblastoma Atlas Project (Ivy GAP) 公共 API 提取
指定基因的 ISH 定量结果，按 GBM 解剖结构（Leading Edge / Infiltrating Tumor /
Cellular Tumor / Microvascular Proliferation 等）汇总，并在"同一张切片内"做
配对比较，用于为转录组层面的结论提供独立模态的佐证。

为什么要配对比较
----------------
Allen 的 expression_energy / expression_density 是按结构面积归一化的，而不同
解剖结构的细胞密度差异极大，直接跨结构比较绝对值会产生假象（例如泛内皮标记
CD34 在 Leading Edge 反而"更高"）。因此本脚本的主分析是：在同一个
SectionDataSet 内部，把目标结构与参考结构（默认 CT）配对，报告比值与精确符号
检验，这样切片间的批次差异会被抵消掉。

数据来源
--------
Allen Brain Map RMA API: https://api.brain-map.org/api/v2/data/query.json
Ivy GAP 资源论文: Puchalski RB, et al. An anatomic transcriptional atlas of
human glioblastoma. Science. 2018.

依赖
----
仅需 pandas（Anaconda 自带）。matplotlib 可选，用于出图。网络请求用标准库
urllib，无需安装 requests。

用法
----
    # 最小用法
    python ivygap_ish.py --genes PLVAP CAV1

    # 完整用法：指定结构、参考结构、输出目录、缓存目录
    python ivygap_ish.py \
        --genes PLVAP CAV1 CD34 ABCB1 ABCG2 ESM1 \
        --structures LE IT CT CTmvp CThbv CTpan \
        --ref CT \
        --outdir results \
        --cache cache \
        --plot

    # 查某个基因在不在 Ivy GAP 的探针库里（不下载数据）
    python ivygap_ish.py --genes MFSD2A --probe-only

输出
----
results/ivygap_ish_long.csv            逐条记录（基因×切片×结构）
results/ivygap_ish_structure_summary.csv  每个基因每个结构的 n / 中位数 / IQR
results/ivygap_ish_paired_vs_CT.csv    同切片内配对比较 + 精确符号检验
results/ivygap_ish_<gene>.png          可选，配对比较图

注意
----
- expression_energy 未按细胞密度归一化，跨结构绝对值比较不可靠，见上。
- 每张 ISH 切片只覆盖部分解剖结构，所以配对样本量通常远小于切片总数。
  结果里的 n_pairs 必须如实报告。
- density 接近 1.0 的记录通常是饱和/异常，默认标记为 qc_fail 并排除。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any

import pandas as pd

API = "https://api.brain-map.org/api/v2/data/query.json"
PRODUCT = "Glioblastoma"  # Ivy GAP 在 Allen API 里的 product 名

# Ivy GAP 解剖结构本体中与"核心—边缘"轴相关的结构
DEFAULT_STRUCTURES = ["LE", "IT", "CT", "CTmvp", "CThbv", "CTpan", "CTpnz", "CTne"]

STRUCTURE_LABELS = {
    "LE": "Leading Edge",
    "IT": "Infiltrating Tumor",
    "CT": "Cellular Tumor",
    "CTmvp": "Microvascular proliferation",
    "CThbv": "Hyperplastic blood vessels in cellular tumor",
    "CTpan": "Pseudopalisading cells around necrosis",
    "CTpnz": "Perinecrotic zone",
    "CTne": "Necrosis",
    "LEhbv": "Hyperplastic blood vessels in leading edge",
}


# --------------------------------------------------------------------------
# API 访问
# --------------------------------------------------------------------------
def _get(criteria: str, cache_dir: str | None, retries: int = 3,
         pause: float = 0.2) -> list[dict[str, Any]]:
    """执行一次 RMA 查询，返回 msg 列表。可选磁盘缓存，便于离线重跑。"""
    # 注意：不要把带反斜杠转义的表达式写进 f-string，Python < 3.12 会直接语法错误
    safe_chars = "[]$,:'"
    quoted = urllib.parse.quote(criteria, safe=safe_chars)
    url = f"{API}?criteria={quoted}&num_rows=all"

    cache_path = None
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        key = hashlib.sha1(url.encode("utf-8")).hexdigest()[:16]
        cache_path = os.path.join(cache_dir, f"{key}.json")
        if os.path.exists(cache_path):
            with open(cache_path, encoding="utf-8") as fh:
                return json.load(fh)

    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "ivygap-ish-script/1.0"}
            )
            with urllib.request.urlopen(req, timeout=60) as resp:
                payload = json.loads(resp.read().decode("utf-8"))
            if not payload.get("success"):
                raise RuntimeError(f"API returned success=false: {payload.get('msg')}")
            msg = payload.get("msg", [])
            if cache_path:
                with open(cache_path, "w", encoding="utf-8") as fh:
                    json.dump(msg, fh)
            time.sleep(pause)  # 对公共 API 客气一点
            return msg
        except (urllib.error.URLError, TimeoutError, RuntimeError, ValueError) as err:
            last_err = err
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError(f"查询失败（已重试 {retries} 次）: {url}\n{last_err}")


def list_section_datasets(gene: str, cache_dir: str | None) -> list[dict[str, Any]]:
    """列出某个基因在 Ivy GAP 里的所有 ISH SectionDataSet。"""
    criteria = (
        f"model::SectionDataSet,rma::criteria,"
        f"products[name$eq'{PRODUCT}'],genes[acronym$eq'{gene}'],"
        f"rma::include,genes,specimen"
    )
    return _get(criteria, cache_dir)


def structure_unionizes(dataset_id: int, cache_dir: str | None) -> list[dict[str, Any]]:
    """取某套切片按解剖结构汇总的 ISH 定量值。"""
    criteria = (
        f"model::StructureUnionize,rma::criteria,"
        f"section_data_set[id$eq{dataset_id}],rma::include,structure"
    )
    return _get(criteria, cache_dir)


# --------------------------------------------------------------------------
# 统计
# --------------------------------------------------------------------------
def sign_test_p(n_lower: int, n_pairs: int) -> float | None:
    """双侧精确符号检验（p=0.5 的二项检验）。不依赖 scipy。"""
    if n_pairs == 0:
        return None
    tail = min(n_lower, n_pairs - n_lower)
    cumulative = sum(math.comb(n_pairs, i) for i in range(tail + 1))
    return min(1.0, 2 * cumulative / (2 ** n_pairs))


# --------------------------------------------------------------------------
# 主流程
# --------------------------------------------------------------------------
def collect(genes: list[str], structures: list[str],
            cache_dir: str | None, verbose: bool = True) -> pd.DataFrame:
    """抓取并整理成 tidy 长表。"""
    rows: list[dict[str, Any]] = []
    wanted = set(structures)

    for gene in genes:
        datasets = list_section_datasets(gene, cache_dir)
        if verbose:
            print(f"[{gene}] Ivy GAP 中找到 {len(datasets)} 套 ISH 切片", flush=True)
        if not datasets:
            if verbose:
                print(f"[{gene}] 该基因不在 Ivy GAP 的探针库里，跳过", flush=True)
            continue

        for ds in datasets:
            ds_id = ds["id"]
            for unionize in structure_unionizes(ds_id, cache_dir):
                structure = unionize.get("structure") or {}
                acronym = structure.get("acronym")
                if acronym not in wanted:
                    continue
                density = unionize.get("expression_density")
                rows.append(
                    {
                        "gene": gene,
                        "section_data_set_id": ds_id,
                        "specimen_id": ds.get("specimen_id"),
                        "structure_acronym": acronym,
                        "structure_name": structure.get("name"),
                        "expression_density": density,
                        "expression_energy": unionize.get("expression_energy"),
                        "sum_expressing_pixels": unionize.get("sum_expressing_pixels"),
                        "sum_pixels": unionize.get("sum_pixels"),
                        # density 顶到 1.0 基本是饱和/分割异常
                        "qc_fail": bool(
                            density is None
                            or density >= 0.999
                            or not unionize.get("sum_pixels")
                        ),
                    }
                )

    return pd.DataFrame(rows)


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    """每个基因每个结构的描述统计。"""
    good = df[~df["qc_fail"]]
    out = (
        good.groupby(["gene", "structure_acronym"])["expression_energy"]
        .agg(
            n="count",
            median="median",
            q1=lambda s: s.quantile(0.25),
            q3=lambda s: s.quantile(0.75),
        )
        .reset_index()
    )
    out["structure_name"] = out["structure_acronym"].map(STRUCTURE_LABELS)
    return out.sort_values(["gene", "median"], ascending=[True, False])


def paired_vs_reference(df: pd.DataFrame, ref: str) -> pd.DataFrame:
    """同一张切片内，各结构 vs 参考结构的配对比较。"""
    good = df[~df["qc_fail"]]
    wide = good.pivot_table(
        index=["gene", "section_data_set_id"],
        columns="structure_acronym",
        values="expression_energy",
        aggfunc="first",
    )

    records = []
    for gene, block in wide.groupby(level="gene"):
        if ref not in block.columns:
            continue
        for structure in block.columns:
            if structure == ref:
                continue
            pair = block[[structure, ref]].dropna()
            n_pairs = len(pair)
            if n_pairs == 0:
                continue
            n_lower = int((pair[structure] < pair[ref]).sum())
            ratios = (pair[structure] / pair[ref]).replace(
                [float("inf"), float("-inf")], pd.NA
            ).dropna()
            records.append(
                {
                    "gene": gene,
                    "structure_acronym": structure,
                    "structure_name": STRUCTURE_LABELS.get(structure),
                    "reference": ref,
                    "n_pairs": n_pairs,
                    f"n_lower_than_{ref}": n_lower,
                    "median_ratio": float(ratios.median()) if len(ratios) else None,
                    "sign_test_p": sign_test_p(n_lower, n_pairs),
                }
            )

    out = pd.DataFrame(records)
    return out.sort_values(["gene", "median_ratio"]) if len(out) else out


def plot(df: pd.DataFrame, ref: str, outdir: str) -> None:
    """每个基因画一张同切片配对折线图。缺 matplotlib 时安静跳过。"""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("未安装 matplotlib，跳过出图", file=sys.stderr)
        return

    good = df[~df["qc_fail"]]
    for gene, block in good.groupby("gene"):
        wide = block.pivot_table(
            index="section_data_set_id",
            columns="structure_acronym",
            values="expression_energy",
            aggfunc="first",
        )
        order = [s for s in DEFAULT_STRUCTURES if s in wide.columns]
        if len(order) < 2:
            continue

        fig, ax = plt.subplots(figsize=(1.1 * len(order) + 2, 4))
        for _, row in wide[order].iterrows():
            ax.plot(order, row.values, marker="o", linewidth=1,
                    alpha=0.65, color="#777777")
        ax.plot(order, wide[order].median(), marker="s", linewidth=2.5,
                color="#1f77b4", label="median")
        ax.set_yscale("log")
        # 每个结构的 n 直接标在轴上：各结构覆盖的切片数不同，这一点必须让读者看见
        counts = wide[order].notna().sum()
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([f"{s}\nn={counts[s]}" for s in order])
        ax.set_ylabel("ISH expression energy (log scale)")
        ax.set_title(f"{gene} — Ivy GAP ISH by anatomic structure\n"
                     f"(lines = individual section datasets; ref = {ref})")
        ax.legend(frameon=False)
        fig.tight_layout()
        path = os.path.join(outdir, f"ivygap_ish_{gene}.png")
        fig.savefig(path, dpi=200)
        plt.close(fig)
        print(f"图已保存: {path}")


def plot_publication(df: pd.DataFrame, genes: list[str], ref: str, outdir: str,
                     structures: tuple[str, ...] = ("LE", "IT", "CT", "CTmvp"),
                     focus: str = "IT") -> None:
    """投稿版图：每个基因一个 panel，只画核心—边缘轴上的结构，同时输出 PNG 和 PDF。

    灰线 = 单张切片（同切片内配对），蓝线 = 中位数，每个结构下标注 n，
    panel 内标注 focus 结构相对参考结构的配对统计。
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("未安装 matplotlib，跳过出图", file=sys.stderr)
        return

    good = df[~df["qc_fail"]]
    panels = [g for g in genes if g in set(good["gene"])]
    if not panels:
        return

    plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.9,
                         "xtick.direction": "out", "ytick.direction": "out"})
    fig, axes = plt.subplots(1, len(panels), figsize=(4.3 * len(panels), 4.3))
    if len(panels) == 1:
        axes = [axes]

    for idx, (ax, gene) in enumerate(zip(axes, panels)):
        block = good[good["gene"] == gene]
        wide = block.pivot_table(
            index="section_data_set_id", columns="structure_acronym",
            values="expression_energy", aggfunc="first",
        )
        order = [s for s in structures if s in wide.columns]
        if len(order) < 2:
            continue
        xs = list(range(len(order)))

        for _, row in wide[order].iterrows():
            ax.plot(xs, row.values.astype(float), marker="o", markersize=4,
                    linewidth=0.9, color="0.68", alpha=0.85, zorder=1)
        ax.plot(xs, wide[order].median().values, marker="s", markersize=7,
                linewidth=2.4, color="#1f4e79", zorder=3, label="median")

        counts = wide[order].notna().sum()
        ax.set_xticks(xs)
        ax.set_xticklabels([f"{s}\n(n={counts[s]})" for s in order])
        ax.set_yscale("log")
        ax.set_xlim(-0.35, len(order) - 0.65)
        ax.set_title(gene, fontsize=12, fontweight="bold", pad=18)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.yaxis.grid(True, which="major", color="0.9", linewidth=0.7)
        ax.set_axisbelow(True)
        if idx == 0:
            ax.set_ylabel("ISH expression energy (log scale)")
            ax.legend(frameon=False, loc="lower right", fontsize=9)
        # panel 标号 A / B / ...
        ax.text(-0.17, 1.10, chr(ord("A") + idx), transform=ax.transAxes,
                fontsize=14, fontweight="bold", va="top")

        # focus vs ref 的同切片配对统计，直接写进图里
        if focus in wide.columns and ref in wide.columns:
            pair = wide[[focus, ref]].dropna()
            n_pairs = len(pair)
            if n_pairs:
                n_lower = int((pair[focus] < pair[ref]).sum())
                ratio = float((pair[focus] / pair[ref]).median())
                p = sign_test_p(n_lower, n_pairs)
                ax.text(
                    0.5, 1.015,
                    f"{focus}/{ref} median ratio {ratio:.2f}  "
                    f"({n_lower}/{n_pairs} slides lower, sign test P={p:.2f})",
                    transform=ax.transAxes, ha="center", va="bottom", fontsize=8.5,
                    color="0.25",
                )

    fig.tight_layout()
    for ext, kw in (("png", {"dpi": 300}), ("pdf", {})):
        path = os.path.join(outdir, f"ivygap_ish_publication.{ext}")
        fig.savefig(path, bbox_inches="tight", **kw)
        print(f"投稿版图已保存: {path}")
    plt.close(fig)

    print(
        "\n建议图注（按实际数字核对后使用）：\n"
        "  Region-resolved in situ hybridisation (Ivy GAP) for the indicated genes across\n"
        "  annotated glioblastoma anatomic structures. Grey lines, individual section\n"
        "  datasets; blue line, median. n below each structure gives the number of section\n"
        "  datasets covering it; structures are not covered by every section. Statistics are\n"
        "  within-slide paired comparisons against Cellular Tumor (exact two-sided sign test).\n"
        "  Results are shown as directional concordance with the transcriptomic analysis and\n"
        "  are not statistically significant at the available sample size."
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="从 Ivy GAP 提取 ISH 按解剖结构的定量值并做同切片配对比较",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--genes", nargs="+", required=True,
                        help="基因 symbol，例如 PLVAP CAV1 CD34")
    parser.add_argument("--structures", nargs="+", default=DEFAULT_STRUCTURES,
                        help=f"解剖结构缩写，默认 {' '.join(DEFAULT_STRUCTURES)}")
    parser.add_argument("--ref", default="CT",
                        help="配对比较的参考结构，默认 CT (Cellular Tumor)")
    parser.add_argument("--outdir", default="results", help="输出目录")
    parser.add_argument("--cache", default=None,
                        help="缓存目录；给了就把原始 JSON 存下来，重跑可离线")
    parser.add_argument("--plot", action="store_true", help="每个基因出一张探索用图")
    parser.add_argument("--publication", action="store_true",
                        help="出投稿版双 panel 图（PNG + PDF），只画核心—边缘轴上的结构")
    parser.add_argument("--pub-structures", nargs="+",
                        default=["LE", "IT", "CT", "CTmvp"],
                        help="投稿版图上要画的结构及顺序，默认 LE IT CT CTmvp")
    parser.add_argument("--probe-only", action="store_true",
                        help="只报告每个基因有几套 ISH 切片，不下载定量值")
    parser.add_argument("--keep-qc-fail", action="store_true",
                        help="保留饱和/异常记录（默认在统计时排除）")
    args = parser.parse_args()

    if args.probe_only:
        for gene in args.genes:
            n = len(list_section_datasets(gene, args.cache))
            verdict = f"{n} 套 ISH 切片" if n else "不在 Ivy GAP 探针库中"
            print(f"{gene:10s} {verdict}")
        return 0

    os.makedirs(args.outdir, exist_ok=True)

    df = collect(args.genes, args.structures, args.cache)
    if df.empty:
        print("没有取到任何数据：检查基因名拼写，或该基因不在 Ivy GAP 探针库中。",
              file=sys.stderr)
        return 1

    if args.keep_qc_fail:
        df = df.assign(qc_fail=False)

    n_fail = int(df["qc_fail"].sum())
    if n_fail:
        print(f"标记为 QC 异常并排除的记录: {n_fail} 条", file=sys.stderr)

    long_path = os.path.join(args.outdir, "ivygap_ish_long.csv")
    df.to_csv(long_path, index=False)
    print(f"长表已保存: {long_path}  ({len(df)} 行)")

    summary = summarise(df)
    summary_path = os.path.join(args.outdir, "ivygap_ish_structure_summary.csv")
    summary.to_csv(summary_path, index=False)
    print(f"结构汇总已保存: {summary_path}")

    paired = paired_vs_reference(df, args.ref)
    paired_path = os.path.join(args.outdir, f"ivygap_ish_paired_vs_{args.ref}.csv")
    paired.to_csv(paired_path, index=False)
    print(f"配对比较已保存: {paired_path}")

    if len(paired):
        print("\n=== 同切片内配对比较（主结果） ===")
        print(paired.to_string(index=False))
        print("\n提示: n_pairs 往往远小于切片总数，因为每张 ISH 切片只覆盖部分解剖"
              "结构。写论文时请如实报告 n_pairs，不要只报 p 值。")

    if args.plot:
        plot(df, args.ref, args.outdir)

    if args.publication:
        plot_publication(df, args.genes, args.ref, args.outdir,
                         tuple(args.pub_structures))

    return 0


if __name__ == "__main__":
    sys.exit(main())
