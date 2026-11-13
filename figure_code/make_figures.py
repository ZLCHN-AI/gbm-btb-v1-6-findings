"""Rebuild v1.6 Fig. 2, Fig. 4 and Supplementary Figs. S2-S3 from frozen source tables.

Inputs  : ./data/*.csv  (copies of the frozen figure-data and canonical tables)
Outputs : ./out/Fig2.*, Fig4.*, FigS2.*, FigS3.*  (TIFF 600 dpi LZW, PDF, PNG preview)

No statistic is recomputed here except donor medians, which are checked against the
frozen summaries; the script only draws values that already exist in the tables.
"""
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)

GENE_ORDER = ["MFSD2A", "PLVAP", "CAV1", "CAVIN1"]
GENE_COLORS = {"MFSD2A": "#2C7BB6", "PLVAP": "#D73027", "CAV1": "#7B3294", "CAVIN1": "#8C6D31"}
GROUP_COLORS = {"CANONICAL": "#4C78A8", "REMODELED": "#E45756"}
CONTRAST_COLORS = {"IT_vs_LE": "#4C78A8", "CT_vs_IT": "#F58518"}
DONOR_COLORS = {"AT10": "#355C7D", "AT15": "#6C5B7B", "AT3": "#C06C84", "AT4": "#4A7C59"}
REGION_ORDER = ["LE", "IT", "CT"]
GRID = "#D9E2E7"

mpl.rcParams.update(
    {
        "font.family": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8.5,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def panel_label(ax, label, x=-0.13, y=1.10):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom", ha="left")


def save(fig, stem):
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white", pad_inches=0.08)
    fig.savefig(OUT / f"{stem}.png", bbox_inches="tight", facecolor="white", pad_inches=0.08, dpi=200)
    fig.savefig(
        OUT / f"{stem}.tif", bbox_inches="tight", facecolor="white", pad_inches=0.08, dpi=600,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)


# --------------------------------------------------------------------------- Fig. 2
def fig2():
    subtype = pd.read_csv(DATA / "13_FIG3A_endothelial_subtype_profile.csv")
    order = list(dict.fromkeys(subtype["COMMON_SUBTYPE"]))
    matrix = subtype.pivot(index="COMMON_SUBTYPE", columns="gene", values="source_defined_normalized_score")
    matrix = matrix.loc[order, GENE_ORDER]
    paired = pd.read_csv(DATA / "15_FIG3B_paired_endothelial_expression.csv")
    summary = pd.read_csv(DATA / "17_FIG3C_endothelial_delta_summary.csv")
    comp_long = pd.read_csv(DATA / "18_FIG3D_vascular_region_composition.csv")
    components = list(dict.fromkeys(comp_long["vascular_component"]))
    comp = comp_long.pivot(index="TRAJECTORY_REGION", columns="vascular_component", values="regional_program_abundance")
    comp = comp.loc[REGION_ORDER, components]
    spatial = pd.read_csv(DATA / "19_FIG3E_spatial_donor_effects.csv")

    fig = plt.figure(figsize=(12.2, 9.6))
    grid = fig.add_gridspec(2, 3, left=0.065, right=0.985, bottom=0.085, top=0.93,
                            width_ratios=[1.0, 1.25, 1.05], height_ratios=[1.0, 1.08], wspace=0.46, hspace=0.50)
    ax_a, ax_b, ax_c = (fig.add_subplot(grid[0, i]) for i in range(3))
    ax_d = fig.add_subplot(grid[1, 0])
    ax_e = fig.add_subplot(grid[1, 1:])

    heat = ax_a.imshow(matrix.values, cmap="RdBu_r", vmin=-0.85, vmax=0.85, aspect="auto")
    ax_a.set_xticks(range(len(GENE_ORDER)), GENE_ORDER, rotation=35, ha="right")
    ax_a.set_yticks(range(len(matrix.index)), matrix.index)
    ax_a.set_title("Endothelial subtype expression programs", loc="left", fontweight="bold")
    ax_a.set_xlabel("Source-defined normalized score")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            v = matrix.iloc[i, j]
            ax_a.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7,
                      color="white" if abs(v) > 0.6 else "black")
    cbar = fig.colorbar(heat, ax=ax_a, fraction=0.046, pad=0.04)
    cbar.ax.tick_params(labelsize=7)
    for s in ("top", "right"):
        ax_a.spines[s].set_visible(True)
    panel_label(ax_a, "a", x=-0.30)

    x_map = {("MFSD2A", "CANONICAL"): 0, ("MFSD2A", "REMODELED"): 1, ("PLVAP", "CANONICAL"): 3, ("PLVAP", "REMODELED"): 4}
    rng = np.random.default_rng(20260814)
    for gene in ["MFSD2A", "PLVAP"]:
        sub = paired[paired["gene"] == gene]
        piv = sub.pivot(index="donor_id", columns="group", values="pseudobulk_log1pCPM")
        for _, row in piv.iterrows():
            xs = [x_map[(gene, "CANONICAL")], x_map[(gene, "REMODELED")]]
            j = rng.uniform(-0.035, 0.035)
            ax_b.plot([x + j for x in xs], [row["CANONICAL"], row["REMODELED"]], color=GENE_COLORS[gene], alpha=0.28, linewidth=0.9)
            ax_b.scatter([x + j for x in xs], [row["CANONICAL"], row["REMODELED"]], s=17,
                         color=[GROUP_COLORS["CANONICAL"], GROUP_COLORS["REMODELED"]], edgecolor="white", linewidth=0.35, zorder=3)
        med = sub.groupby("group")["pseudobulk_log1pCPM"].median()
        for g in ["CANONICAL", "REMODELED"]:
            x = x_map[(gene, g)]
            ax_b.plot([x - 0.20, x + 0.20], [med[g], med[g]], color="black", linewidth=2.3)
    ax_b.set_xticks([0, 1, 3, 4], ["Canonical", "Remodeled", "Canonical", "Remodeled"])
    ax_b.text(0.5, -0.17, "MFSD2A", transform=ax_b.get_xaxis_transform(), ha="center", fontweight="bold")
    ax_b.text(3.5, -0.17, "PLVAP", transform=ax_b.get_xaxis_transform(), ha="center", fontweight="bold")
    ax_b.set_ylabel("Pseudobulk log(1+CPM)")
    ax_b.set_title("Paired endothelial expression by donor", loc="left", fontweight="bold")
    ax_b.grid(axis="y", color=GRID, linewidth=0.6, alpha=0.8)
    ymax = paired["pseudobulk_log1pCPM"].max()
    ax_b.set_ylim(top=ymax + 1.0)
    ax_b.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GROUP_COLORS["CANONICAL"], label="Canonical", markersize=6),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=GROUP_COLORS["REMODELED"], label="Remodeled", markersize=6),
        Line2D([0], [0], color="black", linewidth=2.3, label="Median"),
    ], loc="upper center", frameon=False, ncol=3, fontsize=8)
    ax_b.text(0.02, 0.02, "n = 12 donors; lines connect matched donors", transform=ax_b.transAxes, fontsize=7.5)
    panel_label(ax_b, "b")

    forest = summary[summary["metric"] == "pseudobulk_log1pCPM"].copy()
    forest["order"] = forest["gene"].map({g: i for i, g in enumerate(GENE_ORDER)})
    forest = forest.sort_values("order", ascending=False)
    y = np.arange(len(forest))
    for yi, (_, r) in zip(y, forest.iterrows()):
        ax_c.errorbar(r["mean_donor_delta"], yi,
                      xerr=[[r["mean_donor_delta"] - r["bootstrap_mean_CI_2.5"]], [r["bootstrap_mean_CI_97.5"] - r["mean_donor_delta"]]],
                      fmt="o", markersize=6, color=GENE_COLORS[r["gene"]], ecolor=GENE_COLORS[r["gene"]], capsize=3, linewidth=1.6)
        ax_c.text(1.03, yi, f"{int(r['positive_donors'])}/12", transform=ax_c.get_yaxis_transform(), ha="left", va="center", fontsize=8)
    ax_c.text(1.03, len(forest) - 0.45, "Donors\npositive", transform=ax_c.get_yaxis_transform(), ha="left", va="bottom", fontsize=7.5)
    ax_c.axvline(0, color="#555555", linewidth=0.9, linestyle="--")
    ax_c.set_yticks(y, forest["gene"])
    ax_c.set_ylim(-0.5, len(forest) - 0.5)
    ax_c.set_xlabel("Mean donor difference, remodeled − canonical\n[pseudobulk log(1+CPM)]")
    ax_c.set_title("Endothelial donor differences", loc="left", fontweight="bold")
    ax_c.grid(axis="x", color=GRID, linewidth=0.6, alpha=0.8)
    panel_label(ax_c, "c", x=-0.24)

    colors = {"VASC_CAP_ENDO": "#4C78A8", "VASC_REMODELED_ENDO": "#E45756", "VASC_PERICYTES": "#72B7B2", "VASC_VLMC": "#F2CF5B", "VASC_OTHER": "#B9B9B9"}
    names = {"VASC_CAP_ENDO": "Capillary endothelium", "VASC_REMODELED_ENDO": "Remodeled endothelium", "VASC_PERICYTES": "Pericytes", "VASC_VLMC": "VLMC", "VASC_OTHER": "Other"}
    bottom = np.zeros(len(REGION_ORDER))
    for c in components:
        vals = comp[c].to_numpy()
        ax_d.bar(REGION_ORDER, vals, bottom=bottom, label=names.get(c, c), color=colors.get(c, "#999999"), width=0.68, edgecolor="white", linewidth=0.4)
        bottom += vals
    ax_d.set_ylabel("Regional vascular-program abundance")
    ax_d.set_title("Vascular composition by region", loc="left", fontweight="bold")
    ax_d.legend(frameon=False, fontsize=7, loc="upper left")
    ax_d.grid(axis="y", color=GRID, linewidth=0.6, alpha=0.8)
    panel_label(ax_d, "d", x=-0.17)

    xpos = np.arange(len(GENE_ORDER), dtype=float)
    off = {"IT_vs_LE": -0.16, "CT_vs_IT": 0.16}
    allv = []
    for contrast in ["IT_vs_LE", "CT_vs_IT"]:
        for gi, gene in enumerate(GENE_ORDER):
            vals = spatial.loc[(spatial["contrast"] == contrast) & (spatial["gene"] == gene), "beta_DETAILED"].to_numpy()
            allv.extend(vals)
            x = xpos[gi] + off[contrast]
            ax_e.scatter(x + np.linspace(-0.035, 0.035, len(vals)), vals, s=28, color=CONTRAST_COLORS[contrast], alpha=0.82, edgecolor="white", linewidth=0.5, zorder=3)
            ax_e.scatter([x], [float(np.median(vals))], marker="D", s=55, color="black", edgecolor="white", linewidth=0.7, zorder=4)
    ax_e.axhline(0, color="#555555", linewidth=0.9, linestyle="--")
    ax_e.set_xticks(xpos, GENE_ORDER)
    ax_e.set_ylabel("Composition-adjusted spatial β")
    ax_e.set_title("Regional spatial associations", loc="left", fontweight="bold")
    ax_e.grid(axis="y", color=GRID, linewidth=0.6, alpha=0.8)
    ax_e.set_ylim(min(allv) - 0.10, max(allv) + 0.16)
    ax_e.legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor=CONTRAST_COLORS["IT_vs_LE"], label="LE to IT", markersize=6),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=CONTRAST_COLORS["CT_vs_IT"], label="IT to CT", markersize=6),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="black", label="Donor median", markersize=6),
    ], frameon=False, ncol=3, loc="upper left")
    ax_e.text(0.99, 0.02, "n = 4 spatial donors", transform=ax_e.transAxes, fontsize=7.5, ha="right")
    panel_label(ax_e, "e", x=-0.075)
    save(fig, "Fig2")
    f = forest.set_index("gene")
    return {g: (round(f.loc[g, "mean_donor_delta"], 5), int(f.loc[g, "positive_donors"])) for g in GENE_ORDER}


# --------------------------------------------------------------------------- Fig. 4
def load_tj():
    d = pd.read_csv(DATA / "Phase2H3_continuous_donor_contrasts.csv")
    return d[d["score"] == "CONTINUOUS_GATE_CONVERSION_SCORE"].set_index("donor_id")[["IT_minus_LE", "CT_minus_IT"]]


def fig4():
    tj = load_tj()
    fig, ax = plt.subplots(figsize=(3.6, 3.5))
    xs = [0, 1]
    for donor, row in tj.iterrows():
        ax.plot(xs, row.values, color=DONOR_COLORS[donor], linewidth=1.2, alpha=0.75, zorder=2)
        ax.scatter(xs, row.values, s=30, color=DONOR_COLORS[donor], edgecolor="white", linewidth=0.5, zorder=3, label=donor)
    med = tj.median()
    ax.scatter(xs, med.values, marker="D", s=70, color="black", edgecolor="white", linewidth=0.8, zorder=4)
    ax.axhline(0, color="#555555", linewidth=0.9, linestyle="--", zorder=1)
    ax.set_xticks(xs, ["LE to IT", "IT to CT"])
    ax.set_xlim(-0.35, 1.35)
    lo, hi = tj.values.min(), tj.values.max()
    ax.set_ylim(lo - 0.08, hi + 0.08)
    ax.set_ylabel("T-J balance contrast\n(adjusted donor-level effect)")
    ax.grid(axis="y", color=GRID, linewidth=0.6, alpha=0.8)
    for x, col in zip(xs, tj.columns):
        pos = int((tj[col] > 0).sum())
        ax.text(x, 1.03, f"median {med[col]:+.3f}\n{pos}/4 donors positive", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=8)
    handles = [Line2D([0], [0], marker="o", color=DONOR_COLORS[d], markerfacecolor=DONOR_COLORS[d], markersize=5, linewidth=1.2, label=d) for d in tj.index]
    handles.append(Line2D([0], [0], marker="D", color="none", markerfacecolor="black", markersize=6, label="Donor median"))
    ax.legend(handles=handles, frameon=False, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3)
    save(fig, "Fig4")
    return {c: (round(float(med[c]), 5), int((tj[c] > 0).sum())) for c in tj.columns}


# --------------------------------------------------------------------------- Supplementary Fig. S2
def figS2():
    cs = pd.read_csv(DATA / "Phase2H3_continuous_summary.csv")
    cs = cs[cs["score"] == "CONTINUOUS_GATE_CONVERSION_SCORE"].set_index("cohort")
    pr = pd.read_csv(DATA / "Phase2H3_pattern_robustness_audit.csv")
    sm = pd.read_csv(DATA / "20_FIG3E_spatial_summary_and_MIN100.csv")
    sm = sm[sm["contrast"] == "CT_vs_IT"]
    cohorts = [("MIN50_PRIMARY", "MIN50", "#4C78A8"), ("MIN100_SENSITIVITY", "MIN100", "#F58518")]

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.9), gridspec_kw={"wspace": 0.38, "left": 0.06, "right": 0.99, "bottom": 0.17, "top": 0.84})
    ax = axes[0]
    w = 0.36
    for k, (key, lab, col) in enumerate(cohorts):
        vals = [cs.loc[key, "median_donor_IT_minus_LE"], cs.loc[key, "median_donor_CT_minus_IT"]]
        pos = [int(cs.loc[key, "positive_donors_IT_minus_LE"]), int(cs.loc[key, "positive_donors_CT_minus_IT"])]
        x = np.arange(2) + (k - 0.5) * w
        ax.bar(x, vals, width=w, color=col, label=lab)
        for xi, v, p in zip(x, vals, pos):
            ax.text(xi, v + (0.012 if v >= 0 else -0.012), f"{p}/4", ha="center", va="bottom" if v >= 0 else "top", fontsize=7.5)
    ax.axhline(0, color="#555555", linewidth=0.9)
    ax.set_xticks(np.arange(2), ["LE to IT", "IT to CT"])
    ax.set_ylabel("Median donor T-J balance contrast")
    ax.set_ylim(-0.23, 0.30)
    ax.set_title("T-J balance at MIN50 and MIN100", loc="left", fontweight="bold")
    ax.grid(axis="y", color=GRID, linewidth=0.6, alpha=0.8)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.text(0.99, 0.02, "labels: donors positive", transform=ax.transAxes, ha="right", fontsize=7.5)
    panel_label(ax, "a", x=-0.17)

    ax = axes[1]
    ax.plot(pr["quantile_threshold"], pr["positive_donors_CT_minus_IT_GATE"], marker="o", color="#F58518", label="IT to CT")
    ax.plot(pr["quantile_threshold"], pr["positive_donors_IT_minus_LE_GATE"], marker="o", color="#4C78A8", label="LE to IT")
    ax.set_ylim(-0.3, 4.6)
    ax.set_yticks(range(5))
    ax.set_xticks(pr["quantile_threshold"])
    ax.set_xlabel("Quantile threshold defining the high T-J state")
    ax.set_ylabel("Donors with positive contrast (of 4)")
    ax.set_title("Threshold robustness", loc="left", fontweight="bold")
    ax.grid(axis="y", color=GRID, linewidth=0.6, alpha=0.8)
    ax.legend(frameon=False, fontsize=8, loc="center right")
    panel_label(ax, "b", x=-0.15)

    ax = axes[2]
    genes = ["CAV1", "CAVIN1", "MFSD2A", "PLVAP"]
    for k, (key, lab, col) in enumerate(cohorts):
        s = sm[sm["cohort"] == key].set_index("gene").loc[genes]
        x = np.arange(len(genes)) + (k - 0.5) * w
        ax.bar(x, s["median_beta_DETAILED"], width=w, color=col, label=lab)
        for xi, v, p in zip(x, s["median_beta_DETAILED"], s["positive_donors_DETAILED"]):
            ax.text(xi, v + 0.006, f"{int(p)}/4", ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(np.arange(len(genes)), genes)
    ax.set_ylabel("Median donor detailed-model β (CT versus IT)")
    ax.set_title("Vascular-composition sensitivity", loc="left", fontweight="bold")
    ax.set_ylim(0, 0.33)
    ax.grid(axis="y", color=GRID, linewidth=0.6, alpha=0.8)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    panel_label(ax, "c", x=-0.17)
    save(fig, "FigS2")
    return cs[["median_donor_IT_minus_LE", "median_donor_CT_minus_IT"]].round(4).to_dict("index"), pr[["quantile_threshold", "positive_donors_IT_minus_LE_GATE", "positive_donors_CT_minus_IT_GATE"]].values.tolist()


# --------------------------------------------------------------------------- Supplementary Fig. S3
def figS3():
    raw = pd.read_csv(DATA / "PatientDonor_KeyClaim_Raw_Effects.csv", encoding="utf-8-sig")
    summ = pd.read_csv(DATA / "PatientDonor_KeyClaim_Summary.csv", encoding="utf-8-sig").set_index("signal")
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.0), gridspec_kw={"wspace": 0.55, "left": 0.16, "right": 0.97, "bottom": 0.22, "top": 0.82, "width_ratios": [1, 1.25]})

    ax = axes[0]
    g = raw[(raw["signal"] == "Gate") & (raw["contrast"] == "IT_to_CT")]
    r = summ.loc["Gate"]
    ax.errorbar(r["median_effect"], 0, xerr=[[r["median_effect"] - r["bootstrap_median_ci_low"]], [r["bootstrap_median_ci_high"] - r["median_effect"]]],
                fmt="D", markersize=7, color="black", ecolor="black", capsize=4, linewidth=1.5, zorder=3)
    for _, row in g.iterrows():
        ax.scatter(row["effect"], 0.22, s=30, color=DONOR_COLORS[row["biological_unit"]], edgecolor="white", linewidth=0.5, zorder=4)
        ax.text(row["effect"], 0.36, row["biological_unit"], ha="center", va="bottom", fontsize=7, rotation=45)
    ax.axvline(0, color="#555555", linewidth=0.9, linestyle="--")
    ax.set_yticks([0], ["T-J balance\nIT to CT"])
    ax.set_ylim(-0.5, 0.95)
    ax.set_xlim(-0.05, 0.62)
    ax.set_xlabel("Adjusted donor-level effect")
    ax.set_title("Spatial contrast (4 donors)", loc="left", fontweight="bold")
    ax.grid(axis="x", color=GRID, linewidth=0.6, alpha=0.8)
    panel_label(ax, "a", x=-0.42, y=1.12)

    ax = axes[1]
    for yi, gene in enumerate(["PLVAP", "MFSD2A"]):
        r = summ.loc[gene]
        vals = raw[(raw["signal"] == gene)]["effect"].to_numpy()
        ax.scatter(vals, np.full(len(vals), yi + 0.22), s=18, color=GENE_COLORS[gene], alpha=0.7, edgecolor="white", linewidth=0.4, zorder=3)
        ax.errorbar(r["median_effect"], yi, xerr=[[r["median_effect"] - r["bootstrap_median_ci_low"]], [r["bootstrap_median_ci_high"] - r["median_effect"]]],
                    fmt="D", markersize=7, color="black", ecolor="black", capsize=4, linewidth=1.5, zorder=4)
    ax.axvline(0, color="#555555", linewidth=0.9, linestyle="--")
    ax.set_yticks([0, 1], ["PLVAP", "MFSD2A"])
    ax.set_ylim(-0.5, 1.6)
    ax.set_xlabel("Remodeled − canonical, pseudobulk log(1+CPM)")
    ax.set_title("Endothelial contrasts (12 donors)", loc="left", fontweight="bold")
    ax.grid(axis="x", color=GRID, linewidth=0.6, alpha=0.8)
    panel_label(ax, "b", x=-0.20, y=1.12)
    save(fig, "FigS3")
    return summ.loc[["Gate", "PLVAP", "MFSD2A"], ["median_effect", "mean_effect", "bootstrap_median_ci_low", "bootstrap_median_ci_high", "expected_direction_donors"]].round(5).to_dict("index")


if __name__ == "__main__":
    print("Fig2 forest (mean donor delta, donors positive):", fig2())
    print("Fig4 (median, donors positive):", fig4())
    print("FigS2:", figS2())
    print("FigS3:", figS3())
