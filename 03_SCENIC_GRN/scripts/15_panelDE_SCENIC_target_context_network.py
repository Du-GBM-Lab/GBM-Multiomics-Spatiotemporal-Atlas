#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path
from math import comb, log10

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import networkx as nx
import numpy as np
import pandas as pd


SUBTYPE_COLORS = {
    "MES-V": "#C85B4B",
    "ATAC": "#D58B2A",
    "SCENIC": "#526D9E",
    "FOSL1": "#9C2F2F",
    "Other": "#B9C0C9",
}


FUNCTIONAL_SETS = {
    "ECM remodeling": {
        "COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL6A1", "COL6A2",
        "FN1", "POSTN", "SPARC", "THBS1", "THBS2", "TNC", "VCAN", "DCN", "LOX", "LOXL2",
        "MMP2", "MMP9", "MMP14", "ADAM12", "ADAMTS1", "ADAMTS4", "ADAMTS9",
    },
    "Cell adhesion / invasion": {
        "CD44", "CDH2", "VIM", "ITGA5", "ITGA6", "ITGAV", "ITGB1", "ITGB3", "CAV1",
        "ANXA2", "ACTA2", "ACTN1", "RAC1", "RHOA", "LAMB1", "LAMC1", "LAMA4", "SERPINE1",
        "PLAUR", "PLAU", "MMP2", "MMP9", "MMP14",
    },
    "Hypoxia / vascular niche": {
        "VEGFA", "CA9", "NDRG1", "ADM", "HILPDA", "ANGPT2", "EDN1", "PDGFB", "EPAS1",
        "HIF1A", "SLC2A1", "BNIP3", "LOX", "CXCR4", "KDR", "FLT1", "VWF", "PECAM1",
    },
    "Inflammatory / stress response": {
        "IL6", "IL1B", "TNF", "CXCL8", "CXCL2", "CXCL3", "CCL2", "PTGS2", "NFKBIA",
        "NFKB1", "STAT3", "IRF1", "SOCS3", "JUN", "JUNB", "FOS", "FOSL1", "FOSL2",
        "ATF3", "DUSP1", "DUSP5", "EGR1",
    },
    "AP-1 / immediate-early": {
        "FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "ATF3", "BATF",
        "EGR1", "EGR2", "DUSP1", "DUSP5", "IER2", "IER3", "NR4A1", "KLF6",
    },
    "PLAU-PLAUR axis context": {
        "PLAUR", "PLAU", "SERPINE1", "SERPINB2", "ANXA2", "MMP2", "MMP9", "MMP14",
        "ITGB1", "CD44", "VIM", "FN1", "POSTN", "CDH2",
    },
}


def hypergeom_sf(x: int, n_category: int, n_draw: int, n_background: int) -> float:
    """P[X >= x] for hypergeometric; pure Python to avoid extra dependency assumptions."""
    if x <= 0:
        return 1.0
    max_k = min(n_category, n_draw)
    denom = comb(n_background, n_draw)
    return sum(
        comb(n_category, k) * comb(n_background - n_category, n_draw - k) / denom
        for k in range(x, max_k + 1)
    )


def bh_fdr(pvals: list[float]) -> list[float]:
    arr = np.asarray(pvals, dtype=float)
    order = np.argsort(arr)
    ranked = arr[order]
    out = np.empty_like(ranked)
    n = len(arr)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        val = ranked[i] * n / (i + 1)
        prev = min(prev, val)
        out[i] = prev
    res = np.empty_like(out)
    res[order] = out
    return res.tolist()


def add_box(ax, xy, text, fc, ec="#333333", w=1.25, h=0.34, fs=8.2, weight="normal"):
    x, y = xy
    box = FancyBboxPatch(
        (x - w / 2, y - h / 2),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.05",
        fc=fc,
        ec=ec,
        lw=0.9,
        mutation_aspect=1,
    )
    ax.add_patch(box)
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, fontweight=weight, color="#1F2428")


def add_arrow(ax, start, end, color="#8A8A8A", lw=1.0, alpha=0.9, ls="-", rad=0.0):
    arr = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=8,
        lw=lw,
        color=color,
        alpha=alpha,
        linestyle=ls,
        connectionstyle=f"arc3,rad={rad}",
        shrinkA=12,
        shrinkB=12,
    )
    ax.add_patch(arr)


def main() -> None:
    scenic_dir = Path(__file__).resolve().parents[1]
    project_root = scenic_dir.parents[2]
    final_root = (
        project_root
        / "图片表格"
        / "05_发育时间_TF_ATAC验证"
        / "02_TF_regulon_SCENIC"
    )

    outputs = scenic_dir / "outputs"
    script_name = Path(__file__).name

    candidate_dir = final_root / "正文图候选"
    main_dir = final_root / "正文图"
    source_dir = final_root / "source_data" / "DE_SCENIC_target_context_network"
    record_dir = final_root / "记录"
    for d in [candidate_dir, main_dir, source_dir, record_dir]:
        d.mkdir(parents=True, exist_ok=True)

    grn = pd.read_csv(outputs / "scenic_grn_edgelist.csv")
    grn["TF"] = grn["TF"].astype(str)
    grn["target"] = grn["target"].astype(str)
    grn["importance"] = pd.to_numeric(grn["importance"], errors="coerce").fillna(0)

    overlap = pd.read_csv(
        final_root
        / "source_data"
        / "S_SCENIC_ATAC_FOSL1_prioritization_intersection"
        / "PLAUR_upstream_and_MESV_top20_overlap.csv"
    )
    candidate_tfs = overlap.sort_values("MESV_RSS_rank")["TF"].tolist()

    # Use top targets per TF to make target context focused and reproducible.
    top_n = 500
    rows = []
    target_rows = []
    background = set(grn["target"].dropna().unique())
    for tf in candidate_tfs:
        sub = grn.loc[grn["TF"].eq(tf)].sort_values("importance", ascending=False)
        top_targets = sub.head(top_n)["target"].dropna().astype(str).tolist()
        top_set = set(top_targets)
        for rank, (_, r) in enumerate(sub.head(top_n).iterrows(), start=1):
            target_rows.append(
                {
                    "TF": tf,
                    "target": r["target"],
                    "target_rank_by_importance": rank,
                    "importance": r["importance"],
                }
            )
        for context, genes in FUNCTIONAL_SETS.items():
            genes_bg = set(genes) & background
            hits = sorted(top_set & genes_bg)
            p = hypergeom_sf(len(hits), len(genes_bg), len(top_set), len(background))
            rows.append(
                {
                    "TF": tf,
                    "functional_context": context,
                    "n_hits": len(hits),
                    "n_context_genes_in_background": len(genes_bg),
                    "n_top_targets": len(top_set),
                    "p_value": p,
                    "hit_genes": ";".join(hits),
                }
            )

    enrich = pd.DataFrame(rows)
    enrich["FDR"] = bh_fdr(enrich["p_value"].tolist())
    enrich["minus_log10_FDR"] = -np.log10(enrich["FDR"].clip(lower=1e-12))
    enrich.to_csv(source_dir / "D_candidate_TF_target_functional_context.csv", index=False)
    pd.DataFrame(target_rows).to_csv(source_dir / "D_candidate_TF_top500_targets.csv", index=False)

    # Panel D: target functional context dot plot.
    tf_order = candidate_tfs
    context_order = [
        "ECM remodeling",
        "Cell adhesion / invasion",
        "Hypoxia / vascular niche",
        "Inflammatory / stress response",
        "AP-1 / immediate-early",
        "PLAU-PLAUR axis context",
    ]
    plot_df = enrich.copy()
    plot_df["x"] = plot_df["TF"].map({v: i for i, v in enumerate(tf_order)})
    plot_df["y"] = plot_df["functional_context"].map({v: i for i, v in enumerate(context_order[::-1])})
    plot_df["plot_score"] = plot_df["minus_log10_FDR"].clip(0, 6)
    plot_df["size"] = 18 + plot_df["n_hits"] * 13

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
        }
    )

    fig, ax = plt.subplots(figsize=(7.4, 3.85))
    sc = ax.scatter(
        plot_df["x"],
        plot_df["y"],
        s=plot_df["size"],
        c=plot_df["plot_score"],
        cmap="YlOrRd",
        vmin=0,
        vmax=max(3, plot_df["plot_score"].max()),
        edgecolor="#4F4F4F",
        linewidth=0.35,
    )
    ax.set_xticks(range(len(tf_order)))
    labels = []
    for tf in tf_order:
        if tf == "FOSL1":
            labels.append(r"$\bf{FOSL1}$")
        else:
            labels.append(tf)
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8.6)
    ax.set_yticks(range(len(context_order)))
    ax.set_yticklabels(context_order[::-1], fontsize=8.6)
    ax.set_xlim(-0.55, len(tf_order) - 0.45)
    ax.set_ylim(-0.55, len(context_order) - 0.45)
    ax.grid(axis="both", color="#E8E8E8", lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title("Candidate TF target-gene functional context", fontsize=10.5, pad=10, weight="bold")
    ax.set_xlabel("PLAUR-upstream TFs also ranked in MES-V top20 RSS", fontsize=8.8)
    ax.text(
        0.99,
        -0.29,
        "Curated contexts tested on top 500 SCENIC targets per TF; exploratory context, not causal proof.",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.8,
        color="#5B5B5B",
    )
    cbar = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.025)
    cbar.set_label("-log10(FDR)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    legend_hits = [2, 5, 10]
    handles = [
        ax.scatter([], [], s=18 + h * 13, facecolor="white", edgecolor="#555555", linewidth=0.4)
        for h in legend_hits
    ]
    ax.legend(
        handles,
        [str(h) for h in legend_hits],
        title="hit genes",
        loc="upper left",
        bbox_to_anchor=(1.22, 0.44),
        frameon=False,
        fontsize=7,
        title_fontsize=7,
        labelspacing=0.8,
    )
    fig.tight_layout()
    for out_dir in [candidate_dir, main_dir]:
        fig.savefig(out_dir / "D_候选TF_target功能语义点图.pdf", bbox_inches="tight")
    fig.savefig(candidate_dir / "D_候选TF_target功能语义点图.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Panel E: evidence network, not causal GRN.
    net_rows = []
    for _, r in overlap.iterrows():
        net_rows.append(
            {
                "TF": r["TF"],
                "PLAUR_upstream_rank": r["rank"],
                "PLAUR_edge_importance": r["importance"],
                "MESV_RSS_rank": r["MESV_RSS_rank"],
                "MESV_RSS": r["MESV_RSS"],
                "ATAC_motif_context": bool(r["has_ATAC_MESlineage_characteristic_motif"]),
                "prioritized": r["TF"] == "FOSL1",
            }
        )
    pd.DataFrame(net_rows).to_csv(source_dir / "E_prioritization_network_nodes.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.35, 4.35))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")

    # Evidence boxes
    add_box(ax, (1.35, 5.25), "SCENIC\nPLAUR upstream TFs\nn = 86", "#E8EEF9", "#526D9E", w=1.65, h=0.78, fs=8)
    add_box(ax, (1.35, 3.05), "MES-V top20\nRSS regulons", "#F6E9E5", "#C85B4B", w=1.65, h=0.66, fs=8)
    add_box(ax, (7.3, 3.05), "ATAC MES-lineage\nmotif program", "#FAEBD1", "#D58B2A", w=1.75, h=0.66, fs=8)
    add_box(ax, (8.85, 1.25), "PLAUR promoter-proximal\nAP-1 anchor", "#FAEBD1", "#D58B2A", w=2.0, h=0.66, fs=7.6)
    add_box(ax, (8.85, 5.25), "FOSL1 prioritized\nfor wet-lab validation", "#F4D6D1", "#9C2F2F", w=2.0, h=0.66, fs=8.2, weight="bold")

    y_positions = np.linspace(4.65, 1.75, len(tf_order))
    tf_pos = {}
    for tf, y in zip(tf_order, y_positions):
        x = 4.0
        tf_pos[tf] = (x, y)
        row = overlap.loc[overlap["TF"].eq(tf)].iloc[0]
        color = "#F4D6D1" if tf == "FOSL1" else "#F3F4F6"
        edge = "#9C2F2F" if tf == "FOSL1" else "#A9AEB6"
        add_box(
            ax,
            (x, y),
            f"{tf}\nPLAUR r{int(float(row['rank']))}; MES-V r{int(float(row['MESV_RSS_rank']))}",
            color,
            edge,
            w=1.65,
            h=0.46,
            fs=7.5,
            weight="bold" if tf == "FOSL1" else "normal",
        )

        add_arrow(ax, (2.25, 5.15), (x - 0.85, y + 0.06), color="#526D9E", lw=0.75, alpha=0.55, rad=0.08)
        add_arrow(ax, (2.25, 3.05), (x - 0.85, y - 0.06), color="#C85B4B", lw=0.75, alpha=0.55, rad=-0.05)

    # FOSL1 uniquely passes ATAC filter.
    fpos = tf_pos["FOSL1"]
    add_arrow(ax, (fpos[0] + 0.92, fpos[1]), (6.35, 3.05), color="#D58B2A", lw=1.8, alpha=0.95)
    add_arrow(ax, (7.3, 3.4), (8.5, 5.0), color="#9C2F2F", lw=1.4, alpha=0.9)
    add_arrow(ax, (7.6, 2.75), (8.45, 1.58), color="#D58B2A", lw=1.2, alpha=0.9, ls="--")

    ax.text(5.7, 3.52, "only FOSL1\npasses ATAC filter", ha="center", va="center", fontsize=7.2, color="#9C2F2F")
    ax.text(5.0, 5.9, "SCENIC + ATAC evidence convergence prioritizes FOSL1", ha="center", fontsize=10.5, weight="bold")
    ax.text(5.0, 5.62, "Edges indicate evidence links for prioritization, not causal regulation.", ha="center", fontsize=7.3, color="#555555")

    legend_elements = [
        Line2D([0], [0], color="#526D9E", lw=1.4, label="SCENIC PLAUR edge"),
        Line2D([0], [0], color="#C85B4B", lw=1.4, label="MES-V RSS context"),
        Line2D([0], [0], color="#D58B2A", lw=1.4, label="ATAC motif/anchor context"),
    ]
    ax.legend(handles=legend_elements, loc="lower left", bbox_to_anchor=(0.02, 0.05), frameon=False, fontsize=7.2)
    fig.tight_layout()
    for out_dir in [candidate_dir, main_dir]:
        fig.savefig(out_dir / "E_FOSL1证据收敛网络.pdf", bbox_inches="tight")
    fig.savefig(candidate_dir / "E_FOSL1证据收敛网络.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    manifest = pd.DataFrame(
        [
            {
                "panel": "D",
                "figure": "D_候选TF_target功能语义点图.pdf",
                "role": "Functional context of targets for the six TFs passing SCENIC PLAUR-upstream and MES-V RSS filters",
                "claim_boundary": "Exploratory curated target-gene context; not formal GO and not causal FOSL1-PLAUR evidence",
                "source_data": "source_data/DE_SCENIC_target_context_network/D_candidate_TF_target_functional_context.csv",
                "script": script_name,
            },
            {
                "panel": "E",
                "figure": "E_FOSL1证据收敛网络.pdf",
                "role": "Prioritization network showing SCENIC/MES-V RSS/ATAC evidence convergence to FOSL1",
                "claim_boundary": "Evidence links for prioritization only; not a causal GRN",
                "source_data": "source_data/DE_SCENIC_target_context_network/E_prioritization_network_nodes.csv",
                "script": script_name,
            },
        ]
    )
    manifest.to_csv(record_dir / "PanelDE_SCENIC_target_context_network_manifest.csv", index=False)


if __name__ == "__main__":
    main()
