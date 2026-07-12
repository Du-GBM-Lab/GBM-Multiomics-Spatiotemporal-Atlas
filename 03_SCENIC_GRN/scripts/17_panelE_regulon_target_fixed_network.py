#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


SUBTYPE_TFS = {
    "NPC-P": ["E2F1", "E2F7", "TFDP1"],
    "OPC-M": ["SOX10", "SOX8", "SREBF1"],
    "MES-V": ["FOSL1", "ETS1", "BHLHE40"],
    "MES-I": ["SPI1", "IRF8", "IKZF1"],
}

COLORS = {
    "NPC-P": "#4C78A8",
    "OPC-M": "#4B9A8D",
    "MES-V": "#C85B4B",
    "MES-I": "#8F6BB1",
}

FOCUS_GENES = ["PLAUR", "PLAU", "VIM", "CD44", "COL1A1", "MMP9", "FN1", "GAP43", "TOP2A"]
FORCE_TARGETS = {"FOSL1": ["PLAUR", "PLAU", "VIM", "CD44", "FN1"]}


def add_confidence_ellipse(ax, xy: np.ndarray, color: str, label: str) -> None:
    if xy.shape[0] < 3:
        return
    center = xy.mean(axis=0)
    cov = np.cov(xy.T)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    angle = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    width, height = 2.4 * np.sqrt(np.maximum(vals, 1e-4))
    ell = patches.Ellipse(
        center,
        width=width,
        height=height,
        angle=angle,
        facecolor=color,
        edgecolor=color,
        lw=1.2,
        alpha=0.11,
        linestyle="--",
        zorder=0,
    )
    ax.add_patch(ell)
    ax.text(center[0], center[1] + height * 0.52, label, color=color, fontsize=9, weight="bold", ha="center")


def main() -> None:
    scenic_dir = Path(__file__).resolve().parents[1]
    project_root = scenic_dir.parents[2]
    final_root = project_root / "图片表格" / "05_发育时间_TF_ATAC验证" / "02_TF_regulon_SCENIC"

    candidate_dir = final_root / "正文图候选"
    main_dir = final_root / "正文图"
    source_dir = final_root / "source_data" / "E_四亚型代表regulon_target_network"
    record_dir = final_root / "记录"
    for d in [candidate_dir, main_dir, source_dir, record_dir]:
        d.mkdir(parents=True, exist_ok=True)

    grn = pd.read_csv(scenic_dir / "outputs" / "scenic_grn_edgelist.csv")
    grn["TF"] = grn["TF"].astype(str)
    grn["target"] = grn["target"].astype(str)
    grn["importance"] = pd.to_numeric(grn["importance"], errors="coerce").fillna(0)

    tf_to_subtype = {tf: subtype for subtype, tfs in SUBTYPE_TFS.items() for tf in tfs}
    available_tfs = [tf for tf in tf_to_subtype if tf in set(grn["TF"])]

    G = nx.Graph()
    edge_rows = []
    node_meta = {}

    top_n = 12
    for tf in available_tfs:
        subtype = tf_to_subtype[tf]
        G.add_node(tf)
        node_meta[tf] = {"node_type": "TF", "subtype": subtype, "highlight": tf == "FOSL1"}

        rel = grn.loc[grn["TF"].eq(tf)].sort_values("importance", ascending=False)
        targets = rel.head(top_n)["target"].tolist()
        for forced in FORCE_TARGETS.get(tf, []):
            if forced in set(rel["target"]):
                targets.append(forced)
        targets = list(dict.fromkeys(targets))

        for target in targets:
            imp = float(rel.loc[rel["target"].eq(target), "importance"].iloc[0])
            G.add_node(target)
            if target not in node_meta:
                node_meta[target] = {
                    "node_type": "target",
                    "subtype": subtype,
                    "highlight": target in FOCUS_GENES,
                }
            elif node_meta[target]["node_type"] == "target":
                node_meta[target]["subtype"] = "shared"
                node_meta[target]["highlight"] = node_meta[target]["highlight"] or target in FOCUS_GENES
            G.add_edge(tf, target, weight=max(imp, 1e-4))
            edge_rows.append({"TF": tf, "target": target, "subtype": subtype, "importance": imp})

    # One reproducible force-directed coordinate system. This is the corrected Fig.3E logic:
    # all TFs and targets occupy the same layout; subtype regions are overlays, not separate panels.
    pos = nx.spring_layout(G, k=1.35, iterations=600, seed=42, weight="weight", scale=4.0)

    node_rows = []
    for node, (x, y) in pos.items():
        meta = node_meta.get(node, {"node_type": "target", "subtype": "unknown", "highlight": False})
        node_rows.append({"node": node, "x": x, "y": y, **meta})
    nodes_df = pd.DataFrame(node_rows)
    edges_df = pd.DataFrame(edge_rows)
    nodes_df.to_csv(source_dir / "E_fixed_network_nodes.csv", index=False)
    edges_df.to_csv(source_dir / "E_fixed_network_edges.csv", index=False)
    pd.DataFrame(
        [{"subtype": s, "TF": tf} for s, tfs in SUBTYPE_TFS.items() for tf in tfs if tf in available_tfs]
    ).to_csv(source_dir / "E_representative_TF_selection.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.4, 5.65))
    ax.axis("off")
    ax.set_title("Subtype representative TF-target regulon landscape", fontsize=11.2, weight="bold", pad=7)

    # Draw subtype regions from TF + their immediate targets in the shared coordinate system.
    for subtype, color in COLORS.items():
        sub_nodes = []
        for tf in [t for t in SUBTYPE_TFS[subtype] if t in G.nodes]:
            sub_nodes.append(tf)
            sub_nodes.extend(list(G.neighbors(tf)))
        sub_nodes = [n for n in dict.fromkeys(sub_nodes) if n in pos]
        xy = np.array([pos[n] for n in sub_nodes])
        add_confidence_ellipse(ax, xy, color, subtype)

    weights = [d["weight"] for _, _, d in G.edges(data=True)]
    w_min, w_max = min(weights), max(weights)
    for u, v, d in G.edges(data=True):
        norm = (d["weight"] - w_min) / (w_max - w_min) if w_max > w_min else 0.5
        color = COLORS.get(tf_to_subtype.get(u, tf_to_subtype.get(v, "")), "#B0B0B0")
        ax.plot(
            [pos[u][0], pos[v][0]],
            [pos[u][1], pos[v][1]],
            color=color,
            alpha=0.13 + 0.22 * norm,
            lw=0.35 + 1.1 * norm,
            zorder=1,
        )

    target_nodes = [n for n, m in node_meta.items() if m["node_type"] == "target" and not m["highlight"]]
    focus_nodes = [n for n, m in node_meta.items() if m["node_type"] == "target" and m["highlight"]]
    tf_nodes = [n for n, m in node_meta.items() if m["node_type"] == "TF"]

    for subtype, color in COLORS.items():
        sub_targets = [n for n in target_nodes if node_meta[n]["subtype"] == subtype]
        if sub_targets:
            ax.scatter(
                [pos[n][0] for n in sub_targets],
                [pos[n][1] for n in sub_targets],
                s=35,
                c=color,
                alpha=0.42,
                edgecolors="white",
                linewidths=0.25,
                zorder=2,
            )

    shared_targets = [n for n in target_nodes if node_meta[n]["subtype"] == "shared"]
    if shared_targets:
        ax.scatter(
            [pos[n][0] for n in shared_targets],
            [pos[n][1] for n in shared_targets],
            s=38,
            c="#9AA0A6",
            alpha=0.48,
            edgecolors="white",
            linewidths=0.25,
            zorder=2,
        )

    if focus_nodes:
        ax.scatter(
            [pos[n][0] for n in focus_nodes],
            [pos[n][1] for n in focus_nodes],
            s=86,
            c="#FFD44D",
            alpha=0.96,
            edgecolors="#6A5200",
            linewidths=0.75,
            zorder=4,
        )

    for subtype, color in COLORS.items():
        sub_tfs = [n for n in tf_nodes if node_meta[n]["subtype"] == subtype]
        if sub_tfs:
            ax.scatter(
                [pos[n][0] for n in sub_tfs],
                [pos[n][1] for n in sub_tfs],
                s=[235 if n == "FOSL1" else 170 for n in sub_tfs],
                c=color,
                alpha=0.98,
                edgecolors="#262626",
                linewidths=0.8,
                zorder=5,
            )

    # Labels: all TFs, plus selected focus targets. Ordinary targets remain unlabeled to keep the panel readable.
    for n in tf_nodes:
        x, y = pos[n]
        ax.text(x, y, n, ha="center", va="center", fontsize=6.4, weight="bold", color="white", zorder=6)
    for n in focus_nodes:
        x, y = pos[n]
        ax.text(x, y + 0.17, n, ha="center", va="bottom", fontsize=6.2, weight="bold", color="#3B3100", zorder=6)

    handles = [
        mlines.Line2D([], [], color=color, marker="o", linestyle="None", markersize=6.5, label=subtype)
        for subtype, color in COLORS.items()
    ]
    handles.append(mlines.Line2D([], [], color="#FFD44D", marker="o", markeredgecolor="#6A5200", linestyle="None", markersize=7, label="highlighted target"))
    handles.append(mlines.Line2D([], [], color="#9AA0A6", marker="o", linestyle="None", markersize=6, label="shared target"))
    ax.legend(handles=handles, loc="lower right", fontsize=6.8, frameon=False, borderpad=0.2, handletextpad=0.4)

    ax.text(
        0.01,
        0.015,
        "Shared force-directed layout (seed=42); edges indicate SCENIC-inferred TF-target links.",
        transform=ax.transAxes,
        fontsize=6.4,
        color="#555555",
    )

    fig.tight_layout()
    for out_dir in [candidate_dir, main_dir]:
        fig.savefig(out_dir / "E_四亚型代表regulon固定网络.pdf", bbox_inches="tight")
    fig.savefig(candidate_dir / "E_四亚型代表regulon固定网络.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    manifest = pd.DataFrame(
        [
            {
                "panel": "E",
                "figure": "E_四亚型代表regulon固定网络.pdf",
                "role": "Shared force-directed TF-target network for multiple representative regulons per malignant subtype",
                "layout": "networkx.spring_layout(k=1.35, iterations=600, seed=42, weight='weight', scale=4.0)",
                "target_definition": f"Top {top_n} SCENIC GRN targets per representative TF, with selected FOSL1 targets forced if present",
                "claim_boundary": "Shows subtype-region regulon context and nominated target links; not wet-lab validated causality",
                "source_data": "source_data/E_四亚型代表regulon_target_network/E_fixed_network_nodes.csv",
            }
        ]
    )
    manifest.to_csv(record_dir / "PanelE_四亚型代表regulon固定网络_manifest.csv", index=False)


if __name__ == "__main__":
    main()
