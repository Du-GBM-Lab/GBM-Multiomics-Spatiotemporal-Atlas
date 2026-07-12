#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


SUBTYPE_ORDER = ["NPC-P", "OPC-M", "MES-V", "MES-I"]
SUBTYPE_COLORS = {
    "NPC-P": "#0072B5",
    "OPC-M": "#E18727",
    "MES-V": "#20854E",
    "MES-I": "#BC3C29",
}
TOP_N = 20
AP1_FOCUS_TFS = {"FOSL1", "FOSL2", "ATF3", "JUN", "JUNB", "FOSB"}
PANEL_C_TFS = {
    "NFATC4", "CREB3L2", "HIC1", "ETS1", "NR2F2", "BHLHE40",
    "FOSL2", "FOSL1", "ATF3", "JUN", "JUNB", "FOSB",
}


def main() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    outputs_dir = project_dir / "outputs"
    figures_dir = project_dir / "figures"
    fig_dir = figures_dir / "正文图"
    source_dir = figures_dir / "source_data" / "B_四亚型top_regulon排序"
    record_dir = figures_dir / "记录"
    for path in [fig_dir, source_dir, record_dir]:
        path.mkdir(parents=True, exist_ok=True)

    rank_tbl = pd.read_csv(outputs_dir / "scenic_rss_all_regulons_ranked_by_subtype.csv")
    plot_tbl = (
        rank_tbl.loc[rank_tbl["subtype"].isin(SUBTYPE_ORDER) & rank_tbl["rank"].le(TOP_N)]
        .copy()
        .sort_values(["subtype", "rank"])
    )
    plot_tbl["subtype"] = pd.Categorical(plot_tbl["subtype"], SUBTYPE_ORDER, ordered=True)
    plot_tbl["label"] = plot_tbl["rank"].astype(str) + ". " + plot_tbl["tf"].astype(str)

    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "sans-serif"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    fig, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(7.0, 7.6),
        sharex=True,
        sharey=False,
        constrained_layout=False,
    )
    axes = axes.ravel()
    xmax = plot_tbl["rss"].max() * 1.08
    for ax, subtype in zip(axes, SUBTYPE_ORDER):
        one = plot_tbl.loc[plot_tbl["subtype"].eq(subtype)].sort_values("rank")
        y = list(range(len(one)))
        color = SUBTYPE_COLORS[subtype]
        label_weight = ["bold" if tf in PANEL_C_TFS else "normal" for tf in one["tf"]]
        ax.hlines(y=y, xmin=0, xmax=one["rss"], color="#D8D8D8", linewidth=0.8)
        ax.scatter(one["rss"], y, s=23, color=color, edgecolor="white", linewidth=0.35, zorder=3)

        focus = one["tf"].isin(AP1_FOCUS_TFS)
        if focus.any():
            ax.scatter(
                one.loc[focus, "rss"],
                [y[i] for i, keep in enumerate(focus.tolist()) if keep],
                s=58,
                facecolors="none",
                edgecolors="#111111",
                linewidth=0.8,
                zorder=4,
            )

        ax.set_yticks(y)
        ax.set_yticklabels(one["label"], fontsize=6.8)
        for tick, weight in zip(ax.get_yticklabels(), label_weight):
            tick.set_fontweight(weight)
        ax.invert_yaxis()
        ax.set_xlim(0, xmax)
        ax.set_title(subtype, color=color, fontsize=10.0, fontweight="bold", pad=4)
        ax.grid(axis="x", color="#ECECEC", linewidth=0.45)
        ax.tick_params(axis="x", labelsize=7)
        ax.tick_params(axis="y", length=0, pad=1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_color("#A0A0A0")

    for ax in axes:
        ax.set_xlabel("RSS", fontsize=8)
    fig.suptitle("Top RSS-ranked regulons by malignant subtype", fontsize=11.0, fontweight="bold")
    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="#555555",
                   markeredgecolor="white", markersize=5, label="RSS-ranked regulon"),
        plt.Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="none",
                   markeredgecolor="#111111", markersize=7, label="AP-1/FOS family"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.018),
               ncol=2, frameon=False, fontsize=8)
    fig.subplots_adjust(left=0.10, right=0.985, top=0.91, bottom=0.08, wspace=0.30, hspace=0.25)

    pdf_path = fig_dir / "B_四亚型top_regulon排序.pdf"
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    plt.close(fig)

    plot_tbl.to_csv(source_dir / "top20_regulons_by_subtype_RSS.csv", index=False)
    rank_tbl.to_csv(source_dir / "all_regulon_RSS_ranking.csv", index=False)

    manifest = pd.DataFrame(
        [
            {
                "figure": "Panel B 四亚型top regulon排序",
                "pdf": str(pdf_path),
                "source_data_dir": str(source_dir),
                "selection_rule": f"top {TOP_N} RSS-ranked regulons per subtype",
                "color_quantity": "RSS",
                "n_rows": int(plot_tbl.shape[0]),
                "FOSL1_MESV_rank": int(
                    rank_tbl.loc[
                        rank_tbl["subtype"].eq("MES-V") & rank_tbl["tf"].eq("FOSL1"),
                        "rank",
                    ].iloc[0]
                ),
                "FOSL1_MESV_RSS": float(
                    rank_tbl.loc[
                        rank_tbl["subtype"].eq("MES-V") & rank_tbl["tf"].eq("FOSL1"),
                        "rss",
                    ].iloc[0]
                ),
                "subtype_order": ",".join(SUBTYPE_ORDER),
            }
        ]
    )
    manifest.to_csv(record_dir / "PanelB_四亚型top_regulon排序_manifest.csv", index=False)
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
