#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SUBTYPE_ORDER = ["NPC-P", "OPC-M", "MES-V", "MES-I"]
PANEL_TFS = [
    "NFATC4", "CREB3L2", "HIC1", "ETS1", "NR2F2", "BHLHE40",
    "FOSL2", "FOSL1", "ATF3", "JUN", "JUNB", "FOSB",
]
AP1_FAMILY_TFS = {"FOSL2", "FOSL1", "ATF3", "JUN", "JUNB", "FOSB"}
PANEL_RATIONALE = (
    "MES-V RSS top regulators plus retained AP-1/FOS family regulons"
)
SUBTYPE_COLORS = {
    "NPC-P": "#0072B5",
    "OPC-M": "#E18727",
    "MES-V": "#20854E",
    "MES-I": "#BC3C29",
}


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    values = df.to_numpy(dtype=float)
    mean = values.mean(axis=1, keepdims=True)
    sd = values.std(axis=1, ddof=0, keepdims=True)
    sd[sd == 0] = 1.0
    return pd.DataFrame((values - mean) / sd, index=df.index, columns=df.columns)


def mean_auc_by_subtype(auc: ad.AnnData, regulons: list[str]) -> pd.DataFrame:
    x = auc[:, regulons].X
    if not isinstance(x, np.ndarray):
        x = x.toarray()
    auc_df = pd.DataFrame(x, index=auc.obs_names, columns=regulons)
    subtype = auc.obs["cell_type"].astype(str)
    means = auc_df.groupby(subtype).mean().T
    return means.loc[regulons, SUBTYPE_ORDER]


def main() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    outputs_dir = project_dir / "outputs"
    figures_dir = project_dir / "figures"
    fig_dir = figures_dir / "正文图"
    source_dir = figures_dir / "source_data" / "C_MES_regulator_panel"
    record_dir = figures_dir / "记录"
    for path in [fig_dir, source_dir, record_dir]:
        path.mkdir(parents=True, exist_ok=True)

    rank_tbl = pd.read_csv(outputs_dir / "scenic_rss_all_regulons_ranked_by_subtype.csv")
    auc = ad.read_h5ad(outputs_dir / "scenic_AUCell_regulon_activity.h5ad")
    retained_tfs = set(rank_tbl["tf"].astype(str))
    tf_to_regulon = (
        rank_tbl.drop_duplicates("tf")
        .set_index("tf")["regulon"]
        .to_dict()
    )
    retained_panel = [tf for tf in PANEL_TFS if tf in retained_tfs]
    retained_regulons = [tf_to_regulon[tf] for tf in retained_panel]

    mean_auc = mean_auc_by_subtype(auc, retained_regulons)
    mean_auc.index = retained_panel
    z_auc = zscore_rows(mean_auc)

    rows = []
    for tf in PANEL_TFS:
        retained = tf in retained_tfs
        for subtype in SUBTYPE_ORDER:
            row = {
                "tf": tf,
                "subtype": subtype,
                "retained_regulon": retained,
                "regulon": tf_to_regulon.get(tf, pd.NA),
                "rss": pd.NA,
                "rss_rank": pd.NA,
                "mean_auc": pd.NA,
                "mean_auc_z": pd.NA,
            }
            if retained:
                hit = rank_tbl.loc[
                    rank_tbl["tf"].eq(tf) & rank_tbl["subtype"].eq(subtype)
                ].iloc[0]
                row.update(
                    {
                        "rss": float(hit["rss"]),
                        "rss_rank": int(hit["rank"]),
                        "mean_auc": float(mean_auc.loc[tf, subtype]),
                        "mean_auc_z": float(z_auc.loc[tf, subtype]),
                    }
                )
            rows.append(row)
    plot_tbl = pd.DataFrame(rows)

    retained_summary = (
        plot_tbl.drop_duplicates("tf")
        .loc[:, ["tf", "retained_regulon", "regulon"]]
        .copy()
    )
    retained_summary["note"] = np.where(
        retained_summary["retained_regulon"],
        "retained in SCENIC regulons",
        "not retained in SCENIC regulons; not manually added",
    )

    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "sans-serif"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    fig, ax = plt.subplots(figsize=(5.10, 5.45))
    x_pos = {s: i for i, s in enumerate(SUBTYPE_ORDER)}
    y_order = PANEL_TFS[::-1]
    y_pos = {tf: i for i, tf in enumerate(y_order)}

    retained_rows = plot_tbl.loc[plot_tbl["retained_regulon"]].copy()
    rss = retained_rows["rss"].astype(float)
    size = 60 + 520 * (rss - rss.min()) / (rss.max() - rss.min())
    sc = ax.scatter(
        retained_rows["subtype"].map(x_pos),
        retained_rows["tf"].map(y_pos),
        s=size,
        c=retained_rows["mean_auc_z"],
        cmap="RdBu_r",
        vmin=-1.5,
        vmax=1.5,
        edgecolor="black",
        linewidth=0.35,
        alpha=0.95,
    )

    missing_rows = plot_tbl.loc[~plot_tbl["retained_regulon"]].copy()
    if not missing_rows.empty:
        ax.scatter(
            missing_rows["subtype"].map(x_pos),
            missing_rows["tf"].map(y_pos),
            marker="x",
            s=34,
            color="#9A9A9A",
            linewidth=0.9,
        )

    ax.set_xticks(range(len(SUBTYPE_ORDER)))
    ax.set_xticklabels(SUBTYPE_ORDER, rotation=30, ha="right")
    ax.set_yticks(range(len(y_order)))
    ax.set_yticklabels(y_order)
    for label in ax.get_yticklabels():
        if label.get_text() in AP1_FAMILY_TFS:
            label.set_fontweight("bold")

    ax.set_xlim(-0.55, len(SUBTYPE_ORDER) - 0.45)
    ax.set_ylim(-0.55, len(y_order) - 0.45)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("MES-V top and AP-1/FOS regulon activity", fontsize=10.5, fontweight="bold")
    ax.axhline(5.5, color="#CFCFCF", linewidth=0.8, linestyle="--", zorder=0)
    ax.grid(color="#E8E8E8", linewidth=0.7)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)

    cbar = fig.colorbar(sc, ax=ax, fraction=0.05, pad=0.02)
    cbar.set_label("Mean AUCell z-score", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    size_values = [0.25, 0.40, 0.55]
    handles = [
        plt.scatter([], [], s=60 + 520 * (v - rss.min()) / (rss.max() - rss.min()),
                    facecolor="white", edgecolor="black", linewidth=0.35)
        for v in size_values
    ]
    labels = [f"RSS {v:.2f}" for v in size_values]
    if not missing_rows.empty:
        handles.append(plt.scatter([], [], marker="x", s=34, color="#9A9A9A", linewidth=0.9))
        labels.append("not retained")
    ax.legend(
        handles,
        labels,
        title="RSS",
        loc="lower center",
        bbox_to_anchor=(0.50, -0.24),
        ncol=3,
        frameon=False,
        fontsize=7,
        title_fontsize=8,
        handletextpad=0.4,
        columnspacing=0.9,
    )

    fig.subplots_adjust(left=0.16, right=0.86, top=0.90, bottom=0.22)
    pdf_path = fig_dir / "C_MES_regulator_panel点图.pdf"
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
    plt.close(fig)

    plot_tbl.to_csv(source_dir / "MES_regulator_panel_RSS_AUCell.csv", index=False)
    retained_summary.to_csv(source_dir / "MES_regulator_retained_summary.csv", index=False)

    manifest = pd.DataFrame(
        [
            {
                "figure": "Panel C MES-V top and AP-1/FOS regulon activity",
                "pdf": str(pdf_path),
                "source_data_dir": str(source_dir),
                "panel_tfs": ",".join(PANEL_TFS),
                "panel_rationale": PANEL_RATIONALE,
                "retained_tfs": ",".join(retained_panel),
                "missing_tfs": ",".join([tf for tf in PANEL_TFS if tf not in retained_tfs]),
                "dot_color": "mean AUCell z-score by subtype",
                "dot_size": "RSS",
                "FOSL1_MESV_rank": int(
                    rank_tbl.loc[
                        rank_tbl["tf"].eq("FOSL1") & rank_tbl["subtype"].eq("MES-V"),
                        "rank",
                    ].iloc[0]
                ),
                "FOSL1_MESV_RSS": float(
                    rank_tbl.loc[
                        rank_tbl["tf"].eq("FOSL1") & rank_tbl["subtype"].eq("MES-V"),
                        "rss",
                    ].iloc[0]
                ),
            }
        ]
    )
    manifest.to_csv(record_dir / "PanelC_MES_regulator_panel_manifest.csv", index=False)
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
