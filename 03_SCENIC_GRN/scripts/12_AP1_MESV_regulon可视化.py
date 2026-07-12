#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


SUBTYPE_ORDER = ["NPC-P", "OPC-M", "MES-V", "MES-I"]
AP1_PREFIXES = ("FOS", "JUN", "ATF", "BATF")
MANDATORY_TFS = ["FOSL1", "FOSL2", "JUNB"]


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
    fig_dir = figures_dir / "补充图"
    source_dir = figures_dir / "source_data" / "AP1_MESV_regulon"
    record_dir = figures_dir / "记录"
    for path in [fig_dir, source_dir, record_dir]:
        path.mkdir(parents=True, exist_ok=True)

    rank_tbl = pd.read_csv(outputs_dir / "scenic_rss_all_regulons_ranked_by_subtype.csv")
    edgelist = pd.read_csv(outputs_dir / "scenic_grn_edgelist.csv")
    auc = ad.read_h5ad(outputs_dir / "scenic_AUCell_regulon_activity.h5ad")

    retained_tfs = sorted(rank_tbl["tf"].astype(str).unique())
    ap1_tfs = [tf for tf in retained_tfs if tf.startswith(AP1_PREFIXES)]
    missing_mandatory = [tf for tf in MANDATORY_TFS if tf not in ap1_tfs]
    if missing_mandatory:
        raise ValueError(f"Mandatory AP-1 TFs missing from retained regulons: {missing_mandatory}")

    n_targets = (
        edgelist.loc[edgelist["TF"].astype(str).isin(ap1_tfs)]
        .groupby("TF")["target"]
        .nunique()
        .rename("n_targets_from_grn_edgelist")
    )
    mesv = rank_tbl.loc[rank_tbl["subtype"].eq("MES-V")].set_index("tf")
    ap1_summary = (
        pd.DataFrame({"tf": ap1_tfs})
        .assign(regulon=lambda d: d["tf"] + "(+)")
        .merge(n_targets.reset_index(), left_on="tf", right_on="TF", how="left")
        .drop(columns=["TF"])
    )
    ap1_summary["n_targets_from_grn_edgelist"] = (
        ap1_summary["n_targets_from_grn_edgelist"].fillna(0).astype(int)
    )
    ap1_summary["MESV_RSS_rank"] = ap1_summary["tf"].map(mesv["rank"]).astype(int)
    ap1_summary["MESV_RSS"] = ap1_summary["tf"].map(mesv["rss"]).astype(float)
    ap1_summary = ap1_summary.sort_values(["MESV_RSS_rank", "tf"]).reset_index(drop=True)

    regulons = ap1_summary["regulon"].tolist()
    mean_auc = mean_auc_by_subtype(auc, regulons)
    mean_auc.index = ap1_summary["tf"].tolist()
    z_auc = zscore_rows(mean_auc)

    # Plot 1: AP-1 activity across subtypes.
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "sans-serif"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    fig, ax = plt.subplots(figsize=(4.3, 3.2))
    sns.heatmap(
        z_auc,
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-1.5,
        vmax=1.5,
        linewidths=0.35,
        linecolor="white",
        cbar_kws={"label": "Mean AUCell z-score"},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("AP-1 family regulon activity", fontsize=10.5, fontweight="bold")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=35, ha="right", fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
    fig.subplots_adjust(left=0.20, right=0.95, top=0.86, bottom=0.20)
    activity_pdf = fig_dir / "S_AP1_regulon跨亚型活性热图.pdf"
    fig.savefig(activity_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)

    # Plot 2: AP-1 members on the MES-V RSS ranking background.
    all_mesv = rank_tbl.loc[rank_tbl["subtype"].eq("MES-V")].sort_values("rank")
    fig, ax = plt.subplots(figsize=(5.7, 3.3))
    ax.plot(all_mesv["rank"], all_mesv["rss"], color="#C7C7C7", linewidth=1.0, zorder=1)
    ax.scatter(all_mesv["rank"], all_mesv["rss"], s=10, color="#C7C7C7", zorder=2)
    focus = all_mesv.loc[all_mesv["tf"].isin(ap1_tfs)].copy()
    ax.scatter(
        focus["rank"],
        focus["rss"],
        s=46,
        color="#BC3C29",
        edgecolor="black",
        linewidth=0.35,
        zorder=3,
    )
    for _, row in focus.iterrows():
        ax.text(
            row["rank"] + 1.0,
            row["rss"] + 0.004,
            f"{row['tf']} ({int(row['rank'])})",
            fontsize=7,
            ha="left",
            va="center",
        )
    ax.invert_xaxis()
    ax.set_xlabel("MES-V RSS rank (left = higher rank)", fontsize=8)
    ax.set_ylabel("RSS", fontsize=8)
    ax.set_title("AP-1 members in MES-V RSS ranking", fontsize=10.5, fontweight="bold")
    ax.grid(color="#E8E8E8", linewidth=0.6)
    ax.tick_params(labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.subplots_adjust(left=0.13, right=0.98, top=0.86, bottom=0.18)
    rank_pdf = fig_dir / "S_AP1_MESV_RSS排序标注.pdf"
    fig.savefig(rank_pdf, format="pdf", bbox_inches="tight")
    plt.close(fig)

    ap1_summary.to_csv(source_dir / "retained_AP1_regulon清单.csv", index=False)
    mean_auc.to_csv(source_dir / "AP1_mean_AUCell_by_subtype.csv")
    z_auc.to_csv(source_dir / "AP1_mean_AUCell_zscore_by_subtype.csv")
    all_mesv.to_csv(source_dir / "MESV_all_regulon_RSS_ranking.csv", index=False)
    focus.to_csv(source_dir / "MESV_AP1_regulon_RSS_ranking.csv", index=False)

    manifest = pd.DataFrame(
        [
            {
                "figure": "S_AP1_regulon跨亚型活性热图",
                "pdf": str(activity_pdf),
                "source_data_dir": str(source_dir),
                "included_tfs": ",".join(ap1_summary["tf"]),
                "selection_rule": "retained regulons with TF prefix FOS/JUN/ATF/BATF",
                "color_quantity": "subtype-level mean AUCell z-score",
                "interpretation": "AP-1 family regulons are associated with MES-V regulon program; nomination only",
            },
            {
                "figure": "S_AP1_MESV_RSS排序标注",
                "pdf": str(rank_pdf),
                "source_data_dir": str(source_dir),
                "included_tfs": ",".join(ap1_summary["tf"]),
                "selection_rule": "same AP-1 retained regulons highlighted on MES-V RSS ranking",
                "color_quantity": "RSS",
                "interpretation": "FOSL2 and FOSL1 ranks are shown as observed; no implication that FOSL1 is stronger than FOSL2",
            },
        ]
    )
    manifest.to_csv(record_dir / "AP1_MESV_regulon_manifest.csv", index=False)

    print(ap1_summary.to_string(index=False))
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
