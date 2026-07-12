#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path
import os

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap


SUBTYPE_ORDER = ["NPC-P", "OPC-M", "MES-V", "MES-I"]
SUBTYPE_COLORS = {
    "NPC-P": "#0072B5",
    "OPC-M": "#E18727",
    "MES-V": "#20854E",
    "MES-I": "#BC3C29",
}
TOP_N_PER_SUBTYPE = 7


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    values = df.to_numpy(dtype=float)
    mean = values.mean(axis=1, keepdims=True)
    sd = values.std(axis=1, ddof=0, keepdims=True)
    sd[sd == 0] = 1.0
    return pd.DataFrame((values - mean) / sd, index=df.index, columns=df.columns)


def select_regulons(rank_tbl: pd.DataFrame) -> pd.DataFrame:
    selected = rank_tbl.loc[rank_tbl["subtype"].isin(SUBTYPE_ORDER)].copy()
    selected = selected.loc[selected["rank"].le(TOP_N_PER_SUBTYPE)].copy()

    all_rss = rank_tbl.pivot(index="regulon", columns="subtype", values="rss")
    rows = []
    for regulon, one in selected.groupby("regulon"):
        rss_row = all_rss.loc[regulon, SUBTYPE_ORDER]
        selected_subtypes = one["subtype"].tolist()
        owner_subtype = rss_row.loc[selected_subtypes].astype(float).idxmax()
        owner_rank = int(
            rank_tbl.loc[
                rank_tbl["subtype"].eq(owner_subtype) & rank_tbl["regulon"].eq(regulon),
                "rank",
            ].iloc[0]
        )
        rows.append(
            {
                "regulon": regulon,
                "tf": str(one["tf"].iloc[0]),
                "assigned_subtype": owner_subtype,
                "assigned_subtype_rank": owner_rank,
                "assigned_subtype_rss": float(rss_row[owner_subtype]),
                "selected_from_subtypes": ";".join(one["subtype"].tolist()),
            }
        )

    out = pd.DataFrame(rows)
    out["assigned_subtype"] = pd.Categorical(
        out["assigned_subtype"], categories=SUBTYPE_ORDER, ordered=True
    )
    out = out.sort_values(
        ["assigned_subtype", "assigned_subtype_rank", "assigned_subtype_rss", "regulon"],
        ascending=[True, True, False, True],
    ).reset_index(drop=True)
    return out


def mean_auc_by_subtype(auc: ad.AnnData, regulons: list[str]) -> pd.DataFrame:
    if "cell_type" not in auc.obs.columns:
        raise ValueError("AUCell h5ad missing obs['cell_type']")
    missing = [x for x in regulons if x not in auc.var_names]
    if missing:
        raise ValueError(f"Selected regulons missing from AUCell matrix: {missing}")

    x = auc[:, regulons].X
    if not isinstance(x, np.ndarray):
        x = x.toarray()
    auc_df = pd.DataFrame(x, index=auc.obs_names, columns=regulons)
    subtype = auc.obs["cell_type"].astype(str)
    means = auc_df.groupby(subtype).mean().T
    return means.loc[regulons, SUBTYPE_ORDER]


def save_heatmap(z: pd.DataFrame, selected: pd.DataFrame, pdf_path: Path) -> None:
    sns.set_theme(style="white", font="Arial", font_scale=0.9)
    display = z.copy()
    display.index = selected["tf"].tolist()

    fig_height = max(6.2, 0.23 * display.shape[0] + 1.2)
    fig = plt.figure(figsize=(7.8, fig_height))
    gs = fig.add_gridspec(
        nrows=2,
        ncols=4,
        width_ratios=[0.16, 4.25, 0.9, 1.25],
        height_ratios=[4.8, 0.22],
        wspace=0.04,
        hspace=0.32,
    )
    ax_group = fig.add_subplot(gs[0, 0])
    ax_heat = fig.add_subplot(gs[0, 1])
    ax_spacer = fig.add_subplot(gs[0, 2])
    ax_legend = fig.add_subplot(gs[0, 3])
    ax_cbar = fig.add_subplot(gs[1, 1])

    group_codes = selected["assigned_subtype"].map({x: i for i, x in enumerate(SUBTYPE_ORDER)})
    ax_group.imshow(
        group_codes.to_numpy()[:, None],
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap([SUBTYPE_COLORS[x] for x in SUBTYPE_ORDER]),
        vmin=0,
        vmax=len(SUBTYPE_ORDER) - 1,
    )
    ax_group.set_xticks([])
    ax_group.set_yticks([])
    for spine in ax_group.spines.values():
        spine.set_visible(False)

    hm = sns.heatmap(
        display,
        ax=ax_heat,
        cmap="RdBu_r",
        center=0,
        vmin=-1.8,
        vmax=1.8,
        linewidths=0.25,
        linecolor="white",
        cbar=True,
        cbar_ax=ax_cbar,
        cbar_kws={"label": "Mean AUCell activity z-score", "orientation": "horizontal"},
    )
    ax_heat.set_xlabel("")
    ax_heat.set_ylabel("")
    ax_heat.yaxis.tick_right()
    ax_heat.tick_params(axis="both", length=0)
    ax_heat.set_xticklabels(ax_heat.get_xticklabels(), rotation=35, ha="right")
    ax_heat.set_yticklabels(ax_heat.get_yticklabels(), rotation=0)
    ax_heat.tick_params(axis="y", pad=3)

    ax_cbar.tick_params(labelsize=8, length=2)
    ax_cbar.set_xlabel("Mean AUCell activity z-score", fontsize=8)

    ax_spacer.axis("off")
    ax_legend.axis("off")
    legend_handles = [
        Patch(facecolor=SUBTYPE_COLORS[x], edgecolor="none", label=x) for x in SUBTYPE_ORDER
    ]
    ax_legend.legend(
        handles=legend_handles,
        title="Selected by RSS",
        loc="upper left",
        frameon=False,
        fontsize=8,
        title_fontsize=8,
    )
    fig.suptitle(
        "Subtype-associated TF regulon programs",
        x=0.42,
        y=0.99,
        fontsize=11,
        fontweight="bold",
    )
    fig.subplots_adjust(top=0.94, left=0.06, right=0.98, bottom=0.09)
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    out_format = pdf_path.suffix.lstrip(".") or "pdf"
    fig.savefig(pdf_path, format=out_format, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    outputs_dir = project_dir / "outputs"
    figures_dir = project_dir / "figures"
    main_fig_dir = figures_dir / "正文图"
    source_dir = figures_dir / "source_data" / "A_四亚型regulon热图"
    record_dir = figures_dir / "记录"
    for path in [main_fig_dir, source_dir, record_dir]:
        path.mkdir(parents=True, exist_ok=True)

    rank_tbl = pd.read_csv(outputs_dir / "scenic_rss_all_regulons_ranked_by_subtype.csv")
    selected = select_regulons(rank_tbl)
    auc = ad.read_h5ad(outputs_dir / "scenic_AUCell_regulon_activity.h5ad")
    mean_auc = mean_auc_by_subtype(auc, selected["regulon"].tolist())
    z = zscore_rows(mean_auc)

    pdf_path = main_fig_dir / "A_四亚型regulon热图.pdf"
    save_heatmap(z, selected, pdf_path)
    preview_png = os.environ.get("PANELA_PREVIEW_PNG")
    if preview_png:
        save_heatmap(z, selected, Path(preview_png))

    selected.to_csv(source_dir / "selected_regulons_by_RSS.csv", index=False)
    mean_auc.to_csv(source_dir / "mean_AUCell_by_subtype.csv")
    z.to_csv(source_dir / "mean_AUCell_zscore_by_regulon.csv")
    rank_tbl.to_csv(source_dir / "all_regulon_RSS_ranking.csv", index=False)

    manifest = pd.DataFrame(
        [
            {
                "figure": "Panel A 四亚型regulon程序热图",
                "pdf": str(pdf_path),
                "source_data_dir": str(source_dir),
                "selection_rule": f"top {TOP_N_PER_SUBTYPE} RSS-ranked regulons per subtype, duplicates assigned to subtype with highest RSS",
                "color_quantity": "mean AUCell activity by subtype, row z-score",
                "n_selected_regulons": int(selected.shape[0]),
                "subtype_order": ",".join(SUBTYPE_ORDER),
            }
        ]
    )
    manifest.to_csv(record_dir / "PanelA_四亚型regulon热图_manifest.csv", index=False)
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
