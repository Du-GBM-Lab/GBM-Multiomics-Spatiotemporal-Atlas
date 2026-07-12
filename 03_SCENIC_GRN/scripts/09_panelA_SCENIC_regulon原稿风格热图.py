#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


SUBTYPE_ORDER = ["NPC-P", "OPC-M", "MES-V", "MES-I"]
SUBTYPE_COLORS = {
    "NPC-P": "#0072B5",
    "OPC-M": "#E18727",
    "MES-V": "#20854E",
    "MES-I": "#BC3C29",
}
TOP_N_PER_SUBTYPE = 20
Z_VMIN = -1.2
Z_VMAX = 1.2
SUBTYPE_GAP_LINEWIDTH = 5.0
REGULON_GAP_LINEWIDTH = 3.0


def select_regulons(rank_tbl: pd.DataFrame) -> pd.DataFrame:
    selected = []
    rows = []
    for subtype in SUBTYPE_ORDER:
        top = rank_tbl.loc[rank_tbl["subtype"].eq(subtype)].sort_values("rank").head(
            TOP_N_PER_SUBTYPE
        )
        for _, row in top.iterrows():
            regulon = row["regulon"]
            if regulon in selected:
                continue
            selected.append(regulon)
            rows.append(
                {
                    "regulon": regulon,
                    "tf": row["tf"],
                    "assigned_subtype": subtype,
                    "assigned_subtype_rank": int(row["rank"]),
                    "assigned_subtype_rss": float(row["rss"]),
                    "selection_reason": f"top{TOP_N_PER_SUBTYPE}_RSS",
                }
            )
    out = pd.DataFrame(rows)
    out["assigned_subtype"] = pd.Categorical(
        out["assigned_subtype"], categories=SUBTYPE_ORDER, ordered=True
    )
    return out.sort_values(
        ["assigned_subtype", "assigned_subtype_rank", "regulon"]
    ).reset_index(drop=True)


def make_zscore_auc(auc: ad.AnnData, var_names: list[str]) -> ad.AnnData:
    subset = auc[:, var_names].copy()
    x = subset.X
    if not isinstance(x, np.ndarray):
        x = x.toarray()
    x = np.asarray(x, dtype=float)
    mean = x.mean(axis=0, keepdims=True)
    sd = x.std(axis=0, ddof=0, keepdims=True)
    sd[sd == 0] = 1.0
    z = (x - mean) / sd

    z_auc = ad.AnnData(X=z, obs=subset.obs.copy(), var=subset.var.copy())
    z_auc.obs_names = subset.obs_names.copy()
    z_auc.var_names = subset.var_names.copy()
    z_auc.uns["cell_type_colors"] = [SUBTYPE_COLORS[x] for x in SUBTYPE_ORDER]
    return z_auc


def main() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    outputs_dir = project_dir / "outputs"
    figures_dir = project_dir / "figures"
    fig_dir = figures_dir / "正文图"
    source_dir = figures_dir / "source_data" / "A_四亚型regulon原稿风格热图"
    record_dir = figures_dir / "记录"
    for path in [fig_dir, source_dir, record_dir]:
        path.mkdir(parents=True, exist_ok=True)

    rank_tbl = pd.read_csv(outputs_dir / "scenic_rss_all_regulons_ranked_by_subtype.csv")
    auc = ad.read_h5ad(outputs_dir / "scenic_AUCell_regulon_activity.h5ad")
    auc.obs["cell_type"] = pd.Categorical(
        auc.obs["cell_type"].astype(str), categories=SUBTYPE_ORDER, ordered=True
    )
    auc.uns["cell_type_colors"] = [SUBTYPE_COLORS[x] for x in SUBTYPE_ORDER]

    selected = select_regulons(rank_tbl)
    var_names = selected["regulon"].tolist()
    z_auc = make_zscore_auc(auc, var_names)

    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "sans-serif"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    # This mirrors the old manuscript heatmap style: cells on x-axis,
    # regulons on y-axis, no clustering, visible row labels, and stronger color contrast.
    with plt.rc_context({"font.size": 6.4}):
        ax_dict = sc.pl.heatmap(
            z_auc,
            var_names=var_names,
            groupby="cell_type",
            cmap="RdBu_r",
            swap_axes=True,
            figsize=(16, 10.2),
            show_gene_labels=True,
            dendrogram=False,
            show=False,
            vmin=Z_VMIN,
            vmax=Z_VMAX,
        )
        main_ax = ax_dict.get("heatmap_ax") or ax_dict.get("mainplot_ax")
        if main_ax is not None:
            main_ax.set_xlabel("")
            main_ax.set_ylabel("")
            main_ax.set_xticks([])
            subtype_counts = (
                z_auc.obs["cell_type"].value_counts(sort=False).reindex(SUBTYPE_ORDER).fillna(0)
            )
            for boundary in np.cumsum(subtype_counts.to_numpy())[:-1]:
                main_ax.axvline(
                    boundary - 0.5,
                    color="white",
                    linewidth=SUBTYPE_GAP_LINEWIDTH,
                    solid_capstyle="butt",
                    zorder=10,
                )
            regulon_counts = (
                selected["assigned_subtype"]
                .value_counts(sort=False)
                .reindex(SUBTYPE_ORDER)
                .fillna(0)
            )
            for boundary in np.cumsum(regulon_counts.to_numpy())[:-1]:
                main_ax.axhline(
                    boundary - 0.5,
                    color="white",
                    linewidth=REGULON_GAP_LINEWIDTH,
                    solid_capstyle="butt",
                    zorder=10,
                )
            main_ax.tick_params(axis="y", labelsize=6.4, length=0, pad=1)
            for label in main_ax.get_yticklabels():
                label.set_fontweight("bold")
                label.set_rotation(22)
                label.set_rotation_mode("anchor")
                label.set_ha("right")
                label.set_va("center")

        groupby_ax = ax_dict.get("groupby_ax")
        if groupby_ax is not None:
            groupby_ax.tick_params(axis="x", labelsize=12, length=0)
            subtype_counts = (
                z_auc.obs["cell_type"].value_counts(sort=False).reindex(SUBTYPE_ORDER).fillna(0)
            )
            for boundary in np.cumsum(subtype_counts.to_numpy())[:-1]:
                groupby_ax.axvline(
                    boundary - 0.5,
                    color="white",
                    linewidth=SUBTYPE_GAP_LINEWIDTH,
                    solid_capstyle="butt",
                    zorder=10,
                )

        colorbar_ax = ax_dict.get("heatmap_cbar_ax")
        if colorbar_ax is not None:
            colorbar_ax.tick_params(labelsize=7, length=2)

        plt.suptitle("Subtype-associated TF regulon programs", fontsize=11, fontweight="bold")
        pdf_path = fig_dir / "A_四亚型regulon原稿风格热图.pdf"
        plt.savefig(pdf_path, bbox_inches="tight", dpi=300)
        plt.close("all")

    selected = selected.copy()
    selected["row_position_1based"] = range(1, len(selected) + 1)
    selected.to_csv(source_dir / "selected_regulons_top20_by_RSS.csv", index=False)
    rank_tbl.to_csv(source_dir / "all_regulon_RSS_ranking.csv", index=False)
    pd.DataFrame(
        [
            {
                "color_quantity": "single-cell AUCell z-score per regulon",
                "z_vmin": Z_VMIN,
                "z_vmax": Z_VMAX,
                "cmap": "RdBu_r",
                "note": "z-score computed across cells for each regulon; values outside the displayed range are color-saturated",
            }
        ]
    ).to_csv(source_dir / "heatmap_color_scale.csv", index=False)

    manifest = pd.DataFrame(
        [
            {
                "figure": "Panel A 四亚型regulon原稿风格热图",
                "pdf": str(pdf_path),
                "source_data_dir": str(source_dir),
                "selection_rule": f"top {TOP_N_PER_SUBTYPE} RSS-ranked regulons per subtype, duplicates removed",
                "color_quantity": f"single-cell AUCell z-score per regulon, vmin={Z_VMIN}, vmax={Z_VMAX}",
                "n_selected_regulons": int(len(var_names)),
                "FOSL1_included": bool("FOSL1(+)" in var_names),
                "FOSL1_row_position_1based": int(var_names.index("FOSL1(+)") + 1)
                if "FOSL1(+)" in var_names
                else pd.NA,
                "subtype_order": ",".join(SUBTYPE_ORDER),
                "show_gene_labels": True,
                "dendrogram": False,
                "subtype_gap_linewidth": SUBTYPE_GAP_LINEWIDTH,
                "regulon_gap_linewidth": REGULON_GAP_LINEWIDTH,
            }
        ]
    )
    manifest.to_csv(record_dir / "PanelA_四亚型regulon原稿风格热图_manifest.csv", index=False)
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
