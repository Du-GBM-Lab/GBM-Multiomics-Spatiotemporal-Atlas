#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


SUBTYPE_ORDER = ["NPC-P", "OPC-M", "MES-V", "MES-I"]
TOP_N_PER_SUBTYPE = 20
FORCE_INCLUDE = ["FOSL1(+)"]


def select_regulons(rank_tbl: pd.DataFrame, auc: ad.AnnData) -> pd.DataFrame:
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

    for regulon in FORCE_INCLUDE:
        if regulon in selected or regulon not in auc.var_names:
            continue
        hit = rank_tbl.loc[rank_tbl["regulon"].eq(regulon)].copy()
        if hit.empty:
            continue
        row = hit.loc[hit["subtype"].eq("MES-V")].iloc[0]
        selected.append(regulon)
        rows.append(
            {
                "regulon": regulon,
                "tf": row["tf"],
                "assigned_subtype": row["subtype"],
                "assigned_subtype_rank": int(row["rank"]),
                "assigned_subtype_rss": float(row["rss"]),
                "selection_reason": "force_include_for_MESV_context",
            }
        )

    out = pd.DataFrame(rows)
    out["assigned_subtype"] = pd.Categorical(
        out["assigned_subtype"], categories=SUBTYPE_ORDER, ordered=True
    )
    return out.sort_values(
        ["assigned_subtype", "assigned_subtype_rank", "regulon"]
    ).reset_index(drop=True)


def main() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    outputs_dir = project_dir / "outputs"
    figures_dir = project_dir / "figures"
    fig_dir = figures_dir / "正文图"
    source_dir = figures_dir / "source_data" / "A_四亚型regulon扁大热图"
    record_dir = figures_dir / "记录"
    for path in [fig_dir, source_dir, record_dir]:
        path.mkdir(parents=True, exist_ok=True)

    rank_tbl = pd.read_csv(outputs_dir / "scenic_rss_all_regulons_ranked_by_subtype.csv")
    auc = ad.read_h5ad(outputs_dir / "scenic_AUCell_regulon_activity.h5ad")
    auc.obs["cell_type"] = pd.Categorical(
        auc.obs["cell_type"].astype(str), categories=SUBTYPE_ORDER, ordered=True
    )

    selected = select_regulons(rank_tbl, auc)
    var_names = selected["regulon"].tolist()

    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "sans-serif"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    with plt.rc_context({"font.size": 7}):
        ax_dict = sc.pl.heatmap(
            auc,
            var_names=var_names,
            groupby="cell_type",
            cmap="RdYlBu_r",
            standard_scale="var",
            swap_axes=True,
            figsize=(18, 6),
            show_gene_labels=False,
            dendrogram=False,
            show=False,
        )
        main_ax = ax_dict.get("heatmap_ax") or ax_dict.get("mainplot_ax")
        if main_ax is not None:
            main_ax.set_xlabel("Malignant cells grouped by subtype", fontsize=9)
            main_ax.set_ylabel(f"{len(var_names)} RSS-selected regulons", fontsize=9)
            main_ax.set_xticks([])
        plt.suptitle("Subtype-associated TF regulon programs", fontsize=11, fontweight="bold")
        pdf_path = fig_dir / "A_四亚型regulon扁大热图.pdf"
        plt.savefig(pdf_path, bbox_inches="tight", dpi=300)
        plt.close("all")

    selected = selected.copy()
    selected["row_position_1based"] = range(1, len(selected) + 1)
    selected.to_csv(source_dir / "selected_regulons_top20_by_RSS.csv", index=False)
    rank_tbl.to_csv(source_dir / "all_regulon_RSS_ranking.csv", index=False)

    manifest = pd.DataFrame(
        [
            {
                "figure": "Panel A 四亚型regulon扁大热图",
                "pdf": str(pdf_path),
                "source_data_dir": str(source_dir),
                "selection_rule": f"top {TOP_N_PER_SUBTYPE} RSS-ranked regulons per subtype, duplicates removed; FOSL1 forced only if absent",
                "color_quantity": "single-cell AUCell activity, scanpy standard_scale='var'",
                "n_selected_regulons": int(len(var_names)),
                "FOSL1_included": bool("FOSL1(+)" in var_names),
                "FOSL1_row_position_1based": int(var_names.index("FOSL1(+)") + 1)
                if "FOSL1(+)" in var_names
                else pd.NA,
                "subtype_order": ",".join(SUBTYPE_ORDER),
                "show_gene_labels": False,
            }
        ]
    )
    manifest.to_csv(record_dir / "PanelA_四亚型regulon扁大热图_manifest.csv", index=False)
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
