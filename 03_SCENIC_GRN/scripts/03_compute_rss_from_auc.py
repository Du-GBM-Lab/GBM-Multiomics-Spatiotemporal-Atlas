#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

if not hasattr(np, "float"):
    np.float = float

SUBTYPE_ORDER = ["NPC-P", "OPC-M", "MES-V", "MES-I"]


def regulon_name_to_tf(name: str) -> str:
    return str(name).replace("(+)", "").replace("(-)", "")


def compute_rss(auc_adata: ad.AnnData, outputs_dir: Path) -> pd.DataFrame:
    from pyscenic.rss import regulon_specificity_scores

    auc_x = auc_adata.X
    if not isinstance(auc_x, np.ndarray):
        auc_x = auc_x.toarray()
    auc_df = pd.DataFrame(auc_x, index=auc_adata.obs_names, columns=auc_adata.var_names)
    cell_type = auc_adata.obs["cell_type"].astype(str)
    rss = regulon_specificity_scores(auc_df, cell_type)

    if set(SUBTYPE_ORDER).issubset(rss.columns):
        rss_corrected = rss.loc[:, SUBTYPE_ORDER]
    elif set(SUBTYPE_ORDER).issubset(rss.index):
        rss_corrected = rss.T.loc[:, SUBTYPE_ORDER]
    else:
        raise ValueError(f"RSS output shape not recognized: {rss.shape}")

    rss.to_csv(outputs_dir / "scenic_rss_scores_raw.csv")
    rss_corrected.to_csv(outputs_dir / "scenic_rss_scores_regulon_by_subtype.csv")
    return rss_corrected


def save_top_regulons(rss: pd.DataFrame, outputs_dir: Path, source_data_dir: Path) -> pd.DataFrame:
    rows = []
    for subtype in SUBTYPE_ORDER:
        top = rss[subtype].sort_values(ascending=False)
        for rank, (regulon, score) in enumerate(top.items(), start=1):
            rows.append(
                {
                    "subtype": subtype,
                    "rank": rank,
                    "regulon": regulon,
                    "tf": regulon_name_to_tf(regulon),
                    "rss": score,
                }
            )
    top_tbl = pd.DataFrame(rows)
    top_tbl.to_csv(outputs_dir / "scenic_rss_all_regulons_ranked_by_subtype.csv", index=False)
    top_tbl.groupby("subtype").head(30).to_csv(
        source_data_dir / "scenic_top30_regulons_by_subtype.csv", index=False
    )
    return top_tbl


def main() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    outputs_dir = project_dir / "outputs"
    source_data_dir = project_dir / "figures" / "source_data"
    logs_dir = project_dir / "logs"
    source_data_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    auc = ad.read_h5ad(outputs_dir / "scenic_AUCell_regulon_activity.h5ad")
    rss = compute_rss(auc, outputs_dir)
    top = save_top_regulons(rss, outputs_dir, source_data_dir)

    with open(logs_dir / "03_compute_RSS_from_AUCell_summary.txt", "w", encoding="utf-8") as fh:
        fh.write(f"n_cells={auc.n_obs}\n")
        fh.write(f"n_regulons={auc.n_vars}\n")
        fh.write(f"n_ranked_rows={len(top)}\n")
        fh.write("subtype_order=" + ",".join(SUBTYPE_ORDER) + "\n")


if __name__ == "__main__":
    main()
