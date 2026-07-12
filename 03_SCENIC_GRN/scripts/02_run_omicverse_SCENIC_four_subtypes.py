#!/usr/bin/env python

from __future__ import annotations

import os
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import io

if not hasattr(np, "float"):
    np.float = float


SUBTYPE_ORDER = ["NPC-P", "OPC-M", "MES-V", "MES-I"]


def prepare_tf_names(motif_path: Path, adata: ad.AnnData, outputs_dir: Path) -> list[str]:
    motif_tbl = pd.read_csv(motif_path, sep="\t", dtype=str, low_memory=False)
    if "gene_name" not in motif_tbl.columns:
        raise ValueError(f"Motif annotation missing gene_name column: {motif_path}")

    motif_tfs = sorted(set(motif_tbl["gene_name"].dropna().astype(str)))
    genes = set(adata.var_names.astype(str))
    tf_names = [tf for tf in motif_tfs if tf in genes]
    if len(tf_names) < 100:
        raise ValueError(
            f"Too few TFs after intersecting motif annotation with HVG genes: {len(tf_names)}"
        )

    pd.Series(tf_names, name="tf").to_csv(
        outputs_dir / "scenic_tf_names_from_motif_annotation.csv", index=False
    )
    return tf_names


def regulon_name_to_tf(name: str) -> str:
    return str(name).replace("(+)", "").replace("(-)", "")


def load_inputs(data_dir: Path) -> ad.AnnData:
    with open(data_dir / "gbm_counts.mtx", "rb") as fh:
        counts = io.mmread(fh).T.tocsr()
    genes = pd.read_csv(data_dir / "gbm_genes.csv")["gene"].astype(str).tolist()
    barcodes = pd.read_csv(data_dir / "gbm_barcodes.csv")["cell"].astype(str).tolist()
    meta = pd.read_csv(data_dir / "gbm_metadata.csv", index_col=0)
    umap = pd.read_csv(data_dir / "gbm_umap.csv", index_col=0)

    adata = ad.AnnData(X=counts)
    adata.obs_names = barcodes
    adata.var_names = genes
    adata.var_names_make_unique()
    adata.obs = meta.loc[adata.obs_names].copy()
    adata.obsm["X_umap"] = umap.loc[adata.obs_names].values
    adata.layers["counts"] = adata.X.copy()

    observed = set(adata.obs["cell_type"].dropna().unique())
    expected = set(SUBTYPE_ORDER)
    if observed != expected:
        raise ValueError(f"cell_type is not current four-subtype labels: {observed}")
    adata.obs["cell_type"] = pd.Categorical(adata.obs["cell_type"], categories=SUBTYPE_ORDER)
    return adata


def compute_rss(auc_adata: ad.AnnData, outputs_dir: Path) -> pd.DataFrame:
    try:
        from pyscenic.rss import regulon_specificity_scores
    except ImportError as exc:
        raise ImportError("pyscenic.rss.regulon_specificity_scores is required for explicit RSS.") from exc

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
    module_dir = project_dir.parent
    data_dir = project_dir / "data"
    outputs_dir = project_dir / "outputs"
    figures_dir = project_dir / "figures"
    source_data_dir = figures_dir / "source_data"
    logs_dir = project_dir / "logs"
    for path in [outputs_dir, figures_dir, source_data_dir, logs_dir]:
        path.mkdir(parents=True, exist_ok=True)

    db_dir = data_dir / "scenic_db"
    db_glob = str(db_dir / "*feather")
    motif_path = db_dir / "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    ranking_files = glob.glob(db_glob)
    if len(ranking_files) < 2:
        raise FileNotFoundError(f"Expected at least two ranking feather files under {db_dir}")
    if not motif_path.exists():
        raise FileNotFoundError(f"Missing motif annotation table: {motif_path}")

    import omicverse as ov

    adata = load_inputs(data_dir)
    adata.write_h5ad(outputs_dir / "gbm_assembled四亚型.h5ad")

    sc.pp.filter_genes(adata, min_cells=50)
    non_mt = [name for name in adata.var_names if not name.upper().startswith("MT-")]
    adata = adata[:, non_mt].copy()
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=10000,
        flavor="seurat_v3",
        layer="counts",
        subset=True,
    )
    adata.write_h5ad(outputs_dir / "gbm_SCENIC输入_HVG10000.h5ad")
    tf_names = prepare_tf_names(motif_path, adata, outputs_dir)

    scenic_obj = ov.single.SCENIC(
        adata=adata,
        db_glob=db_glob,
        motif_path=str(motif_path),
        n_jobs=16,
    )
    edgelist = scenic_obj.cal_grn(method="grnboost2", layer="counts", tf_names=tf_names)
    edgelist.to_csv(outputs_dir / "scenic_grn_edgelist.csv", index=False)

    scenic_obj.cal_regulons()
    scenic_obj.ad_auc_mtx.obs["cell_type"] = adata.obs.loc[
        scenic_obj.ad_auc_mtx.obs_names, "cell_type"
    ].astype(str)
    scenic_obj.ad_auc_mtx.obs["cell_type"] = pd.Categorical(
        scenic_obj.ad_auc_mtx.obs["cell_type"], categories=SUBTYPE_ORDER
    )
    scenic_obj.ad_auc_mtx.write_h5ad(outputs_dir / "scenic_AUCell_regulon_activity.h5ad")

    rss = compute_rss(scenic_obj.ad_auc_mtx, outputs_dir)
    save_top_regulons(rss, outputs_dir, source_data_dir)
    edgelist.to_csv(source_data_dir / "scenic_edgelist_final.csv", index=False)

    with open(logs_dir / "02_run_SCENIC_summary.txt", "w", encoding="utf-8") as fh:
        fh.write(f"n_cells={adata.n_obs}\n")
        fh.write(f"n_genes_after_hvg={adata.n_vars}\n")
        fh.write(f"n_tf_candidates={len(tf_names)}\n")
        fh.write("grn_method=grnboost2\n")
        fh.write(f"n_regulons={scenic_obj.ad_auc_mtx.n_vars}\n")
        fh.write(f"n_edges={len(edgelist)}\n")
        fh.write("subtype_order=" + ",".join(SUBTYPE_ORDER) + "\n")


if __name__ == "__main__":
    main()
