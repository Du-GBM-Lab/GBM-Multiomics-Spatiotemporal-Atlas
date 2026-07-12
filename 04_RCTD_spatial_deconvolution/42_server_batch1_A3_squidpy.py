#!/usr/bin/env python3
"""
R9 batch 1 | A3 squidpy neighborhood enrichment / co-occurrence on server.

Boundaries:
- This is a landscape/sensitivity analysis, not the locked v3-B null.
- Visium spots are not single cells.
- Any result following neuron_control is downgraded to internal-audit.
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd

import anndata as ad
import squidpy as sq

BASE = Path("/home/data/t010639/projects/GBM_R9_spatial_RCTD")
score_path = BASE / "tables" / "R9_stage3_signature_scores" / "R9_stage3_IVY_hypoxia_scores_per_spot.csv"
out_dir = BASE / "outputs" / "R9_batch1_unbiased_landscape" / "A3_squidpy"
tab_dir = BASE / "tables" / "R9_batch1_unbiased_landscape" / "A3_squidpy"
out_dir.mkdir(parents=True, exist_ok=True)
tab_dir.mkdir(parents=True, exist_ok=True)

top_frac = 0.10
n_perms = 999
n_neighs = 6

dt = pd.read_csv(score_path)
required = [
    "spot_id", "slice", "image", "x", "y",
    "MES-lineage", "MES-V", "MES-I", "vascular", "myeloid", "neuron_control",
]
missing = [c for c in required if c not in dt.columns]
if missing:
    raise SystemExit(f"Missing required columns: {missing}")

def top_label_by_slice(frame: pd.DataFrame, columns, prefix: str) -> pd.Series:
    out = pd.Series("other", index=frame.index, dtype="object")
    for sl, idx in frame.groupby("slice").groups.items():
        sub = frame.loc[idx, columns]
        q = sub.quantile(1 - top_frac)
        hit = sub.ge(q, axis=1)
        max_col = sub.idxmax(axis=1)
        any_hit = hit.any(axis=1)
        out.loc[idx[any_hit]] = prefix + "_" + max_col.loc[idx[any_hit]].astype(str)
    return out

dt["malignant_high_class"] = top_label_by_slice(dt, ["MES-lineage", "MES-V", "MES-I"], "mal")
dt["context_high_class"] = top_label_by_slice(dt, ["vascular", "myeloid", "neuron_control"], "ctx")
dt["combined_class"] = np.where(
    dt["malignant_high_class"] != "other",
    dt["malignant_high_class"],
    dt["context_high_class"],
)

obs = dt[["spot_id", "slice", "image", "combined_class", "malignant_high_class", "context_high_class"]].copy()
obs.index = obs["spot_id"].astype(str)
for col in ["slice", "combined_class", "malignant_high_class", "context_high_class"]:
    obs[col] = obs[col].astype("category")
X = np.zeros((obs.shape[0], 1), dtype=np.float32)
adata = ad.AnnData(X=X, obs=obs)
adata.obsm["spatial"] = dt[["x", "y"]].to_numpy(dtype=float)

sq.gr.spatial_neighbors(
    adata,
    spatial_key="spatial",
    library_key="slice",
    coord_type="generic",
    n_neighs=n_neighs,
)
sq.gr.nhood_enrichment(
    adata,
    cluster_key="combined_class",
    n_perms=n_perms,
    seed=1,
)
sq.gr.co_occurrence(
    adata,
    cluster_key="combined_class",
)

cats = list(adata.obs["combined_class"].cat.categories)
z = adata.uns["combined_class_nhood_enrichment"]["zscore"]
z_df = pd.DataFrame(z, index=cats, columns=cats)
z_long = z_df.reset_index(names="source_class").melt(
    id_vars="source_class", var_name="target_class", value_name="nhood_enrichment_z"
)
z_long.to_csv(tab_dir / "A3_nhood_enrichment_zscore_long.csv", index=False)

occ = adata.uns["combined_class_co_occurrence"]
with open(tab_dir / "A3_co_occurrence_uns_keys.json", "w", encoding="utf-8") as fh:
    json.dump({k: str(type(v)) for k, v in occ.items()}, fh, indent=2)

occ_arr = occ["occ"]
interval = occ["interval"]
occ_rows = []
for source_idx, source_class in enumerate(cats):
    for target_idx, target_class in enumerate(cats):
        for distance_idx, value in enumerate(occ_arr[source_idx, target_idx, :]):
            occ_rows.append({
                "source_class": source_class,
                "target_class": target_class,
                "distance_index": distance_idx,
                "distance": float(interval[distance_idx]) if distance_idx < len(interval) else np.nan,
                "co_occurrence": float(value),
            })
pd.DataFrame(occ_rows).to_csv(tab_dir / "A3_co_occurrence_long.csv", index=False)

adata.write_h5ad(out_dir / "A3_squidpy_input_and_results.h5ad")

params = pd.DataFrame({
    "parameter": ["input", "top_frac", "n_neighs", "n_perms", "cluster_key", "boundary"],
    "value": [
        str(score_path), top_frac, n_neighs, n_perms, "combined_class",
        "squidpy permutation landscape only; neuron_control required for interpretation",
    ],
})
params.to_csv(tab_dir / "A3_squidpy_parameters.csv", index=False)

print("[STOP A3 squidpy] Neighborhood enrichment and co-occurrence source data written. Review neuron controls before any claim.")
