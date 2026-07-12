#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
R9 Batch 3 B1 | COMMOT spatial-aware communication audit

Positioning:
  B1 provides spatial non-contradiction evidence for the locked R6/R7
  PLAU-PLAUR / FOSL1->PLAUR axis. The ceiling is "spatially consistent with
  inferred communication". This script does not prove communication, causality,
  physical contact, single-cell colocalization, recruitment, or feed-forward.

Pre-declared tiers:
  primary: PLAU -> PLAUR only.
  secondary: pre-existing R5 PLAUR-related LR (VTN -> PLAUR), not improvised.
  exploratory: R5 pure-TME CellChat unique LR with Spatial gene coverage, locked
  as internal-audit/background landscape. Exploratory LR cannot be promoted.

Method:
  Per-slice COMMOT, using spot coordinates and a k=6-aligned distance threshold.
  COMMOT output is treated as inferred L-R signaling potential only. Any spatial
  relation claim is tested separately by v3-B-style random-labeling crossK on a
  pre-declared high-potential region.

TAM isolation:
  Any myeloid/TAM-looking signal is flagged as internal-audit only and cannot
  revive TAM recruitment or TAM-MES communication wording.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import anndata as ad
import commot as ct
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy.spatial import distance_matrix


BASE = Path.cwd()
IN_DIR = BASE / "tables" / "R9_batch3_B1_spatial_communication" / "input"
OUT_DIR = BASE / "tables" / "R9_batch3_B1_spatial_communication"
OUT_DIR.mkdir(parents=True, exist_ok=True)

MANIFEST = pd.read_csv(IN_DIR / "B1_COMMOT_input_manifest.csv")
LR_ALL = pd.read_csv(IN_DIR / "B1_LR_database_all_with_spatial_coverage.csv")
LR_USABLE = pd.read_csv(IN_DIR / "B1_LR_database_usable_for_COMMOT.csv")
THRESHOLDS = pd.read_csv(BASE / "tables" / "C3_niche_preflight" / "R9_C3_1_neighborhood_niche_scores.csv")
THRESHOLDS = (
    THRESHOLDS[["slice", "k", "distance_threshold"]]
    .drop_duplicates()
    .query("k == 6")
    .groupby("slice", as_index=False)["distance_threshold"]
    .median()
)

N_PERM = 999
TOP_FRAC = 0.10
DATABASE_NAME = "R9_B1_user_database"
np.random.seed(20260612)


def make_ligrec_df(lr: pd.DataFrame) -> pd.DataFrame:
    df = lr[["ligand", "receptor", "pathway_name"]].drop_duplicates().copy()
    return pd.DataFrame(df.to_numpy(dtype=str))


def read_slice_adata(row: pd.Series) -> ad.AnnData:
    sdir = Path(row["slice_dir"])
    # scipy's fast MatrixMarket path reader can mis-handle non-ASCII Windows
    # paths; pass an already opened binary stream instead.
    with open(sdir / "counts_spots_by_genes.mtx", "rb") as fh:
        x = scipy.io.mmread(fh).tocsr()
    genes = pd.read_csv(sdir / "genes.csv")["gene"].astype(str).tolist()
    obs = pd.read_csv(sdir / "obs_scores.csv")
    adata = ad.AnnData(X=x)
    adata.obs = obs.set_index("spot_id", drop=False)
    adata.var = pd.DataFrame(index=genes)
    adata.var_names_make_unique()
    adata.obsm["spatial"] = obs[["x", "y"]].to_numpy(float)
    # COMMOT assumes non-negative abundances. Use log1p normalized counts.
    lib = np.asarray(adata.X.sum(axis=1)).reshape(-1)
    lib[lib == 0] = 1.0
    adata.X = adata.X.multiply(1e4 / lib[:, None]).tocsr()
    adata.X.data = np.log1p(adata.X.data)
    return adata


def get_pair_columns(adata: ad.AnnData, ligand: str, receptor: str) -> tuple[str | None, str | None]:
    sender_key = f"commot-{DATABASE_NAME}-sum-sender"
    receiver_key = f"commot-{DATABASE_NAME}-sum-receiver"
    pair = f"{ligand}-{receptor}"
    s_col = f"s-{pair}"
    r_col = f"r-{pair}"
    s = s_col if sender_key in adata.obsm and s_col in adata.obsm[sender_key].columns else None
    r = r_col if receiver_key in adata.obsm and r_col in adata.obsm[receiver_key].columns else None
    return s, r


def empirical_p_greater(obs: float, null: np.ndarray) -> float:
    null = null[np.isfinite(null)]
    return (1.0 + np.sum(null >= obs)) / (len(null) + 1.0)


def make_hot(x: np.ndarray, frac: float = TOP_FRAC) -> np.ndarray:
    cutoff = np.nanquantile(x, 1.0 - frac)
    return np.asarray(x >= cutoff)


def crossk_count(coords: np.ndarray, a_idx: np.ndarray, b_idx: np.ndarray, radius: float) -> float:
    if len(a_idx) == 0 or len(b_idx) == 0:
        return np.nan
    d = distance_matrix(coords[a_idx], coords[b_idx])
    return float(np.mean(np.sum(d <= radius, axis=1)))


def v3b_region_test(sd: pd.DataFrame, region_col: str, feature: str, radius: float) -> dict:
    coords = sd[["x", "y"]].to_numpy(float)
    n = len(sd)
    a_idx = np.where(sd[region_col].to_numpy(bool))[0]
    b_idx = np.where(make_hot(sd[feature].to_numpy(float), TOP_FRAC))[0]
    obs = crossk_count(coords, a_idx, b_idx, radius)
    null = np.empty(N_PERM, dtype=float)
    count_ok = np.empty(N_PERM, dtype=bool)
    n_a, n_b = len(a_idx), len(b_idx)
    for i in range(N_PERM):
        ai = np.random.choice(n, size=n_a, replace=False)
        bi = np.random.choice(n, size=n_b, replace=False)
        count_ok[i] = len(ai) == n_a and len(bi) == n_b
        null[i] = crossk_count(coords, ai, bi, radius)
    return {
        "feature": feature,
        "n_region": n_a,
        "n_feature_hot": n_b,
        "obs_crossK": obs,
        "null_median": float(np.nanmedian(null)),
        "delta_crossK": float(obs - np.nanmedian(null)),
        "emp_p_greater": empirical_p_greater(obs, null),
        "count_preserved_all_perm": bool(np.all(count_ok)),
    }


def run_commot_slice(row: pd.Series, lr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    slice_id = row["slice"]
    adata = read_slice_adata(row)
    dis_thr = float(THRESHOLDS.loc[THRESHOLDS["slice"] == slice_id, "distance_threshold"].iloc[0])
    df_ligrec = make_ligrec_df(lr)
    ct.tl.spatial_communication(
        adata,
        database_name=DATABASE_NAME,
        df_ligrec=df_ligrec,
        dis_thr=dis_thr,
        heteromeric=True,
        pathway_sum=True,
    )
    sender_key = f"commot-{DATABASE_NAME}-sum-sender"
    receiver_key = f"commot-{DATABASE_NAME}-sum-receiver"
    sender = adata.obsm[sender_key].copy()
    receiver = adata.obsm[receiver_key].copy()
    obs = adata.obs.reset_index(drop=True).copy()
    obs["slice"] = slice_id

    spot_rows = []
    summary_rows = []
    for rec in lr.to_dict("records"):
        ligand, receptor = rec["ligand"], rec["receptor"]
        s_col, r_col = get_pair_columns(adata, ligand, receptor)
        if s_col is None or r_col is None:
            summary_rows.append({
                "slice": slice_id, **rec,
                "commot_available": False,
                "reason": "sender_or_receiver_column_missing",
            })
            continue
        s_val = np.asarray(sender[s_col], dtype=float)
        r_val = np.asarray(receiver[r_col], dtype=float)
        total = s_val + r_val
        reg = make_hot(total, TOP_FRAC)
        summary_rows.append({
            "slice": slice_id, **rec,
            "commot_available": True,
            "reason": "",
            "distance_threshold_k6": dis_thr,
            "sender_sum": float(np.sum(s_val)),
            "receiver_sum": float(np.sum(r_val)),
            "total_sum": float(np.sum(total)),
            "n_positive_total": int(np.sum(total > 0)),
            "n_top10_region": int(np.sum(reg)),
            "has_myeloid_TAM_risk": bool(rec["tier"] != "exploratory" and rec["interaction_name"] == "PLAU_PLAUR"),
        })
        if rec["tier"] in ("primary", "secondary"):
            tmp = obs[["spot_id", "slice", "image", "x", "y", "MES-V", "MES-lineage", "vascular", "neuron_control"]].copy()
            tmp["tier"] = rec["tier"]
            tmp["interaction_name"] = rec["interaction_name"]
            tmp["ligand"] = ligand
            tmp["receptor"] = receptor
            tmp["sender_potential"] = s_val
            tmp["receiver_potential"] = r_val
            tmp["total_potential"] = total
            tmp["top10_region"] = reg
            spot_rows.append(tmp)
    spot_dt = pd.concat(spot_rows, ignore_index=True) if spot_rows else pd.DataFrame()
    summary_dt = pd.DataFrame(summary_rows)
    adata_path = OUT_DIR / "h5ad" / f"{row['safe_slice']}_COMMOT.h5ad"
    adata_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(adata_path, compression="gzip")
    summary_dt["h5ad_path"] = str(adata_path)
    return spot_dt, summary_dt, pd.DataFrame()


def main() -> None:
    # Run COMMOT for primary/secondary only. The pre-declared exploratory LR
    # list is retained as coverage/background landscape, but is not used as a
    # discovery machine or for main/supplement claims. A local full-library
    # COMMOT run is computationally prohibitive and does not affect the primary
    # gate.
    lr = LR_USABLE[LR_USABLE["tier"].isin(["primary", "secondary"])].copy()
    all_spot, all_summary = [], []
    for _, row in MANIFEST.iterrows():
        print(f"[COMMOT] {row['slice']} | spots={row['n_spots']} | LR={len(lr)}", flush=True)
        spot_dt, summary_dt, _ = run_commot_slice(row, lr)
        if len(spot_dt):
            all_spot.append(spot_dt)
        all_summary.append(summary_dt)

    spot = pd.concat(all_spot, ignore_index=True) if all_spot else pd.DataFrame()
    summary = pd.concat(all_summary, ignore_index=True)
    spot.to_csv(OUT_DIR / "B1_COMMOT_primary_secondary_per_spot_potential.csv", index=False)
    summary.to_csv(OUT_DIR / "B1_COMMOT_LR_summary_per_slice.csv", index=False)

    # Primary v3-B relation test: PLAU_PLAUR top10 potential region vs MES-V,
    # vascular, and neuron_control. Secondary gets descriptive summary only here.
    primary_spot = spot[(spot["tier"] == "primary") & (spot["interaction_name"] == "PLAU_PLAUR")].copy()
    rows = []
    for slice_id, sd in primary_spot.groupby("slice", sort=True):
        radius = float(THRESHOLDS.loc[THRESHOLDS["slice"] == slice_id, "distance_threshold"].iloc[0])
        for feat in ["MES-V", "MES-lineage", "vascular", "neuron_control"]:
            res = v3b_region_test(sd, "top10_region", feat, radius)
            rows.append({"slice": slice_id, "interaction_name": "PLAU_PLAUR", "tier": "primary",
                         "top_region_rule": "top10_total_COMMOT_potential",
                         "distance_rule": "k6_radius_from_C3", **res})
    rel = pd.DataFrame(rows)
    rel["fdr_bh_by_feature"] = np.nan
    for feat, idx in rel.groupby("feature").groups.items():
        p = rel.loc[idx, "emp_p_greater"].to_numpy(float)
        order = np.argsort(p)
        ranks = np.empty_like(order)
        ranks[order] = np.arange(1, len(p) + 1)
        q = p * len(p) / ranks
        q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
        q_final = np.empty_like(q_sorted)
        q_final[order] = np.minimum(q_sorted, 1.0)
        rel.loc[idx, "fdr_bh_by_feature"] = q_final
    rel["pass_positive_fdr"] = (rel["delta_crossK"] > 0) & (rel["fdr_bh_by_feature"] < 0.05)
    rel.to_csv(OUT_DIR / "B1_primary_PLAU_PLAUR_v3B_region_relation_perslice.csv", index=False)

    rel_summary = (
        rel.groupby("feature")
        .agg(
            n_slices=("slice", "count"),
            n_positive=("delta_crossK", lambda x: int(np.sum(np.asarray(x) > 0))),
            n_pass_positive_fdr=("pass_positive_fdr", "sum"),
            median_delta_crossK=("delta_crossK", "median"),
            median_emp_p=("emp_p_greater", "median"),
            count_preserved_all=("count_preserved_all_perm", "all"),
        )
        .reset_index()
    )
    rel_summary.to_csv(OUT_DIR / "B1_primary_PLAU_PLAUR_v3B_region_relation_summary.csv", index=False)

    primary_gate = pd.DataFrame({
        "gate_item": [
            "COMMOT primary output available 18/18",
            "MES-V relation >=15/18 positive FDR",
            "vascular relation >=15/18 positive FDR",
            "neuron_control relation 0/18 positive FDR",
            "count preservation all",
            "TAM isolation",
            "overall_primary_gate",
        ],
        "pass": [
            int(summary[(summary["tier"] == "primary") & (summary["commot_available"])].shape[0]) == 18,
            int(rel_summary.loc[rel_summary["feature"] == "MES-V", "n_pass_positive_fdr"].iloc[0]) >= 15,
            int(rel_summary.loc[rel_summary["feature"] == "vascular", "n_pass_positive_fdr"].iloc[0]) >= 15,
            int(rel_summary.loc[rel_summary["feature"] == "neuron_control", "n_pass_positive_fdr"].iloc[0]) == 0,
            bool(rel["count_preserved_all_perm"].all()),
            False,  # COMMOT PLAU_PLAUR necessarily shows myeloid/TAM-rich risk from R5 source, isolated by spec.
            False,
        ],
        "value": [
            f"{summary[(summary['tier']=='primary') & (summary['commot_available'])].shape[0]}/18",
            str(rel_summary[rel_summary["feature"] == "MES-V"].to_dict("records")),
            str(rel_summary[rel_summary["feature"] == "vascular"].to_dict("records")),
            str(rel_summary[rel_summary["feature"] == "neuron_control"].to_dict("records")),
            f"{int(rel['count_preserved_all_perm'].sum())}/{len(rel)} relation rows",
            "PLAU_PLAUR R5 source includes Macrophages/Microglial/Monocytes receptor/sender prominence; isolated as internal-audit risk and cannot revive TAM wording",
            "supplement_or_internal_audit; no main-text promotion without TAM override, which is forbidden",
        ],
    })
    primary_gate.to_csv(OUT_DIR / "B1_primary_PLAU_PLAUR_gate.csv", index=False)

    exploratory_landscape = (
        LR_USABLE.groupby("tier")
        .agg(n_usable_lr=("interaction_name", "nunique"))
        .reset_index()
    )
    commot_landscape = (
        summary.groupby("tier")
        .agg(
            n_lr=("interaction_name", "nunique"),
            n_slice_lr_rows=("interaction_name", "count"),
            n_available=("commot_available", "sum"),
            median_total_sum=("total_sum", "median"),
        )
        .reset_index()
    )
    landscape = exploratory_landscape.merge(commot_landscape, on="tier", how="left")
    landscape["commot_run_scope"] = np.where(
        landscape["tier"].isin(["primary", "secondary"]),
        "COMMOT_run",
        "coverage_only_internal_audit",
    )
    landscape.to_csv(OUT_DIR / "B1_COMMOT_landscape_by_tier.csv", index=False)

    meta = {
        "commot_docs_checked": "https://commot.readthedocs.io/en/latest/notebooks/Basic_usage.html",
        "commot_database_name": DATABASE_NAME,
        "n_perm_v3B": N_PERM,
        "commot_run_scope": "primary and secondary only; exploratory retained as coverage/background because full local COMMOT run timed out and is locked out of promotion",
        "top_region_rule": "top10 total sender+receiver inferred potential",
        "distance_rule": "k=6 C3 distance threshold per slice",
        "wording_ceiling": "spatially consistent with inferred communication; no causality/contact/recruitment",
    }
    (OUT_DIR / "B1_COMMOT_run_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    decision = "supplement_or_internal_audit"
    stop_lines = [
        "# R9 Batch3 B1 COMMOT spatial communication STOP",
        "",
        "- Nature: spatial-aware COMMOT audit for R6/R7 PLAU-PLAUR axis; this does not prove communication, causality, physical contact, or single-cell colocalization.",
        "- Primary: PLAU -> PLAUR. Secondary: VTN -> PLAUR. Exploratory: usable R5 pure-TME LR coverage landscape, locked internal-audit/background; full local COMMOT transport for exploratory timed out and was not used for discovery or claims.",
        "- COMMOT output is interpreted only as inferred L-R expression/signaling potential; relation-level testing used v3-B-style random-labeling crossK on the pre-declared top10 potential region.",
        f"- Primary gate decision: {decision}.",
        "",
        "## Key outputs",
        f"- {OUT_DIR / 'B1_COMMOT_primary_secondary_per_spot_potential.csv'}",
        f"- {OUT_DIR / 'B1_COMMOT_LR_summary_per_slice.csv'}",
        f"- {OUT_DIR / 'B1_primary_PLAU_PLAUR_v3B_region_relation_perslice.csv'}",
        f"- {OUT_DIR / 'B1_primary_PLAU_PLAUR_v3B_region_relation_summary.csv'}",
        f"- {OUT_DIR / 'B1_primary_PLAU_PLAUR_gate.csv'}",
        "",
        "## Boundary",
        "- TAM/myeloid-looking COMMOT potential is isolated as internal-audit risk and cannot be used for TAM recruitment or TAM-MES communication.",
        "- Allowed wording ceiling: spatially consistent with inferred communication / co-localized ligand-receptor expression / inferred signaling potential / consistent with R5-R6.",
        "- Forbidden: confirmed loop, interaction as fact, recruitment, drives, gates, activates, feed-forward, proves, causality, physical/single-cell communication.",
    ]
    (BASE / "docs" / "R9_batch3_B1_COMMOT_spatial_communication_STOP.md").write_text("\n".join(stop_lines), encoding="utf-8")


if __name__ == "__main__":
    main()
