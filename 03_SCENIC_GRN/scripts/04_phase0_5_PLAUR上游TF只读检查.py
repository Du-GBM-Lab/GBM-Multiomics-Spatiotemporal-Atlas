#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path

import pandas as pd


def main() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    outputs_dir = project_dir / "outputs"
    source_data_dir = project_dir / "figures" / "source_data" / "Phase0_5_PLAUR上游TF"
    source_data_dir.mkdir(parents=True, exist_ok=True)

    edgelist = pd.read_csv(outputs_dir / "scenic_grn_edgelist.csv")
    required = {"TF", "target", "importance"}
    missing = required.difference(edgelist.columns)
    if missing:
        raise ValueError(f"Missing columns in GRN edgelist: {sorted(missing)}")

    plaur_edges = edgelist.loc[edgelist["target"].astype(str).eq("PLAUR")].copy()
    plaur_edges["importance"] = pd.to_numeric(plaur_edges["importance"], errors="coerce")
    plaur_by_tf = (
        plaur_edges.dropna(subset=["importance"])
        .groupby("TF", as_index=False)["importance"]
        .max()
        .sort_values(["importance", "TF"], ascending=[False, True])
        .reset_index(drop=True)
    )
    plaur_by_tf["rank"] = range(1, len(plaur_by_tf) + 1)
    plaur_by_tf = plaur_by_tf.loc[:, ["rank", "TF", "importance"]]

    fosl1 = plaur_by_tf.loc[plaur_by_tf["TF"].eq("FOSL1")].copy()
    summary = pd.DataFrame(
        [
            {
                "target": "PLAUR",
                "n_upstream_tf": int(plaur_by_tf["TF"].nunique()),
                "FOSL1_found": bool(len(fosl1)),
                "FOSL1_importance": float(fosl1["importance"].iloc[0]) if len(fosl1) else pd.NA,
                "FOSL1_rank_among_PLAUR_upstream_TF": int(fosl1["rank"].iloc[0])
                if len(fosl1)
                else pd.NA,
                "ranking_rule": "unique TF ranked by max GRN importance for target PLAUR",
            }
        ]
    )

    ranking_path = outputs_dir / "scenic_phase0_5_PLAUR_upstream_tf_ranking.csv"
    summary_path = outputs_dir / "scenic_phase0_5_FOSL1_PLAUR_rank_summary.csv"
    plaur_by_tf.to_csv(ranking_path, index=False)
    summary.to_csv(summary_path, index=False)
    plaur_by_tf.to_csv(source_data_dir / ranking_path.name, index=False)
    summary.to_csv(source_data_dir / summary_path.name, index=False)

    print(summary.to_string(index=False))
    print(f"ranking={ranking_path}")
    print(f"summary={summary_path}")


if __name__ == "__main__":
    main()
