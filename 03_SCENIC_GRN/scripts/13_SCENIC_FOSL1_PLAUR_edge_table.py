#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    project_dir = Path(__file__).resolve().parents[1]
    out_dir = project_dir / "outputs"
    fig_dir = project_dir / "figures" / "补充图"
    source_dir = project_dir / "figures" / "source_data" / "FOSL1_PLAUR_SCENIC_edge"
    record_dir = project_dir / "figures" / "记录"
    for path in [fig_dir, source_dir, record_dir]:
        path.mkdir(parents=True, exist_ok=True)

    summary = pd.read_csv(out_dir / "scenic_phase0_5_FOSL1_PLAUR_rank_summary.csv")
    edge = summary.iloc[0]

    table = pd.DataFrame(
        [
            {
                "Inferred edge": "FOSL1 -> PLAUR",
                "SCENIC status": "found",
                "importance": f"{float(edge['FOSL1_importance']):.3f}",
                "rank among PLAUR upstream TFs": int(edge["FOSL1_rank_among_PLAUR_upstream_TF"]),
                "interpretation": "weak nomination",
            },
            {
                "Inferred edge": "FOSL1 -> PLAU",
                "SCENIC status": "not recovered",
                "importance": "NA",
                "rank among PLAUR upstream TFs": "NA",
                "interpretation": "not used as target claim",
            },
        ]
    )

    table.to_csv(source_dir / "FOSL1_PLAUR_SCENIC_edge_table.csv", index=False)
    summary.to_csv(source_dir / "scenic_phase0_5_FOSL1_PLAUR_rank_summary.csv", index=False)
    pd.read_csv(out_dir / "scenic_phase0_5_PLAUR_upstream_tf_ranking.csv").to_csv(
        source_dir / "scenic_phase0_5_PLAUR_upstream_tf_ranking.csv", index=False
    )

    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "sans-serif"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    fig, ax = plt.subplots(figsize=(7.6, 1.8))
    ax.axis("off")
    ax.set_title("SCENIC nominated a weak FOSL1-PLAUR edge", fontsize=10.5, fontweight="bold", pad=6)

    display = table.rename(
        columns={
            "Inferred edge": "edge",
            "SCENIC status": "status",
            "rank among PLAUR upstream TFs": "PLAUR TF rank",
        }
    )
    cols = ["edge", "status", "importance", "PLAUR TF rank", "interpretation"]
    tbl = ax.table(
        cellText=display[cols].values,
        colLabels=cols,
        loc="center",
        cellLoc="center",
        colLoc="center",
        colWidths=[0.22, 0.16, 0.14, 0.18, 0.24],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.0)
    tbl.scale(1, 1.25)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_linewidth(0.3)
        cell.set_edgecolor("#B0B0B0")
        if r == 0:
            cell.set_facecolor("#F2F2F2")
            cell.set_text_props(fontweight="bold")
        elif r == 1:
            cell.set_facecolor("#FFF7F2")
        else:
            cell.set_facecolor("white")

    pdf_path = fig_dir / "S_SCENIC_FOSL1_PLAUR_edge_table.pdf"
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    manifest = pd.DataFrame(
        [
            {
                "figure": "Supplementary SCENIC FOSL1-PLAUR inferred edge table",
                "pdf": str(pdf_path),
                "source_data_dir": str(source_dir),
                "FOSL1_PLAUR_status": "found",
                "FOSL1_PLAUR_importance": float(edge["FOSL1_importance"]),
                "FOSL1_rank_among_PLAUR_upstream_TF": int(edge["FOSL1_rank_among_PLAUR_upstream_TF"]),
                "interpretation": "weak computational nomination; not causal proof",
            }
        ]
    )
    manifest.to_csv(record_dir / "Supp_SCENIC_FOSL1_PLAUR_edge_manifest.csv", index=False)
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
