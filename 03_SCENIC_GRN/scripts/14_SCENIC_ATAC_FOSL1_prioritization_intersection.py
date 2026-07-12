#!/usr/bin/env python

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


AP1_FOS_TFS = {"FOSL1", "FOSL2", "FOSB", "FOS", "JUN", "JUNB", "JUND", "ATF3", "BATF"}
MESV_TOP_N = 20


def main() -> None:
    scenic_dir = Path(__file__).resolve().parents[1]
    project_root = scenic_dir.parents[2]
    atac_dir = project_root / "09_GSE240822_snATAC_最终分析"
    scenic_outputs = scenic_dir / "outputs"

    fig_dir = scenic_dir / "figures" / "补充图"
    source_dir = scenic_dir / "figures" / "source_data" / "FOSL1_prioritization_intersection"
    record_dir = scenic_dir / "figures" / "记录"
    for path in [fig_dir, source_dir, record_dir]:
        path.mkdir(parents=True, exist_ok=True)

    final_atac_root = (
        scenic_dir.parents[2]
        / "图片表格"
        / "05_发育时间_TF_ATAC验证"
        / "03_GSE240822_snATAC_最终图表"
        / "表格与source_data"
    )
    atac_motif_block = pd.read_csv(final_atac_root / "X-B_Subtype3_motif_block_heatmap_source.csv")
    atac_mes_motif_tfs = set(
        atac_motif_block.loc[atac_motif_block["block"].eq("MES-lineage"), "motif"]
        .astype(str)
        .str.replace(r"[_\\.].*$", "", regex=True)
    )

    plaur = pd.read_csv(scenic_outputs / "scenic_phase0_5_PLAUR_upstream_tf_ranking.csv")
    rss = pd.read_csv(scenic_outputs / "scenic_rss_all_regulons_ranked_by_subtype.csv")
    mesv = rss.loc[rss["subtype"].eq("MES-V"), ["tf", "regulon", "rank", "rss"]].rename(
        columns={"rank": "MESV_RSS_rank", "rss": "MESV_RSS"}
    )
    merged = plaur.merge(mesv, left_on="TF", right_on="tf", how="left")
    merged["is_PLAUR_upstream_TF"] = True
    merged["in_MESV_top20_RSS"] = merged["MESV_RSS_rank"].le(MESV_TOP_N)
    merged["is_AP1_FOS_family"] = merged["TF"].isin(AP1_FOS_TFS)
    merged["has_ATAC_MESlineage_characteristic_motif"] = merged["TF"].isin(atac_mes_motif_tfs)
    merged["passes_all_three_filters"] = (
        merged["is_PLAUR_upstream_TF"]
        & merged["in_MESV_top20_RSS"].fillna(False)
        & merged["has_ATAC_MESlineage_characteristic_motif"]
    )

    out_cols = [
        "TF",
        "rank",
        "importance",
        "MESV_RSS_rank",
        "MESV_RSS",
        "is_PLAUR_upstream_TF",
        "in_MESV_top20_RSS",
        "is_AP1_FOS_family",
        "has_ATAC_MESlineage_characteristic_motif",
        "passes_all_three_filters",
    ]
    candidates = merged.loc[merged["MESV_RSS_rank"].notna(), out_cols].copy()
    candidates = candidates.sort_values(
        ["passes_all_three_filters", "in_MESV_top20_RSS", "MESV_RSS_rank", "rank"],
        ascending=[False, False, True, True],
    )
    candidates.to_csv(source_dir / "SCENIC_PLAUR_upstream_MESV_RSS_AP1_intersection.csv", index=False)

    top20_overlap = candidates.loc[candidates["in_MESV_top20_RSS"]].copy()
    top20_overlap = top20_overlap.sort_values(["MESV_RSS_rank", "rank"])
    top20_overlap.to_csv(source_dir / "PLAUR_upstream_and_MESV_top20_overlap.csv", index=False)
    final_hit = candidates.loc[candidates["passes_all_three_filters"]].copy()
    final_hit.to_csv(source_dir / "FOSL1_three_filter_final_hit.csv", index=False)

    # ATAC evidence is motif-level, not direct TF binding evidence.
    atac_rank_path = final_atac_root / "X-C_AP1_family_rank_by_subtype.csv"
    atac_anchor_path = final_atac_root / "X-D_PLAUR_locus_anchor_coordinates.csv"
    atac_block_path = final_atac_root / "X-B_Subtype3_motif_block_heatmap_source.csv"
    if atac_rank_path.exists():
        pd.read_csv(atac_rank_path).to_csv(source_dir / "ATAC_AP1_family_rank_by_subtype.csv", index=False)
    if atac_anchor_path.exists():
        pd.read_csv(atac_anchor_path).to_csv(source_dir / "ATAC_PLAUR_locus_anchor_coordinates.csv", index=False)
    if atac_block_path.exists():
        pd.read_csv(atac_block_path).to_csv(source_dir / "ATAC_Subtype3_motif_block_heatmap_source.csv", index=False)

    summary = pd.DataFrame(
        [
            {
                "filter_step": "1. SCENIC PLAUR upstream TFs",
                "criterion": "TF has inferred GRN edge to PLAUR",
                "n_passing": int(plaur["TF"].nunique()),
                "passing_TFs": "see source table",
            },
            {
                "filter_step": "2. MES-V regulon specificity",
                "criterion": f"TF is among MES-V top {MESV_TOP_N} RSS-ranked regulons",
                "n_passing": int(top20_overlap["TF"].nunique()),
                "passing_TFs": ", ".join(top20_overlap["TF"].tolist()),
            },
            {
                "filter_step": "3. ATAC MES-lineage motif context",
                "criterion": "TF has a corresponding motif in the ATAC MES-lineage characteristic motif block",
                "n_passing": int(final_hit["TF"].nunique()),
                "passing_TFs": ", ".join(final_hit["TF"].tolist()),
            },
        ]
    )
    summary.to_csv(source_dir / "FOSL1_prioritization_filter_summary.csv", index=False)

    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "sans-serif"]
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42

    fig = plt.figure(figsize=(7.8, 3.9))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    ax.set_title("Three-filter prioritization nominates FOSL1 for PLAUR validation", fontsize=11, fontweight="bold", pad=8)

    x = [0.17, 0.50, 0.83]
    labels = [
        "SCENIC PLAUR\nupstream TFs\nn = 86",
        "also MES-V top20\nRSS regulons\nn = 6",
        "also ATAC MES-lineage\ncharacteristic motif\nn = 1",
    ]
    colors = ["#F2F2F2", "#EAF3EC", "#FFF2E8"]
    for i, (xx, lab, col) in enumerate(zip(x, labels, colors)):
        box = matplotlib.patches.FancyBboxPatch(
            (xx - 0.13, 0.58),
            0.26,
            0.23,
            boxstyle="round,pad=0.02,rounding_size=0.02",
            facecolor=col,
            edgecolor="#888888",
            linewidth=0.8,
        )
        ax.add_patch(box)
        ax.text(xx, 0.695, lab, ha="center", va="center", fontsize=8.5)
        if i < 2:
            ax.annotate("", xy=(x[i + 1] - 0.15, 0.695), xytext=(xx + 0.15, 0.695),
                        arrowprops=dict(arrowstyle="->", lw=1.0, color="#555555"))

    ax.text(0.50, 0.46, "PLAUR upstream ∩ MES-V top20:", ha="center", fontsize=8.5, fontweight="bold")
    ax.text(0.50, 0.39, ", ".join(top20_overlap["TF"].tolist()), ha="center", fontsize=8.0)
    ax.text(0.50, 0.27, "Final intersection: FOSL1", ha="center", fontsize=11.5, fontweight="bold", color="#BC3C29")
    ax.text(
        0.50,
        0.15,
        "ATAC provides MES-lineage motif/chromatin context, not direct FOSL1 binding evidence.",
        ha="center",
        fontsize=7.6,
        color="#444444",
    )

    pdf_path = fig_dir / "S_SCENIC_ATAC_FOSL1_prioritization_intersection.pdf"
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)

    manifest = pd.DataFrame(
        [
            {
                "figure": "Supplementary SCENIC-ATAC FOSL1 prioritization intersection",
                "pdf": str(pdf_path),
                "source_data_dir": str(source_dir),
                "filter1_PLAUR_upstream_TFs": int(plaur["TF"].nunique()),
                "filter2_overlap_with_MESV_top20": int(top20_overlap["TF"].nunique()),
                "filter2_TFs": ", ".join(top20_overlap["TF"].tolist()),
                "filter3_final_ATAC_MESlineage_motif_hit": ", ".join(final_hit["TF"].tolist()),
                "interpretation": "prioritization evidence; ATAC is motif/chromatin context and not direct TF binding proof",
            }
        ]
    )
    manifest.to_csv(record_dir / "Supp_SCENIC_ATAC_FOSL1_prioritization_intersection_manifest.csv", index=False)
    print(manifest.to_string(index=False))


if __name__ == "__main__":
    main()
