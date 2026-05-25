#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def prepare_pca(expr: pd.DataFrame, sample_cols: list[str]) -> tuple[pd.DataFrame, np.ndarray]:
    log_expr = np.log2(expr[sample_cols] + 1)
    x = log_expr.T.values
    pca = PCA(n_components=2, random_state=42)
    pcs = pca.fit_transform(x)
    explained = pca.explained_variance_ratio_ * 100
    pca_df = pd.DataFrame(
        {
            "sample": sample_cols,
            "PC1": pcs[:, 0],
            "PC2": pcs[:, 1],
            "group": ["DMD", "DMD", "WT", "WT"],
        }
    )
    return pca_df, explained


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    summary_path = project_root / "results" / "matrix" / "dmd_vs_wt_summary.tsv"
    matrix_path = project_root / "results" / "matrix" / "expression_matrix.tsv"
    up_path = project_root / "results" / "matrix" / "top10_upregulated.tsv"
    down_path = project_root / "results" / "matrix" / "top10_downregulated.tsv"

    output_dir = project_root / "downstream_analysis" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "summary_panel.png"

    sample_cols = ["DMD1", "DMD2", "WT1", "WT2"]

    summary = pd.read_csv(summary_path, sep="\t")
    expr = pd.read_csv(matrix_path, sep="\t")
    up = pd.read_csv(up_path, sep="\t")
    down = pd.read_csv(down_path, sep="\t")

    if "log2FC" not in summary.columns:
        raise ValueError("Missing 'log2FC' in dmd_vs_wt_summary.tsv")

    pca_df, explained = prepare_pca(expr, sample_cols)

    plot_df = summary[["log2FC"]].copy()
    plot_df["abs_log2FC"] = plot_df["log2FC"].abs()
    plot_df["highlight"] = plot_df["abs_log2FC"] >= 2

    up_col = "log2FC" if "log2FC" in up.columns else up.columns[-1]
    down_col = "log2FC" if "log2FC" in down.columns else down.columns[-1]

    up_plot = up.sort_values(up_col, ascending=False).head(10)
    down_plot = down.sort_values(down_col, ascending=True).head(10)

    up_labels = up_plot["Name"] if "Name" in up_plot.columns else up_plot.iloc[:, 0]
    down_labels = (
        down_plot["Name"] if "Name" in down_plot.columns else down_plot.iloc[:, 0]
    )

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=150)

    ax1 = axes[0, 0]
    base = plot_df[~plot_df["highlight"]]
    high = plot_df[plot_df["highlight"]]
    ax1.scatter(base["log2FC"], base["abs_log2FC"], s=14, alpha=0.45, color="#4C78A8")
    ax1.scatter(high["log2FC"], high["abs_log2FC"], s=20, alpha=0.8, color="#E45756")
    ax1.axvline(0, color="#6B6B6B", linestyle="--", linewidth=1)
    ax1.set_title("Volcano-like (descriptive)")
    ax1.set_xlabel("log2FC")
    ax1.set_ylabel("abs(log2FC)")

    ax2 = axes[0, 1]
    colors = {"DMD": "#D1495B", "WT": "#2E86AB"}
    for group, gdf in pca_df.groupby("group"):
        ax2.scatter(
            gdf["PC1"],
            gdf["PC2"],
            s=120,
            alpha=0.9,
            color=colors[group],
            edgecolors="white",
            linewidths=0.8,
            label=group,
        )
        for _, row in gdf.iterrows():
            ax2.text(row["PC1"] + 0.02, row["PC2"] + 0.02, row["sample"], fontsize=8)
    ax2.set_title("PCA (sample-level exploratory)")
    ax2.set_xlabel(f"PC1 ({explained[0]:.1f}% var)")
    ax2.set_ylabel(f"PC2 ({explained[1]:.1f}% var)")
    ax2.legend(frameon=True, title="Group")

    ax3 = axes[1, 0]
    ax3.barh(up_labels.astype(str), up_plot[up_col], color="#D1495B", alpha=0.85)
    ax3.invert_yaxis()
    ax3.set_title("Top 10 upregulated transcripts")
    ax3.set_xlabel("log2FC")

    ax4 = axes[1, 1]
    ax4.barh(down_labels.astype(str), down_plot[down_col], color="#2E86AB", alpha=0.85)
    ax4.invert_yaxis()
    ax4.set_title("Top 10 downregulated transcripts")
    ax4.set_xlabel("log2FC")

    fig.suptitle("DMD vs WT transcriptomics summary panel (TPM-based exploratory)")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
