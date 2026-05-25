#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    matrix_path = project_root / "results" / "matrix" / "expression_matrix.tsv"
    output_dir = project_root / "downstream_analysis" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "pca_plot.png"

    sample_cols = ["DMD1", "DMD2", "WT1", "WT2"]

    df = pd.read_csv(matrix_path, sep="\t")
    if "Name" not in df.columns:
        raise ValueError("Missing 'Name' column in expression_matrix.tsv")

    missing = [c for c in sample_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing sample columns: {missing}")

    expr = df[sample_cols].copy()
    log_expr = np.log2(expr + 1)
    x = log_expr.T.values

    pca = PCA(n_components=2, random_state=42)
    pcs = pca.fit_transform(x)
    explained = pca.explained_variance_ratio_ * 100

    plot_df = pd.DataFrame(
        {
            "sample": sample_cols,
            "PC1": pcs[:, 0],
            "PC2": pcs[:, 1],
            "group": ["DMD", "DMD", "WT", "WT"],
        }
    )

    colors = {"DMD": "#D1495B", "WT": "#2E86AB"}

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(8, 6), dpi=150)

    for group, gdf in plot_df.groupby("group"):
        ax.scatter(
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
            ax.text(row["PC1"] + 0.02, row["PC2"] + 0.02, row["sample"], fontsize=9)

    ax.set_xlabel(f"PC1 ({explained[0]:.2f}% variance)")
    ax.set_ylabel(f"PC2 ({explained[1]:.2f}% variance)")
    ax.set_title("PCA of samples from log2(TPM + 1) transcript matrix")
    ax.legend(title="Group")
    ax.axhline(0, color="#9B9B9B", linewidth=0.8, linestyle="--")
    ax.axvline(0, color="#9B9B9B", linewidth=0.8, linestyle="--")

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
