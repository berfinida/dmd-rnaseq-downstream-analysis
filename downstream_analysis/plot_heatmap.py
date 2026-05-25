#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    means = df.mean(axis=1)
    stds = df.std(axis=1)
    stds_safe = stds.mask(stds == 0, 1.0)
    return df.sub(means, axis=0).div(stds_safe, axis=0)


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    matrix_path = project_root / "results" / "matrix" / "expression_matrix.tsv"
    summary_path = project_root / "results" / "matrix" / "dmd_vs_wt_summary.tsv"
    output_dir = project_root / "downstream_analysis" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "top_transcripts_heatmap.png"

    sample_cols = ["DMD1", "DMD2", "WT1", "WT2"]

    expr = pd.read_csv(matrix_path, sep="\t")
    summary = pd.read_csv(summary_path, sep="\t")

    if "Name" not in expr.columns:
        raise ValueError("Missing 'Name' in expression matrix")
    if "Name" not in summary.columns or "log2FC" not in summary.columns:
        raise ValueError("Missing 'Name' or 'log2FC' in summary file")

    missing_samples = [c for c in sample_cols if c not in expr.columns]
    if missing_samples:
        raise ValueError(f"Missing sample columns: {missing_samples}")

    top_names = (
        summary.assign(abs_log2FC=summary["log2FC"].abs())
        .sort_values("abs_log2FC", ascending=False)
        .head(20)["Name"]
        .tolist()
    )

    filtered = expr[expr["Name"].isin(top_names)].copy()
    filtered = filtered.set_index("Name").reindex(top_names)

    log_tpm = np.log2(filtered[sample_cols] + 1)
    scaled = zscore_rows(log_tpm).fillna(0)

    plt.style.use("seaborn-v0_8-white")
    fig, ax = plt.subplots(figsize=(8, 10), dpi=150)
    sns.heatmap(
        scaled,
        cmap="vlag",
        center=0,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"label": "Row-scaled log2(TPM + 1)"},
        ax=ax,
    )

    ax.set_title("Top 20 transcripts by |log2FC| (descriptive heatmap)")
    ax.set_xlabel("Samples")
    ax.set_ylabel("Transcripts")

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
