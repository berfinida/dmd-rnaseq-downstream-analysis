#!/usr/bin/env python3
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    input_path = project_root / "results" / "matrix" / "dmd_vs_wt_summary.tsv"
    output_dir = project_root / "downstream_analysis" / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "volcano_like_plot.png"

    df = pd.read_csv(input_path, sep="\t")
    required_cols = ["Name", "DMD_mean_TPM", "WT_mean_TPM", "log2FC"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    plot_df = df[required_cols].copy()
    plot_df["abs_log2FC"] = plot_df["log2FC"].abs()
    plot_df["highlight"] = plot_df["abs_log2FC"] >= 2

    plt.style.use("seaborn-v0_8-whitegrid")
    fig, ax = plt.subplots(figsize=(10, 7), dpi=150)

    base = plot_df[~plot_df["highlight"]]
    high = plot_df[plot_df["highlight"]]

    ax.scatter(
        base["log2FC"],
        base["abs_log2FC"],
        s=18,
        alpha=0.5,
        color="#4C78A8",
        edgecolors="none",
        label="|log2FC| < 2",
    )
    ax.scatter(
        high["log2FC"],
        high["abs_log2FC"],
        s=26,
        alpha=0.85,
        color="#E45756",
        edgecolors="white",
        linewidths=0.3,
        label="|log2FC| >= 2",
    )

    ax.axvline(0, color="#6B6B6B", linestyle="--", linewidth=1)
    ax.set_xlabel("log2FC")
    ax.set_ylabel("abs(log2FC)")
    ax.set_title("Descriptive volcano-like plot (TPM-based, no p-values)")
    ax.legend(frameon=True)

    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
