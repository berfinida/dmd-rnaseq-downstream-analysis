#!/usr/bin/env python3
from pathlib import Path

import pandas as pd


def to_markdown_table(df: pd.DataFrame) -> str:
    cols = df.columns.tolist()
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    rows = []
    for _, row in df.iterrows():
        rows.append("| " + " | ".join(str(row[c]) for c in cols) + " |")
    return "\n".join([header, sep] + rows)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    matrix_path = root / "results" / "matrix" / "expression_matrix.tsv"
    summary_path = root / "results" / "matrix" / "dmd_vs_wt_summary.tsv"
    out_dir = root / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "report_summary.md"

    expr = pd.read_csv(matrix_path, sep="\t")
    summary = pd.read_csv(summary_path, sep="\t")

    top_up = summary.sort_values("log2FC", ascending=False).head(10)
    top_down = summary.sort_values("log2FC", ascending=True).head(10)

    lines = []
    lines.append("# Exploratory Transcriptomics Report Summary")
    lines.append("")
    lines.append("## Dataset Context")
    lines.append("- Source dataset: GSE156496 / SRP278118")
    lines.append("- Groups: WT vs DMD Delta51")
    lines.append("- This report is TPM-based descriptive exploration")
    lines.append("")
    lines.append("## Data Dimensions")
    lines.append(f"- Number of transcripts: {expr.shape[0]}")
    lines.append("- Number of samples: 4 (DMD1, DMD2, WT1, WT2)")
    lines.append("")
    lines.append("## Top 10 Upregulated (by log2FC)")
    lines.append(to_markdown_table(top_up[["Name", "DMD_mean_TPM", "WT_mean_TPM", "log2FC"]]))
    lines.append("")
    lines.append("## Top 10 Downregulated (by log2FC)")
    lines.append(to_markdown_table(top_down[["Name", "DMD_mean_TPM", "WT_mean_TPM", "log2FC"]]))
    lines.append("")
    lines.append("## Figure Paths")
    lines.append("- downstream_analysis/figures/volcano_like_plot.png")
    lines.append("- downstream_analysis/figures/top_transcripts_heatmap.png")
    lines.append("- downstream_analysis/figures/pca_plot.png")
    lines.append("- downstream_analysis/figures/summary_panel.png")
    lines.append("")
    lines.append("## Limitations")
    lines.append("- Small sample size (2 WT, 2 DMD)")
    lines.append("- Descriptive TPM-based comparison")
    lines.append("- No p-values/FDR")
    lines.append("- Not for clinical interpretation")

    out_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
