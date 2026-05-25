# Downstream Exploratory Analysis

## Purpose
This module provides descriptive downstream exploratory visualizations from TPM-derived transcript expression summaries in a DMD vs WT comparison.

## Input Files
- `results/matrix/expression_matrix.tsv`
- `results/matrix/dmd_vs_wt_summary.tsv`

## Scripts
- `plot_volcano_like.py`: descriptive volcano-like scatter using `log2FC` and `abs(log2FC)`
- `plot_heatmap.py`: heatmap of top 20 transcripts by absolute `log2FC` across `DMD1`, `DMD2`, `WT1`, `WT2`

## Run
```bash
python downstream_analysis/plot_volcano_like.py
python downstream_analysis/plot_heatmap.py
```

## Output Figures
- `downstream_analysis/figures/volcano_like_plot.png`
- `downstream_analysis/figures/top_transcripts_heatmap.png`

## Important Note
This analysis is descriptive and exploratory. It is not a formal statistical differential expression analysis because no p-values or FDR values were calculated.
