# Methods

## 1) Matrix Input Preparation
This repository starts from precomputed transcript-level matrices:
- `results/matrix/expression_matrix.tsv`
- `results/matrix/dmd_vs_wt_summary.tsv`

The input matrices are generated in a separate RNA-seq pipeline and imported here for downstream exploratory visualization.

## 2) TPM-based Descriptive Comparison
Sample-level expression values are TPM-based. Downstream comparisons in this repository are descriptive and exploratory.

## 3) log2 Fold-change Calculation
`log2FC` values are read from `dmd_vs_wt_summary.tsv` (generated upstream). This repository visualizes those values and does not recompute formal DE statistics.

## 4) Volcano-like Plot Generation
Script: `downstream_analysis/plot_volcano_like.py`
- Uses `log2FC` and `abs(log2FC)`
- Highlights transcripts with large-magnitude fold changes
- Produces a descriptive volcano-like scatter

## 5) Heatmap Generation
Script: `downstream_analysis/plot_heatmap.py`
- Selects top 20 transcripts by absolute `log2FC`
- Uses `DMD1`, `DMD2`, `WT1`, `WT2`
- Applies `log2(TPM + 1)`
- Applies row-wise z-score scaling for visual comparison
- Uses transcript IDs by default; uses gene symbols when optional annotation is available

## 6) PCA Generation
Script: `downstream_analysis/plot_pca.py`
- Uses `DMD1`, `DMD2`, `WT1`, `WT2`
- Applies `log2(TPM + 1)` transform
- Transposes matrix to samples-as-rows, transcripts-as-features
- Runs PCA with scikit-learn for sample-level exploratory structure

## 7) Summary Panel Generation
Script: `downstream_analysis/plot_summary_panel.py`
- Combines four visual components in one panel:
  - volcano-like plot
  - PCA scatter
  - top 10 upregulated transcript bar chart
  - top 10 downregulated transcript bar chart
