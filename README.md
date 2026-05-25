# DMD RNA-seq Downstream Exploratory Analysis

This repository focuses on downstream, descriptive exploration of an existing DMD vs WT transcript expression output.

## Source Context
Input expression outputs are based on an existing Nextflow RNA-seq pipeline run and are used here only for post hoc exploratory visualization.

## Input Data
- `results/matrix/expression_matrix.tsv`
- `results/matrix/dmd_vs_wt_summary.tsv`

## Downstream Exploratory Analysis
This project includes two descriptive visualizations:

1. Volcano-like plot
- x-axis: `log2FC`
- y-axis: `abs(log2FC)`
- Highlights transcripts where `abs(log2FC) >= 2`

2. Top transcript heatmap
- Selects top 20 transcripts by absolute `log2FC`
- Uses sample TPM columns: `DMD1`, `DMD2`, `WT1`, `WT2`
- Applies `log2(TPM + 1)` then row-wise z-score scaling

## Run Commands
```bash
python downstream_analysis/plot_volcano_like.py
python downstream_analysis/plot_heatmap.py
```

## Generated Figures
### Volcano-like Plot
![Descriptive volcano-like plot](downstream_analysis/figures/volcano_like_plot.png)

### Top Transcript Heatmap
![Top transcripts heatmap](downstream_analysis/figures/top_transcripts_heatmap.png)

## Important Note
This is descriptive exploratory analysis and not a formal statistical differential expression workflow.
No p-values or FDR-adjusted significance estimates are computed in this repository.
