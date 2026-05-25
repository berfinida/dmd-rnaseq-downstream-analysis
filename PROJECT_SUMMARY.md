# Project Summary: dmd-rnaseq-downstream-analysis

## Biological Question
This project explores transcript-level expression differences between DMD and WT muscle RNA-seq samples using TPM-derived matrices from an upstream quantification pipeline.

## Dataset Context
- Condition groups: DMD vs WT
- Sample columns used downstream: `DMD1`, `DMD2`, `WT1`, `WT2`
- Inputs are transcript-level tables exported from a separate RNA-seq processing workflow.

## What Each Figure Shows
- `volcano_like_plot.png`: descriptive `log2FC` vs `abs(log2FC)` view to highlight large-magnitude fold changes.
- `top_transcripts_heatmap.png`: row-scaled heatmap of top 20 transcripts by absolute `log2FC`, based on `log2(TPM + 1)`.
- `pca_plot.png`: sample-level PCA projection from transcript features after `log2(TPM + 1)` transform.
- `summary_panel.png`: integrated dashboard combining volcano-like view, PCA, and top 10 up/down transcript bar charts.

## Why This Is Exploratory
All outputs here are TPM-based descriptive visualizations for pattern discovery and communication. This repository does not perform formal hypothesis testing or estimate statistical significance.

## Limitations
- Small sample size (2 DMD and 2 WT) limits robust inference.
- No explicit control for batch effects or covariates.
- Fold-change summaries are sensitive to low-expression transcripts.
- Transcript-to-gene symbol annotation is optional and depends on external mapping availability.

## Connection To Original RNA-seq Pipeline
This repository is downstream-only and uses matrix outputs generated from the upstream pipeline:
- https://github.com/berfinida/pipeline
