# Downstream Analysis Module

Exploratory transcriptomics visualization scripts for DMD vs WT RNA-seq TPM outputs.

## Scripts
- `plot_volcano_like.py`: descriptive volcano-like scatter (`log2FC` vs `abs(log2FC)`)
- `plot_heatmap.py`: top-transcript heatmap from `log2(TPM + 1)` values
- `plot_pca.py`: PCA projection of samples (`DMD1`, `DMD2`, `WT1`, `WT2`)
- `annotate_transcripts.py`: optional transcript-to-gene-symbol merge helper

## Annotation Behavior
If `annotation/transcript_to_gene_symbol.tsv` exists with columns
`transcript_id` and `gene_symbol`, annotation output is written to:

- `results/matrix/dmd_vs_wt_summary_annotated.tsv`

If the mapping file is missing, the script exits without inventing annotations.

## Run
```bash
python downstream_analysis/annotate_transcripts.py
python downstream_analysis/plot_volcano_like.py
python downstream_analysis/plot_heatmap.py
python downstream_analysis/plot_pca.py
```

## Outputs
- `downstream_analysis/figures/volcano_like_plot.png`
- `downstream_analysis/figures/top_transcripts_heatmap.png`
- `downstream_analysis/figures/pca_plot.png`
