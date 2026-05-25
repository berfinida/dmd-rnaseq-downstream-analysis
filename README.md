# dmd-rnaseq-downstream-analysis

![Python](https://img.shields.io/badge/Python-3.x-blue)
![Streamlit](https://img.shields.io/badge/Streamlit-Dashboard-FF4B4B)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Status](https://img.shields.io/badge/Status-Portfolio%20Project-success)

A reproducible computational biology project for **TPM-based descriptive exploratory** transcriptomics visualization in a DMD vs WT context.

## Live Demo
- Streamlit Live Demo: https://dmd-rnaseq-downstream-analysis-eipcctanb8evu3u57frpva.streamlit.app

## Dataset Provenance
- Source dataset accession: `GSE156496` / `SRP278118`
- Original raw data: SRA FASTQ files
- Upstream processing pipeline (separate repo): https://github.com/berfinida/pipeline
- This repository consumes processed outputs:
  - `results/matrix/expression_matrix.tsv`
  - `results/matrix/dmd_vs_wt_summary.tsv`

## Scientific Scope
This repository is for descriptive exploratory analysis.
- It does **not** perform formal differential expression statistics.
- It does **not** compute p-values or FDR.
- It does **not** provide clinical interpretation.

## App Features
Main app: `app.py`

- Data provenance and project context panel
- Transcript-level explorer by Ensembl transcript ID
- TPM display across `DMD1`, `DMD2`, `WT1`, `WT2`
- Summary display for `DMD_mean_TPM`, `WT_mean_TPM`, and `log2FC`
- Top transcript explorer with multiple real sort modes
- Download button for filtered top transcript tables
- Sample-level PCA on `log2(TPM + 1)`
- Figure gallery for generated static plots
- Reproducibility and limitations sections inside the dashboard

## Reproducibility
```bash
pip install -r requirements.txt
streamlit run app.py
```

## Figures
- `downstream_analysis/figures/volcano_like_plot.png`
- `downstream_analysis/figures/top_transcripts_heatmap.png`
- `downstream_analysis/figures/pca_plot.png`
- `downstream_analysis/figures/summary_panel.png`

## Documentation
- [docs/data_provenance.md](docs/data_provenance.md)
- [docs/dashboard_usage.md](docs/dashboard_usage.md)
- [docs/scientific_limitations.md](docs/scientific_limitations.md)
- [docs/project_architecture.md](docs/project_architecture.md)
- [docs/figure_interpretation.md](docs/figure_interpretation.md)

## Annotation Note
If gene-symbol annotation is unavailable, Ensembl transcript IDs are shown. Gene-symbol labeling requires an external mapping file (for example `annotation/transcript_to_gene_symbol.tsv`).

## Not a Formal DE Workflow
This project remains a data-grounded exploratory visualization layer and should not be interpreted as formal inferential differential expression analysis.
