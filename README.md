# dmd-rnaseq-downstream-analysis

A scientifically cautious, reproducible transcriptomics explorer focused on **descriptive, exploratory, TPM-based** analysis.

## Live Demo
- https://dmd-rnaseq-downstream-analysis-eipcctanb8evu3u57frpva.streamlit.app

## Original Pipeline
- https://github.com/berfinida/pipeline

## Dataset Provenance
- Source dataset: `GSE156496 / SRP278118`
- Original raw data: SRA FASTQ files
- This repository uses processed TSV outputs from the upstream pipeline.

## Input Files
- `results/matrix/expression_matrix.tsv`
- `results/matrix/dmd_vs_wt_summary.tsv`
- `metadata/samples.tsv`
- Optional: `annotation/transcript_to_gene_symbol.tsv`

## Interactive Upload Mode
The dashboard supports optional upload of `expression_matrix.tsv`.

Required columns:
- `Name`, `DMD1`, `DMD2`, `WT1`, `WT2`

For uploaded matrices, the app recomputes:
- `DMD_mean_TPM`
- `WT_mean_TPM`
- `log2FC = log2((DMD_mean_TPM + 1)/(WT_mean_TPM + 1))`

If upload is invalid, the app shows a clear validation error and keeps default data.

## Optional Annotation Support
If `annotation/transcript_to_gene_symbol.tsv` exists with columns:
- `transcript_id`
- `gene_symbol`

the app merges available gene symbols for display. If absent/invalid, Ensembl transcript IDs remain in use.

## App Features
- Data provenance panel
- Data quality/context metrics
- Transcript-level explorer
- Top transcript explorer with multiple sort modes
- Exportable TSV outputs:
  - filtered table
  - top upregulated
  - top downregulated
  - current descriptive summary
- Volcano-like exploratory plot
- Heatmap explorer
- PCA sample-level exploratory view
- Static figure gallery
- Reproducibility and limitations sections

## Exportable Summary Report
Script:
- `downstream_analysis/export_report_summary.py`

Output:
- `outputs/report_summary.md`

Report includes dataset context, transcript/sample counts, top 10 up/down transcripts, limitations, and figure paths.

## Metadata
Sample metadata is stored in:
- `metadata/samples.tsv`

The app uses this file for group labels where available.

## Future Pathway Enrichment
Pathway enrichment is **not currently implemented**.

See plan:
- `docs/pathway_enrichment_plan.md`

## Reproducibility
```bash
pip install -r requirements.txt
streamlit run app.py
```

## Scientific Limitations
- Small sample size
- TPM-based descriptive comparison
- No p-values/FDR
- No formal differential expression statistics in this repo
- No pathway enrichment results currently implemented
- Not for clinical interpretation

## Documentation
- [docs/data_provenance.md](docs/data_provenance.md)
- [docs/scientific_limitations.md](docs/scientific_limitations.md)
- [docs/dashboard_usage.md](docs/dashboard_usage.md)
- [docs/pathway_enrichment_plan.md](docs/pathway_enrichment_plan.md)
- [docs/annotation.md](docs/annotation.md)
