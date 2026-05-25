# Dashboard Usage

## Start the App
```bash
pip install -r requirements.txt
streamlit run app.py
```

## Dashboard Sections
- Overview & Data Provenance
- Data Quality & Context
- Transcript-Level Explorer
- Top Transcript Explorer
- PCA (Sample-Level Exploration)
- Figure Gallery
- Reproducibility
- Limitations

## Input Files
- `results/matrix/expression_matrix.tsv`
- `results/matrix/dmd_vs_wt_summary.tsv`

## Scientific Scope Note
This dashboard visualizes TPM-based descriptive summaries and does not perform formal statistical differential expression testing.
