# Streamlit App

Optional interactive dashboard for the `dmd-rnaseq-downstream-analysis` project.

## Install And Run
```bash
pip install -r requirements.txt
streamlit run app.py
```

## Dashboard Features
- Wide layout with sidebar navigation
- Metric cards for transcript/sample summary
- Transcript search and top-transcript selector
- Interactive Plotly charts:
  - volcano-like scatter
  - top upregulated/downregulated bars
  - selected transcript TPM charts
- Figure gallery for generated PNG outputs

## Safety Note
This app visualizes TPM-based descriptive summaries and does not perform formal statistical differential expression testing.
