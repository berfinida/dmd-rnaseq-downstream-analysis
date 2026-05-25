# Dashboard Usage

## Run
```bash
pip install -r requirements.txt
streamlit run app.py
```

## Key Features
- Optional upload mode for `expression_matrix.tsv`
- Transcript-level explorer by Ensembl transcript ID
- Top transcript sorting and TSV export
- Volcano-like, heatmap, and PCA exploratory views
- Optional annotation support (if mapping file exists)

## Upload Validation Requirements
Uploaded matrix must include:
- `Name`
- `DMD1`
- `DMD2`
- `WT1`
- `WT2`

## Scientific Scope
Dashboard outputs are descriptive and exploratory (TPM-based), not formal inferential differential expression analysis.
