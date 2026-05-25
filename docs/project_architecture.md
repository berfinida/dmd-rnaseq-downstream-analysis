# Project Architecture

## Overview
- The original Nextflow RNA-seq pipeline performs upstream processing and produces matrix outputs.
- This downstream repository consumes processed TSV outputs from that upstream workflow.
- Python scripts in `downstream_analysis/` generate static figures.
- `app.py` provides interactive exploration via Streamlit.

## Data/Processing Flow
```text
Nextflow RNA-seq Pipeline
        ↓
expression_matrix.tsv + dmd_vs_wt_summary.tsv
        ↓
Python downstream scripts
        ↓
Static figures + Streamlit dashboard
```

## Repository Scope
This repository is downstream-only and exploratory. It visualizes TPM-based summaries and does not perform formal statistical differential expression testing.
