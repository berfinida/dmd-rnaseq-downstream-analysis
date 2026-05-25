# Data Provenance

## Source Dataset
- GEO: `GSE156496`
- SRA: `SRP278118`

## Raw Data Origin
The original sequencing inputs are SRA FASTQ files (upstream layer).

## Upstream Pipeline
Raw FASTQ processing and matrix generation are performed in a separate repository:
- https://github.com/berfinida/pipeline

## Downstream Repository Scope
This repository consumes processed TSV outputs and provides TPM-based descriptive exploratory visualization.

Files used directly by the dashboard:
- `results/matrix/expression_matrix.tsv`
- `results/matrix/dmd_vs_wt_summary.tsv`
