# Data Provenance

## Source Dataset
- GEO: `GSE156496`
- SRA: `SRP278118`

## Raw Data Origin
Original reads are SRA FASTQ files processed upstream.

## Upstream Pipeline Repository
- https://github.com/berfinida/pipeline

## Downstream Scope In This Repository
This repository uses processed matrix outputs for TPM-based descriptive exploratory analysis.

Primary files used:
- `results/matrix/expression_matrix.tsv`
- `results/matrix/dmd_vs_wt_summary.tsv`
- `metadata/samples.tsv`

Optional annotation file:
- `annotation/transcript_to_gene_symbol.tsv`
