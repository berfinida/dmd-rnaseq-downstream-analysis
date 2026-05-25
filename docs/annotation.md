# Transcript Annotation

## Current Identifier Type
The summary tables use Ensembl transcript identifiers from the `Name` column.

## Optional Annotation Mapping
If an external file exists at:

- `annotation/transcript_to_gene_symbol.tsv`

with required columns:

- `transcript_id`
- `gene_symbol`

then you can run:

```bash
python downstream_analysis/annotate_transcripts.py
```

and produce:

- `results/matrix/dmd_vs_wt_summary_annotated.tsv`

## Matching Logic
- `transcript_id` is matched to summary `Name`.
- Only real mappings are merged.
- No gene symbols are invented.

## If Mapping Is Unavailable
When `annotation/transcript_to_gene_symbol.tsv` is missing:
- the script prints a clear message,
- no annotated output is created,
- Ensembl transcript IDs remain the primary identifiers.

## Scientific Scope Note
This annotation step only adds label metadata when available. It does not perform statistical testing or change the exploratory TPM-based nature of the analysis.
