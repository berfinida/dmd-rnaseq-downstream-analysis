# Annotation Input

Gene-symbol annotation is optional and requires a real mapping file at:

- `annotation/transcript_to_gene_symbol.tsv`

Expected tab-separated columns:

```text
transcript_id	gene_symbol
```

Notes:
- `transcript_id` must match transcript IDs in the summary `Name` column.
- If this file is missing, the app keeps Ensembl transcript IDs.
- No gene symbols are invented.
