# Figure Interpretation

## Volcano-like Plot
Shows transcript-level `log2FC` versus `abs(log2FC)` to highlight large-magnitude contrasts between DMD and WT in a descriptive way.

## Heatmap
Displays top transcripts by absolute `log2FC` with row-scaled `log2(TPM + 1)` values across samples (`DMD1`, `DMD2`, `WT1`, `WT2`) to visualize relative expression patterns.

## PCA Plot
Projects samples into low-dimensional space from transcript features after `log2(TPM + 1)` transformation, enabling qualitative assessment of sample-level similarity.

## Summary Panel
Combines volcano-like view, PCA, and top up/down transcript bars into one communication-focused visual overview.

## Limitations
These figures are exploratory and descriptive. They are not inferential tests and should not be interpreted as formal differential expression significance.
