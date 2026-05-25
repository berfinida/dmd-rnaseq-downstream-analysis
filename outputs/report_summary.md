# Exploratory Transcriptomics Report Summary

## Dataset Context
- Source dataset: GSE156496 / SRP278118
- Groups: WT vs DMD Delta51
- This report is TPM-based descriptive exploration

## Data Dimensions
- Number of transcripts: 115424
- Number of samples: 4 (DMD1, DMD2, WT1, WT2)

## Top 10 Upregulated (by log2FC)
| Name | DMD_mean_TPM | WT_mean_TPM | log2FC |
| --- | --- | --- | --- |
| ENSMUST00000031243.15 | 69.444322 | 0.563595 | 5.493545 |
| ENSMUST00000116304.3 | 43.188089 | 0.0 | 5.465586 |
| ENSMUST00000015667.9 | 139.298699 | 2.225894 | 5.442659 |
| ENSMUST00000146468.4 | 101.765183 | 1.456659 | 5.38651 |
| ENSMUST00000207870.2 | 39.788249 | 0.0 | 5.350082 |
| ENSMUST00000134453.2 | 33.905514 | 0.0 | 5.125383 |
| ENSMUST00000226585.2 | 31.560421 | 0.0 | 5.025047 |
| ENSMUST00000212806.2 | 260.295989 | 7.917231 | 4.872945 |
| ENSMUST00000117373.8 | 57.748649 | 1.178707 | 4.753012 |
| ENSMUST00000143511.2 | 25.448803 | 0.0 | 4.725131 |

## Top 10 Downregulated (by log2FC)
| Name | DMD_mean_TPM | WT_mean_TPM | log2FC |
| --- | --- | --- | --- |
| ENSMUST00000229280.2 | 0.0 | 55.903163 | -5.830437 |
| ENSMUST00000211494.2 | 0.0 | 36.742198 | -5.238107 |
| ENSMUST00000162550.2 | 0.0 | 25.872584 | -4.748063 |
| ENSMUST00000118324.2 | 0.0 | 22.944367 | -4.581614 |
| ENSMUST00000226671.2 | 0.0 | 14.059048 | -3.912559 |
| ENSMUST00000127425.4 | 0.0 | 13.601706 | -3.868065 |
| ENSMUST00000148569.4 | 0.0 | 12.817679 | -3.788443 |
| ENSMUST00000206012.2 | 0.0 | 12.17586 | -3.719825 |
| ENSMUST00000129170.8 | 0.0 | 11.576398 | -3.652647 |
| ENSMUST00000168856.2 | 0.0 | 11.040612 | -3.589837 |

## Figure Paths
- downstream_analysis/figures/volcano_like_plot.png
- downstream_analysis/figures/top_transcripts_heatmap.png
- downstream_analysis/figures/pca_plot.png
- downstream_analysis/figures/summary_panel.png

## Limitations
- Small sample size (2 WT, 2 DMD)
- Descriptive TPM-based comparison
- No p-values/FDR
- Not for clinical interpretation