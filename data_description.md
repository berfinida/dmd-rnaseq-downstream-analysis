# Data Description

## Source Dataset
- GEO: `GSE156496`
- SRA: `SRP278118`

## Biological Groups
- WT (wild-type)
- DMD ΔEx51

## Selected Samples Used In This Downstream Repository
- `WT1`
- `WT2`
- `DMD1`
- `DMD2`

## Processing Context
The original raw FASTQ files were processed in a separate upstream Nextflow RNA-seq pipeline repository.

- Upstream pipeline: https://github.com/berfinida/pipeline

This repository is downstream-only and uses processed matrix outputs (for example expression and summary TSV files), not raw sequencing files.
