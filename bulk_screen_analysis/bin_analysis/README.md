# Linear Mixed Model Analysis on VPR and KRAB Datasets

This script performs linear mixed model analysis on VPR and KRAB datasets. It reads guide and code information, normalizes the data, and calculates p-values for each bin using parallel processing.

## Required Libraries

- `tidyverse`: for data manipulation and visualization
- `lme4`: for linear mixed-effects models
- `plyr`: for data manipulation
- `parallel`: for parallel processing

## Input Files

- `../tab/phox2bLocus3.tsv`: contains guide information with columns `NAME` and `GUIDE`
- `../codes.txt`: contains code information with columns `GUIDE` and `CODE`
- `../named_counts/KRAB_table.txt`: contains KRAB count data
- `../named_counts/VPR_table.txt`: contains VPR count data

## Output Files

- `num_VPR_pvalues.txt`: contains p-values for VPR data
- `num_KRAB_pvalues.txt`: contains p-values for KRAB data
- `speed_VPR_pvalues.txt`: contains early vs late classification for VPR data (commented out)
- `speed_KRAB_pvalues.txt`: contains early vs late classification for KRAB data (commented out)

## Main Steps

1. Read and preprocess guide and code data.
2. Normalize KRAB and VPR count data.
3. Perform linear mixed model analysis to calculate p-values for each bin.
4. Save the results to output files.