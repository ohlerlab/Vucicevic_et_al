# Scripts Overview

## `filtering_sig_bins_to_bed.R`

This script processes genomic data to identify significant bins based on p-values.

### Libraries:
- **tidyverse**: A collection of R packages for data manipulation and visualization.

### Steps:
1. **Load the necessary data:**
    - `num_VPR_pvalues.txt`: Contains bin, estimate, and p-value information.
    - `bin_coords.txt`: Contains chromosome, start, end, and bin coordinates.
2. **Adjust p-values** for multiple testing using the False Discovery Rate (FDR) method.
3. **Merge the coordinates data** with the p-value data based on the bin identifier.
4. **Filter the merged data** to retain only bins with an FDR-adjusted p-value less than 0.05.
5. **Select relevant columns** (chr, start, end, bin) from the filtered data.
6. **Write the filtered data** to a BED file (`sig_window_bins.bed`) without column names.

## `merging.sh`

This script uses bedtools to merge overlapping or nearby intervals in a BED file.

### Usage:
```sh
./merging.sh
```

### Requirements:
- **bedtools** must be installed and accessible in the system's PATH.
- The input BED file (`sig_window_bins.bed`) must be present in the same directory as this script.

### Description:
The script merges intervals in the input BED file (`sig_window_bins.bed`) that are within 500 base pairs of each other. It collapses the fourth column values of the merged intervals into a comma-separated list. The output is written to a new BED file (`merged_windowed_500.bed`).

### Options:
- `-d 500`: Merge intervals that are within 500 base pairs of each other.
- `-c 4`: Specify the column to operate on (fourth column).
- `-o collapse`: Collapse the values in the specified column into a comma-separated list.