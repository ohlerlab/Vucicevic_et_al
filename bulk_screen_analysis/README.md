# Bulk Screen Analysis

This folder contains all scripts necessary to perform analysis on CRISPRi and CRISPRa screens.

1. `/tab` contains all metadata for sgRNAs and bins for the targeted chromosomal region.
2. `/MLM` contains all functions necessary to run the MLM 
3. `/bin_analysis` contains the scripts used to perform linear mixed model analysis on VPR and KRAB datasets. It reads guide and code information, normalizes the data, and calculates p-values for each bin using parallel processing.
4. `/region_generation` contains scripts used to identify significant bins from the previous step, as well as a bash script to merge overlapping or nearby intervals in a BED file.
5. `/visualization` includes scripts used for most of the plots produced by the bulk screen analysis.

After identifying significant bins and merging regions, most descriptive and exploratory analyses are performed in `new_guide_selection_visualization.Rmd`. In this file we score each of the CaRE, output final results as well as intermediate files used for plotting scripts in the `/visualization` folder

## new_guide_selection_visualization.Rmd

This file performs the following tasks:

1. **Setup**: Loads necessary libraries and sets global options.
2. **Data Loading**: Loads and processes raw region data and results from various analyses (VPR, KRAB).
3. **Data Preparation**: Merges and processes data to integrate all information for further analysis.
4. **Statistical Analysis**: Runs linear models and other statistical analyses on the data.
5. **Visualization**: Generates various plots to visualize the results, including region and bin results.
6. **Guide Selection**: Identifies and selects the top guides based on specific criteria.
7. **Output**: Saves the final selection of guides and results to files for further use.

### Output Files

#### MaGeCK Input Files

- `mageck/input_mageck_vpr.txt`
- `mageck/input_mageck_krab.txt`

#### Region and Bin Results

- `regions_results.txt`
- `regions_bins_results.txt`
- `regions_bins_guides_results.txt`

#### Files for Visualization Scripts

- `visualization/VPR_guide_results_plot.csv`
- `visualization/KRAB_guide_results_plot.csv`
- `visualization/VPR_window_results_plot.csv`
- `visualization/KRAB_window_results_plot.csv`
- `visualization/VPR_region_results_plot.csv`
- `visualization/KRAB_region_results_plot.csv`
