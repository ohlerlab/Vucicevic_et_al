# Track Plots

`track_plot.R` and`heatmap_track_plost.R` are the same script but with different annotation output

This scripts generates heatmap track plots for visualizing genomic data. It uses various libraries including tidyverse, rtracklayer, AnnotationHub, Gviz, Rsamtools, and biomaRt. The script is divided into several sections:

1. **Loading results**: Reads in various CSV and TXT files containing genomic data.
2. **Colors**: Defines color schemes for the plots.
3. **Annotation tracks**: Imports bigwig files and creates DataTrack objects for different genomic features.
4. **Fixed parameters**: Sets up fixed parameters like genome information and annotation tracks.
5. **Regions**: Processes and visualizes data for genomic regions.
6. **Windows**: Processes and visualizes data for genomic windows.
7. **Single bins**: Processes and visualizes data for single bins.
8. **Actual plots**: Generates and saves the final plots as PDF files.

The script uses Gviz to create various DataTrack objects for different genomic features and combines them into plots. The plots are saved as PDF files, showing different genomic features and their significance.

## Visualization Scripts for Phox2B Project

This script performs the following tasks:

1. Load necessary libraries and set working directory.
2. Load various datasets for VPR and KRAB constructs, including guides, windows, bins, and regions.
3. Define plotting functions and themes for visualizations.
4. Generate volcano plots for VPR and KRAB regions, windows, and guides.
5. Create cross-tabulation plots for regions, windows, and guides.
6. Identify and plot shared bins and guides between VPR and KRAB constructs.
7. Save all generated plots to PDF files.

The script uses the following libraries:
- **tidyverse**: for data manipulation and plotting.
- **GGally**: for creating cross-tabulation plots.

The script generates the following output files:
- `region_images.pdf`: contains plots for VPR and KRAB regions.
- `window_images.pdf`: contains plots for VPR and KRAB windows.
- `significant_bins_both.pdf`: contains plots for shared bins between VPR and KRAB constructs.
- `guide_images.pdf`: contains plots for VPR and KRAB guides.
- `significant_GUIDE.pdf`: contains plots for shared guides between VPR and KRAB constructs.
