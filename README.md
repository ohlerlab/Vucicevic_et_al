# Vucicevic_et_al
Author: Dubravka Vucicevic

This repository contains scripts used in the manuscript titled: 'Dissecting transcriptional regulatory landscapes using bulk and targeted single-cell CRISPR activation screening'. 
Data generated for this study are accessible at GEO under accession numbers: GSE274254 for ATAC-seq, GSE274255 for bulk CRISPR screen, GSE274256 for TESLA-seq, GSE274257 for nuCapture-C and GSE274258 for chromatin RNA-seq.

## Buld screen alaysis
- bulk_screen_analysis/

    - bin_analysis/
        
        **num_cluster_lm.R**: File performs parallelized linear mixed-effects modeling on a set of bins and outputs the results to a text file.
    - mageck/
       
        Inputs and outputs from mageck analys
    
    - MLM/
        
        Multiple functions to perform MLM analysis used in (new_guide_selection_visualization.Rmd)
    
    - named_counts/
    
    - region_generation/
        
        **filtering_sig_bins_to_bed.R**: This script processes genomic data to identify significant bins based on p-values
        
        **merging.sh**: shell script that uses `bedtools` to merge overlapping or nearby intervals in a BED file. The intervals are merged if they are within 500 base pairs of each other, and the fourth column of the merged intervals will contain a comma-separated list of values from the original intervals. The input file is `sig_window_bins.bed`, and the output file is `merged_windowed_500.bed`
    
    - tab/

    - visualization/
        
        **heatmap_track_plots.R** : The script uses Gviz to create various DataTrack objects for different genomic features and combines them into plots. The plots are saved as PDF files, showing different genomic features and their significance. Region information is visualized as a heatmap.
        
        **plots_scripts.R** : 
             
             This script performs the following tasks:
                
                1. Load necessary libraries and set working directory.
                
                2. Load various datasets for VPR and KRAB constructs, including guides, windows, bins, and regions.
                
                3. Define plotting functions and themes for visualizations.
                
                4. Generate volcano plots for VPR and KRAB regions, windows, and guides.
                
                5. Create cross-tabulation plots for regions, windows, and guides.
                
                6. Identify and plot shared bins and guides between VPR and KRAB constructs.
                
                7. Save all generated plots to PDF files.
        
        **track_plots.R**: The script uses Gviz to create various DataTrack objects for different genomic features and combines them into plots. The plots are saved as PDF files, showing different genomic features and their significance. Region information is visualized as a dotplot.

    ### new_guide_selection_visualization.Rmd

        File performs the following tasks:

        1. Setup: Loads necessary libraries and sets global options.
        
        2. Data Loading: Loads and processes raw region data and results from various analyses (VPR, KRAB).
        
        3. Data Preparation: Merges and processes data to integrate all information for further analysis.
        
        4. Statistical Analysis: Runs linear models and other statistical analyses on the data.
        
        5. Visualization: Generates various plots to visualize the results, including region and bin results. 
        
        6. Guide Selection: Identifies and selects the top guides based on specific criteria.
        
        7. Output: Saves the final selection of guides and results to files for further use.

        Output files:
        
            MaGeCK Input Files:
        
                1. `mageck/input_mageck_vpr.txt`
        
                2. `mageck/input_mageck_krab.txt`
        
            Region and Bin Results:
        
                1. `regions_results.txt`
        
                2. `regions_bins_results.txt`
        
                3. `regions_bins_guides_results.txt`
        
            Files to be used with scripts in visualizations
        
                1. `visualization/VPR_guide_results_plot.csv`
        
                2. `visualization/KRAB_guide_results_plot.csv`
        
                3. `visualization/VPR_window_results_plot.csv`
        
                4. `visualization/KRAB_window_results_plot.csv`
        
                5. `visualization/VPR_region_results_plot.csv`
        
                6. `visualization/KRAB_region_results_plot.csv`
  
## TESLA-seq data analysis 
- Please reach out to Che-Wei Hsu, Ph.D. ([Duke University](che-wei.hsu@duke.edu)/[Caltech](cwhsu@caltech.edu)) for questions regarding the analysis

  ### In the TESLA-seq folder

        You can find three Jupyter notebooks documenting the codes of analysis, one R-script validating the numbers reported in the paper, and necessary input data and supplementary files to reproduce our results:

        1. Jupyter notebook named "Rhapsody-202210_batch_correct-Normalize_respectively-SCT_var_regress_updated_in_Feb_2025" documents most of the analysis, including preprocessing, normalization, quality checks, gRNA calling, DE and plotting.
        
        2. Jupyter notebook named "Rhapsody_Random_Forest_All_features_top_vs_bottom_50_wo_GC-20230216" documents the code for random forest model and logistic regression.
        
        3. Jupyter notebook named "Rhapsody_Sushi-Region" documents the code for plotting one of the figures in the paper.
        
        4. R-script named "Validate_numbers_in_manuscript" provides example of how to access our processed data and query for numbers or information reported in the paper.
        
    
  
