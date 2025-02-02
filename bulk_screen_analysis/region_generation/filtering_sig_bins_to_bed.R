# This script processes genomic data to identify significant bins based on p-values.
# 
# Steps:
# 1. Load necessary libraries.
# 2. Read input data containing p-values and bin coordinates.
# 3. Adjust p-values using the False Discovery Rate (FDR) method.
# 4. Merge the p-value data with bin coordinates.
# 5. Filter for significant bins with FDR < 0.05.
# 6. Select relevant columns (chr, start, end, bin) from the filtered data.
# 7. Write the filtered significant bins to a BED file.
#
# Input files:
# - '../sliding_window/num_VPR_pvalues.txt': Contains bin, estimate, and p-value data.
# - '../bin_coords.txt': Contains chromosome, start, end, and bin coordinates.
#
# Output file:
# - 'sig_window_bins.bed': Contains significant bins with columns chr, start, end, and bin.
library(tidyverse)
# loading the information needed, coords are from the 'funcitonal' ones not the guide design ones. 
input <- read_delim('../sliding_window/num_VPR_pvalues.txt',delim='\t',col_names=c('bin','estimate','pvalue'))
input$fdr <- p.adjust(input$pvalue,method='fdr')
coords <- read_delim('../bin_coords.txt',delim='\t',col_names=c('chr','start','end','bin'))

# merging
full_results <- left_join(coords,input,by='bin')

# filtering
significant  <- dplyr::filter(full_results,fdr<0.05)

# output
output <- significant %>%
	select(chr,start,end,bin)

#writing out
write_delim(output,'sig_window_bins.bed',delim='\t',col_name=FALSE)