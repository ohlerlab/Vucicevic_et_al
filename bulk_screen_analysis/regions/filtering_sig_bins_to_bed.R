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