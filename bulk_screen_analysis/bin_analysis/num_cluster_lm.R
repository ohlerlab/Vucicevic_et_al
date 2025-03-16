# This script performs linear mixed model analysis on VPR and KRAB datasets.
# It reads in guide and code data, normalizes the counts, and fits linear mixed models
# to identify significant bins based on time points.

# Load necessary libraries
library(tidyverse)
library(lme4)
library(plyr)
library(parallel)

# Read in guide and code data
guides <- read_delim('../tab/phox2bLocus3.tsv', delim='\t', col_names=c('NAME', 'GUIDE'))
codes <- read_delim('../codes.txt', delim=' ', col_names=c('GUIDE', 'CODE'))

# Process code data to extract BIN and REPLICATE information
codes$BIN <- sapply(strsplit(codes$CODE, '_'), '[', 1)
codes$REPLICATE <- sapply(strsplit(codes$CODE, '_'), '[', 3)

# Generate a list of unique bins and create sliding windows of bins
all_binx <- unique(codes$BIN)
w_bins <- lapply(1:(length(all_binx) - 2), function(i) {
  all_binx[(i):(i + 2)]
})
names(w_bins) <- all_binx[2:(length(all_binx) - 1)]

# Read in KRAB and VPR count data
KRAB_counts <- read_delim('../named_counts/KRAB_table.txt', delim=',')
VPR_counts <- read_delim('../named_counts/VPR_table.txt', delim=',')

# Merge guide data with count data
KRAB_full <- inner_join(guides, KRAB_counts)
VPR_full <- inner_join(guides, VPR_counts)

# Normalize and log-transform the count data
lorNorm <- function(matrix) {
  n <- apply(matrix, 2, function(x) (x / sum(x)) * 1e6)
  return(n)
}

# Prepare VPR data for linear model analysis
VPR_norm <- cbind(VPR[, 1:2], log(lorNorm(VPR[, 3:6]) + 1))
coded_VPR <- inner_join(VPR_norm, codes)
coded_VPR$REPLICATE <- as.factor(coded_VPR$REPLICATE)
coded_VPR_df <- gather(coded_VPR, key='Group', value='LogCounts', -c(colnames(coded_VPR)[c(1:2, 7:9)]))
coded_VPR_df$Time <- as.numeric(sapply(coded_VPR_df$Group, str_extract, pattern='[0-9]+'))
coded_VPR_df$Time <- ifelse(is.na(coded_VPR_df$Time), 0, coded_VPR_df$Time)
coded_VPR_df$CODE <- as.factor(coded_VPR_df$CODE)
coded_VPR_df <- coded_VPR_df %>% mutate_if(is.factor, fct_explicit_na, na_level="missing")

# Perform parallel linear mixed model analysis on VPR data
no_core <- 8
cl <- makeCluster(no_core, type="FORK")
pvals_guide <- parLapply(cl, names(w_bins), function(bin) {
  df <- filter(coded_VPR_df, BIN %in% w_bins[[bin]])
  LM.mixed <- lmer(LogCounts ~ Time + (1 | CODE), data=df, REML=F)
  LM.null <- lmer(LogCounts ~ (1 | CODE), data=df, REML=F)
  A <- anova(LM.null, LM.mixed)
  results <- data.frame('BIN'=bin, 'coeff'=summary(LM.mixed)$coefficients[2, 1], 'pvalue'=A$`Pr(>Chisq)`[2])
  return(results)
})
stopCluster(cl)
results <- ldply(pvals_guide, data.frame)
write.table(results, "num_VPR_pvalues.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# Prepare KRAB data for linear model analysis
KRAB_norm <- cbind(KRAB[, 1:2], log(lorNorm(KRAB[, 3:6]) + 1))
coded_KRAB <- inner_join(KRAB_norm, codes)
coded_KRAB$REPLICATE <- as.factor(coded_KRAB$REPLICATE)
coded_KRAB_df <- gather(coded_KRAB, key='Group', value='LogCounts', -c(colnames(coded_KRAB)[c(1:2, 7:9)]))
coded_KRAB_df$Time <- as.numeric(sapply(coded_KRAB_df$Group, str_extract, pattern='[0-9]+'))
coded_KRAB_df$Time <- ifelse(is.na(coded_KRAB_df$Time), 0, coded_KRAB_df$Time)
coded_KRAB_df$CODE <- as.factor(coded_KRAB_df$CODE)
coded_KRAB_df <- coded_KRAB_df %>% mutate_if(is.factor, fct_explicit_na, na_level="missing")

# Perform parallel linear mixed model analysis on KRAB data
cl <- makeCluster(no_core, type="FORK")
pvals_guide <- parLapply(cl, names(w_bins), function(bin) {
  df <- filter(coded_KRAB_df, BIN %in% w_bins[[bin]])
  LM.mixed <- lmer(LogCounts ~ Time + (1 | CODE), data=df, REML=F)
  LM.null <- lmer(LogCounts ~ (1 | CODE), data=df, REML=F)
  A <- anova(LM.null, LM.mixed)
  results <- data.frame('BIN'=bin, 'coeff'=summary(LM.mixed)$coefficients[2, 1], 'pvalue'=A$`Pr(>Chisq)`[2])
  return(results)
})
stopCluster(cl)
results <- ldply(pvals_guide, data.frame)
write.table(results, "num_KRAB_pvalues.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
