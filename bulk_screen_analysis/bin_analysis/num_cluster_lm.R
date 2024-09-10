library(tidyverse)
library(lme4)
library(plyr)
#guides<-read_delim('~/Documents/Phox2b/tab/phox2bLocus3.tsv',delim='\t',col_names = c('NAME','GUIDE'))
guides<-read_delim('../tab/phox2bLocus3.tsv',delim='\t',col_names = c('NAME','GUIDE'))
#codes <- read_delim('~/Documents/Phox2b/codes.txt',delim=' ',col_names = c('GUIDE','CODE'))
codes <- read_delim('../codes.txt',delim=' ',col_names = c('GUIDE','CODE'))
codes$BIN <- sapply(strsplit(codes$CODE,'_'),'[',1)
codes$REPLICATE<- sapply(strsplit(codes$CODE,'_'),'[',3)
#codes <- dplyr::select(codes,-CODE)
head(codes)
all_binx <- unique(codes$BIN)
w_bins <- lapply(1:(length(all_binx)-2), function(i){
  all_binx[(i):(i+2)]
}) 
names(w_bins)<-all_binx[2:(length(all_binx)-1)]
#KRAB_counts<-read_delim('~/Documents/Phox2b/named_counts/KRAB_table.txt',delim=',')
KRAB_counts<-read_delim('../named_counts/KRAB_table.txt',delim=',')
KRAB_full <- inner_join(guides,KRAB_counts)
#VPR_counts<-read_delim('~/Documents/Phox2b/named_counts/VPR_table.txt',delim=',')
VPR_counts<-read_delim('../named_counts/VPR_table.txt',delim=',')
VPR_full <- inner_join(guides,VPR_counts)
KRAB <- KRAB_full %>% mutate(K_5D=(`KRAB_D5-1_S1`+`KRAB_D5-2_S2`)/2,K_20D=(`KRAB_D20_S3`),K_29D=(`KRAB_D29-1_S4`+`KRAB_D29-2_S5`)/2,K_33D=(`KRAB_D33-1_S6`+`KRAB_D33-2_S7`)/2) %>% select(NAME,GUIDE,K_5D,K_20D,K_29D,K_33D)
VPR <- VPR_full %>% mutate(V_5D=(`VPR_D5-1_S8`+`VPR_D5-2_S9`)/2,V_20D=(`VPR_D20_S10`),V_29D=(`VPR_D29-1_S11`+`VPR_D29-2_S12`)/2,V_33D=(`VPR_D33-1_S13`+`VPR_D33-2_S14`)/2) %>% select(NAME,GUIDE,V_5D,V_20D,V_29D,V_33D)
lorNorm <- function(matrix){
  n <- apply(matrix,2,function(x) (x/sum(x))*1e6)
return(n)}
##LINEAR MODEL
VPR_norm <-cbind(VPR[,1:2],log(lorNorm(VPR[,3:6])+1))
coded_VPR <- inner_join(VPR_norm,codes)
#coded_VPR$BIN <- as.factor(coded_VPR$BIN)
coded_VPR$REPLICATE <- as.factor(coded_VPR$REPLICATE)
#m_coded_VPR <- coded_VPR %>% group_by(BIN) %>% filter(n()>1)
coded_VPR_df <- gather(coded_VPR,key='Group',value='LogCounts',-c(colnames(coded_VPR)[c(1:2,7:9)]))
#coded_VPR_df$Time <- factor(sapply(coded_VPR_df$Group,str_extract,pattern='[0-9]+'),ordered = T,levels = c(5,20,29,33))
coded_VPR_df$Time <- as.numeric(sapply(coded_VPR_df$Group,str_extract,pattern='[0-9]+'))
coded_VPR_df$Time <- ifelse(is.na(coded_VPR_df$Time),0,coded_VPR_df$Time)
coded_VPR_df$CODE <- as.factor(coded_VPR_df$CODE)
coded_VPR_df = coded_VPR_df %>% mutate_if(is.factor,
                      fct_explicit_na,
                      na_level = "missing")
#PVALUE FINDING
library(parallel)
 no_core <- 8 
 cl <- makeCluster(no_core, type="FORK")
 pvals_guide <- parLapply(cl,names(w_bins),function(bin){
 #pvals_guide <- lapply(names(w_bins)[8408:8508],function(bin){
   df <- filter(coded_VPR_df,BIN %in% w_bins[[bin]])
   #df <- filter(coded_VPR_df,BIN == bin)
   LM.mixed = lmer(LogCounts ~ Time  + (1|CODE),data=df, REML=F)
   LM.null = lmer(LogCounts ~  (1|CODE),data=df, REML=F)
   A=anova(LM.null,LM.mixed)
   #print(bin)
   results <- data.frame('BIN'=bin, 'coeff'=summary(LM.mixed)$coefficients[2,1], 'pvalue'=A$`Pr(>Chisq)`[2])
   #return(A$`Pr(>Chisq)`[2])
   return(results)
 })
 stopCluster(cl)
 results <-ldply (pvals_guide, data.frame) 
 #results <- data.frame('bin'=names(w_bins),'pvalue'=pvals_guide)
 write.table(results,"num_VPR_pvalues.txt",sep = "\t",row.names = FALSE,col.names = TRUE,quote = FALSE)
#EARLY VS LATE
# results_time <- coded_VPR_df %>%
#   mutate(timeearly = Time%in%c('5'),timelate = Time %in% c('20','29','33'))
# no_core <- 8 
# cl <- makeCluster(no_core, type="FORK")
# speed_bin <- parLapply(cl,names(w_bins),function(bin){
#   df <- filter(results_time,BIN %in% w_bins[[bin]])
#   LM.mixed_both = lmer(LogCounts ~ timeearly + timelate + (1|CODE),data = df,REML=F)
#   LM.mixed_early = lmer(LogCounts ~ timeearly  + (1|CODE),data = df,REML=F)
#   LM.mixed_late = lmer(LogCounts ~ timelate + (1|CODE),data = df,REML=F)
#   early_genes = (anova(LM.mixed_both,LM.mixed_late)$`Pr(>Chisq)`[2])<0.05
#   late_genes = (anova(LM.mixed_both,LM.mixed_early)$`Pr(>Chisq)`[2])<0.05
#   results<-ifelse(early_genes&late_genes,'early',ifelse(early_genes,'early',ifelse(late_genes,'late','stable')))
#   return(results)
# })
# stopCluster(cl)
# results_speed <- ldply(speed_bin,data.frame)
# write.table(results_speed,'speed_VPR_pvalues.txt',sep = "\t",row.names = F,col.names = T,quote = FALSE)



#LINEAR MODEL
KRAB_norm <-cbind(KRAB[,1:2],log(lorNorm(KRAB[,3:6])+1))
coded_KRAB <- inner_join(KRAB_norm,codes)
#coded_KRAB$BIN <- as.factor(coded_KRAB$BIN)
coded_KRAB$REPLICATE <- as.factor(coded_KRAB$REPLICATE)
#m_coded_KRAB <- coded_KRAB %>% group_by(BIN) %>% filter(n()>1)
coded_KRAB_df <- gather(coded_KRAB,key='Group',value='LogCounts',-c(colnames(coded_KRAB)[c(1:2,7:9)]))
#coded_KRAB_df$Time <- factor(sapply(coded_KRAB_df$Group,str_extract,pattern='[0-9]+'),ordered = T,levels = c(5,20,29,33))
coded_KRAB_df$Time <- as.numeric(sapply(coded_KRAB_df$Group,str_extract,pattern='[0-9]+'))
coded_KRAB_df$Time <- ifelse(is.na(coded_KRAB_df$Time),0,coded_KRAB_df$Time)
coded_KRAB_df$CODE <- as.factor(coded_KRAB_df$CODE)
coded_KRAB_df = coded_KRAB_df %>% mutate_if(is.factor,
                      fct_explicit_na,
                      na_level = "missing")
# #str(coded_KRAB_df)
# library(parallel)
# no_core <- 8 
# cl <- makeCluster(no_core, type="FORK")
# pvals_guide <- parLapply(cl,names(w_bins),function(bin){
# #pvals_guide <- lapply(names(w_bins)[8408:8508],function(bin){
#   df <- filter(coded_KRAB_df,BIN %in% w_bins[[bin]])
#   #df <- filter(coded_VPR_df,BIN == bin)
#   LM.mixed = lmer(LogCounts ~ Time  + (1|CODE),data=df, REML=F)
#   LM.null = lmer(LogCounts ~  (1|CODE),data=df, REML=F)
#   A=anova(LM.null,LM.mixed)
#   #print(bin)
#   results <- data.frame('BIN'=bin, 'coeff'=summary(LM.mixed)$coefficients[2,1], 'pvalue'=A$`Pr(>Chisq)`[2])
#   #return(A$`Pr(>Chisq)`[2])
#   return(results)
# })
# stopCluster(cl)
# results <-ldply (pvals_guide, data.frame) 
# write.table(results,"num_KRAB_pvalues.txt",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
#EARLY VS LATE

results_time <- coded_KRAB_df %>%
  mutate(timeearly = Time%in%c('5'),timelate = Time %in% c('20','29','33'))
no_core <- 8 
cl <- makeCluster(no_core, type="FORK")
speed_bin <- parLapply(cl,names(w_bins),function(bin){
  df <- filter(results_time,BIN %in% w_bins[[bin]])
  LM.mixed_both = lmer(LogCounts ~ timeearly + timelate + (1|CODE),data = df,REML=F)
  LM.mixed_early = lmer(LogCounts ~ timeearly  + (1|CODE),data = df,REML=F)
  LM.mixed_late = lmer(LogCounts ~ timelate + (1|CODE),data = df,REML=F)
  early_genes = (anova(LM.mixed_both,LM.mixed_late)$`Pr(>Chisq)`[2])<0.05
  late_genes = (anova(LM.mixed_both,LM.mixed_early)$`Pr(>Chisq)`[2])<0.05
  results<-ifelse(early_genes&late_genes,'early',ifelse(early_genes,'early',ifelse(late_genes,'late','stable')))
  return(results)
})
stopCluster(cl)
results_speed <- ldply(speed_bin,data.frame)
write.table(results_speed,'speed_KRAB_pvalues.txt',sep = "\t",row.names = F,col.names = T,quote = FALSE)
