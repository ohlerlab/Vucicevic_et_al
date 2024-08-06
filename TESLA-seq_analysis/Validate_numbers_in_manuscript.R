rm(list=ls())
setwd('~/Downloads/')

library('tidyverse')
## Raw counts from Rhapsody
R1 <- read.table("./Rhapsody_20240103/UA9_UN_BD.csv", header=TRUE, sep=",")
R1 <- t(as.matrix(R1))
R2 <- read.table("./Rhapsody_20240103/UA9_UN_BD_2.csv", header=TRUE, sep=",")
R2 <- t(as.matrix(R2))

#ngenes
nrow(R1)
nrow(R2)
#ncells
ncol(R1)
ncol(R2)

rownames(R1) <- gsub("_gencode","",rownames(R1))
rownames(R2) <- gsub("_gencode","",rownames(R2))


t_id <- rev(rownames(R1))[1:78] ## 78 targeted genes
t_idx <- sort(match(t_id,rownames(R1)));
g_idx <- seq(1,nrow(R1),1)[-t_idx]
tmtx1 <- R1[t_idx,] ## gene matrix
gmtx1 <- R1[g_idx,] ## gRNA matrix


t_id <- rev(rownames(R2))[1:78]
t_idx <- sort(match(t_id,rownames(R2)));
g_idx <- seq(1,nrow(R2),1)[-t_idx]
tmtx2 <- R2[t_idx,]
gmtx2 <- R2[g_idx,]

##
nrow(tmtx1)
# [1] 78
nrow(tmtx2)
# [1] 78

summary(colSums(tmtx1))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 15.0   361.0   528.0   649.8   848.0  3223.0 
summary(colSums(gmtx1))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0     2.0    19.0   157.8   233.0  3063.0

summary(colSums(tmtx2))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3.0   116.0   169.0   199.6   250.0  1582.0
summary(colSums(gmtx2))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00    1.00   15.00   60.71   84.00 1738.00 

## gene annotation files
genes.bed <- read.table("./Rhapsody_20240103/Rhapsody_genes.bed")
colnames(genes.bed) <- c("chrom", "start", "stop","gene","score","strand","type")
genes.bed$gene[which(genes.bed$gene == "AC195454.1")] <- "LINC02265"
genes.bed$gene[which(genes.bed$gene == "UCHL1-AS1")] <- "UCHL1.AS1"
genes.bed$gene[which(genes.bed$gene == "LVCAT1")] <- "LINC02475"
genes.bed$gene[which(genes.bed$gene == "RP11-663P9.2")] <- "NDUFB4P12"
genes.bed$gene[which(genes.bed$gene == "ZBTB12P1")] <- "ZBTB12BP"
genes.bed <- rbind(genes.bed,data.frame(chrom="chr4",start=41748293,stop=41824119,gene="AC105389.3",score=".",strand="+",type="gene"))

br_anno <- read.table("./Rhapsody_20240103/new_rhapsody_regions_bins_guides_results.txt", sep="\t", header=TRUE)
br_anno <- br_anno[which((br_anno$class != "NEG-CONTROL") & (br_anno$class != "non-targctrl")),]

pr_anno <- read.table("./Rhapsody_20240103/rhapsodygRNAs_promoter500.txt", sep="\t", header=FALSE)
colnames(pr_anno) <- c("chr", "start", "end", "strand", "name", "chr_prom", "start_prom", "end_prom", "prom_name", "overlap")
pr_anno$prom_name[which(pr_anno$prom_name == "AC195454.1")] <- "LINC02265"
pr_anno$prom_name[which(pr_anno$prom_name == "UCHL1-AS1")] <- "UCHL1.AS1"
pr_anno$prom_name[which(pr_anno$prom_name == "LVCAT1")] <- "LINC02475"
pr_anno$prom_name[which(pr_anno$prom_name == "RP11-663P9.2")] <- "NDUFB4P12"
pr_anno$prom_name[which(pr_anno$prom_name == "ZBTB12P1")] <- "ZBTB12BP"

## DE results
TAP_DE <- read.csv("./Rhapsody_20240103/TAP_DE_aaf800_SCT_202212.csv")
head(TAP_DE)
TAP_DE <- TAP_DE %>% filter(p_val_adj < 0.05)

TAP_DE_PR <- TAP_DE[TAP_DE$ident %in% pr_anno$name,]
pr_anno$ident <- pr_anno$name

TAP_DE_PR <- left_join(TAP_DE_PR, pr_anno, by='ident')

TAP_DE_PR <- TAP_DE_PR %>% filter(gene!=prom_name)

head(TAP_DE_PR)

write.csv(TAP_DE_PR, 'Potential_promoters_act_as_enhancers.csv')

## 
rhits <- TAP_DE[grep("Region",TAP_DE$ident),]
nrow(rhits)
#[1] 92

summary(as.numeric(table(rhits$ident)))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   1.000   1.533   2.000   5.000 

length(table(rhits$ident)[table(rhits$ident)>1])
#[1] 22

summary(as.numeric(table(rhits$gene)))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   2.000   2.788   3.000  15.000 

ghits <- rhits
chrom1 <- br_anno$Chr[match(ghits$ident, br_anno$Region_name)]
start1 <- br_anno$Region_start[match(ghits$ident, br_anno$Region_name)]
end1 <- br_anno$Region_end[match(ghits$ident, br_anno$Region_name)]
chrom2 <- genes.bed$chrom[match(ghits$gene, genes.bed$gene)]
start2 <- genes.bed$start[match(ghits$gene, genes.bed$gene)]
end2 <- genes.bed$stop[match(ghits$gene, genes.bed$gene)]
ghits$ident <- gsub("^Region_","R",ghits$ident)
name <- paste0(ghits$ident,"_",ghits$gene)
ghits$ident <- gsub("^R","Region_",ghits$ident)
score <- ghits$avg_logFC
strand1 <- br_anno$strand[match(ghits$ident, br_anno$Region_name)]
strand2 <- genes.bed$strand[match(ghits$gene, genes.bed$gene)]
samplenumber <- as.numeric(as.factor(ghits$gene))
gene <- ghits$gene
gRNA <- ghits$ident
log_pvaladj <- -log(ghits$p_val_adj)
type <- ghits$avg_logFC
type[which(type > 0)] <- "up_regulated"
type[which(type < 0)] <- "down_regulated"
br_anno$Region_name <- gsub("^Region_","R",br_anno$Region_name)
gRNA.bed <- data.frame(chrom=br_anno$Chr, start=br_anno$Region_start, stop=br_anno$Region_end, gene=br_anno$Region_name, score=".", strand=br_anno$strand, type="gene")
gRNA.bed <- gRNA.bed[!is.na(gRNA.bed$gene),]
gRNA.bed <- gRNA.bed %>% distinct()
#ALL
ghits.bedpe <- data.frame(chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2, samplenumber, gene, gRNA, log_pvaladj, type)

## Jumped genes
ghits.bedpe$distance.TSS <- apply(ghits.bedpe, 1, function(x){if(x[9]=="+"&x[10]=="+"){abs(as.numeric(x[2])-as.numeric(x[5]))}else if(x[9]=="+"&x[10]=="-"){abs(as.numeric(x[2])-as.numeric(x[6]))}else if(x[9]=="-"&x[10]=="+"){abs(as.numeric(x[3])-as.numeric(x[5]))}else{abs(as.numeric(x[3])-as.numeric(x[6]))}})
genes.bed.detected <- genes.bed[genes.bed$gene %in% unique(TAP_DE$gene),]

jumped_gene_counts <- c()
count <- 0
for (i in 1:nrow(ghits.bedpe)){
  for (j in 1:nrow(genes.bed.detected)){
    if(genes.bed.detected$start[j]>=min(c(ghits.bedpe$start1[i],ghits.bedpe$start2[i],ghits.bedpe$end1[i],ghits.bedpe$end2[i])) & genes.bed.detected$stop[j]<=max(c(ghits.bedpe$start1[i],ghits.bedpe$start2[i],ghits.bedpe$end1[i],ghits.bedpe$end2[i])) & ghits.bedpe$gene[i]!=genes.bed.detected$gene[j]){
      count <- count + 1
    }
    if(j == nrow(genes.bed.detected)){
      jumped_gene_counts <- c(jumped_gene_counts, count)
      count <- 0
    }
  }
}

table(jumped_gene_counts)
#jumped_gene_counts
#0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 23 24 27 33 34 36 
#27 12  6  3  3  4  3  1  1  2  1  1  2  2  2  4  3  1  1  1  2  1  2  1  1  2  2  1 

## Skipped genes
length(jumped_gene_counts)-27
#[1] 65

options(repr.plot.width=8, repr.plot.height=8)
boxplot(jumped_gene_counts, main="number of jumped genes", cex.lab=1.5)
summary(jumped_gene_counts)

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   0.000   3.000   7.663  14.000  36.000
