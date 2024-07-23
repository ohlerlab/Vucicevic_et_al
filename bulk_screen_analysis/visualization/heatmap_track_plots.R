# ============================================================================
# Track plots
# ============================================================================

## Report Plotting Test
library(tidyverse)
library(rtracklayer)
library(AnnotationHub)
library(Gviz)
library(Rsamtools)
library(biomaRt)
# Parameters
setwd("~/Documents/Phox2b/PAPER/visualization/")
# ============================================================================
# Loading results
# ============================================================================
# ==== Guides =====
VPR_guide_results_plot <- read_csv("VPR_guide_results_plot.csv")
KRAB_guide_results_plot <- read_csv("KRAB_guide_results_plot.csv")
# ==== Bins =====
bin_coords <- read_delim('../bin_coords.txt',delim='\t',col_names=c('chr','start','end','name'))
# ==== Bins =====
VPR_bin_results_plot <- read_csv("CT_Results_input_VPR.txt")
KRAB_bin_results_plot <- read_csv("CT_Results_input_KRAB.txt")
VPR_window_results_plot <- read_csv("VPR_window_results_plot.csv")
KRAB_window_results_plot <- read_csv("KRAB_window_results_plot.csv")
# ==== Regions =====
VPR_results_plot <- read_csv("VPR_regions_results_plot.csv")
KRAB_results_plot <- read_csv("KRAB_regions_results_plot.csv")

# ==== colors =====
cold_colors <- c(
  "#003399",
  "#00528A",
  "#007877",
  "#009966"
)
my_colors_tracks <- c("#013766", "white", "#ac0e28")
names(my_colors_tracks)<-c('down','nc','up')

# ==== Annotation tracks =====

ATAC <- import.bw("~/Documents/Phox2b/bigwigs/SH-SY5Y_ATAC_inhouse.steric38.hg38.bw")
PHOX2B <- import.bw("~/Documents/Phox2b/bigwigs/GSM2664369_CLBGA.PHOX2B.rep1.wig.bw")
GATA3 <- import.bw("~/Documents/Phox2b/bigwigs/GSM2664370_CLBGA.GATA3.rep1.wig.bw")
HAND2 <- import.bw("~/Documents/Phox2b/bigwigs/GSM2664371_CLBGA.HAND2.rep1.wig.bw")
H3K27AC <- import.bw("~/Documents/Phox2b/bigwigs/SH-SY5Y_H3K27ac_boeva2017.hg38.bw")

ATACtrack <- DataTrack(
  range = ATAC,
  name = "        ATAC",
  type = "histogram",
  col.histogram = cold_colors[1],
  fill.histogram = cold_colors[1],
  background.title = cold_colors[1]
)
PHOX2Btrack <- DataTrack(
  range = PHOX2B,
  name = "           PHOX2B",
  type = "histogram",
  col.histogram = cold_colors[3],
  fill.histogram = cold_colors[3],
  background.title = cold_colors[3]
)
GATA3track <- DataTrack(
  range = GATA3,
  name = "        GATA3",
  type = "histogram",
  col.histogram = cold_colors[3],
  fill.histogram = cold_colors[3],
  background.title = cold_colors[3]
)
HAND2track <- DataTrack(
  range = HAND2,
  name = "          HAND2",
  type = "histogram",
  col.histogram = cold_colors[4],
  fill.histogram = cold_colors[4],
  background.title = cold_colors[4]
)
H3K27ACtrack <- DataTrack(
  range = H3K27AC,
  name = "             H3K27AC",
  type = "histogram",
  col.histogram = cold_colors[2],
  fill.histogram = cold_colors[2],
  background.title = cold_colors[2]
)

# ==== Fixed parameters =====
info <- Seqinfo(seqname = "chr4", genome = "hg38")
chr <- "chr4"
gen <- "hg38"
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = gen, chromosome = chr)
biomart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
geneTrack <- BiomartGeneRegionTrack(chromosome = chr, genome = gen, biomart = biomart, name = "GeneTrack", background.title = "black", stacking='squish')
sup_df <- data.frame("start" = c(41681648, 41320748), "width" = c(213253, 326103))
superenhancer <- AnnotationTrack(
  start = sup_df$start,
  width = sup_df$width,
  chromosome = chr,
  strand = rep("*", 2),
  group = c("Enh1", "Enh2"),
  genome = gen, name = "Boeva",
  background.title = "#7C0A02",
  color = "#7C0A02"
)
feature(superenhancer) <- c("Boeva")
hl_phox <- HighlightTrack(geneTrack, start = 41746346, width = 5399, chromosome = chr)

# ============================================================================
# Regions
# ============================================================================
regions <- readr::read_delim("~/Documents/Phox2b/merged_windows/merged_windowed_500.bed", delim = "\t", col_names = c("Chr", "start", "end", "bins"))
regions <- dplyr::mutate(regions, name = paste("Region", seq(1, nrow(regions), 1), sep = "_")) %>%
  dplyr::mutate(BINS = str_split(bins, pattern = ",")) %>%
  tidyr::unnest(cols = c(BINS))
regions_gviz <- regions %>%
  dplyr::select(Chr, start, end, name) %>%
  distinct() %>%
  full_join(VPR_results_plot, by = c("name" = "Gene")) %>%
  replace_na(list(pvalue = 1, fdr = 1, effect = "fixed", estimate = 0))
regions_gviz_k <- regions %>%
  dplyr::select(Chr, start, end, name) %>%
  distinct() %>%
  full_join(KRAB_results_plot, by = c("name" = "Gene")) %>%
  replace_na(list(pvalue = 1, fdr = 1, effect = "fixed", estimate = 0))

scores <- -log10(regions_gviz$fdr)[-944]
significant <- as.numeric(regions_gviz$fdr < 0.05)[-944]
coeffs <- regions_gviz$estimate[-944]

df_regions <- regions_gviz %>%
  mutate("score" = -log10(fdr)) %>%
  mutate("color" = ifelse(fdr < 0.05, ifelse(estimate < 0, "#013766", "#ac0e28"), "#D2D2D2FF")) %>%
  mutate("sig" = as.numeric(fdr < 0.05)) %>%
  drop_na("name")

scores_k <- -log10(regions_gviz_k$fdr)[-944]
significant_k <- as.numeric(regions_gviz_k$fdr < 0.05)[-944]
coeffs_k <- regions_gviz_k$estimate[-944]

df_regions_k <- regions_gviz_k %>%
  mutate("score" = -log10(fdr)) %>%
  mutate("color" = ifelse(fdr < 0.05, ifelse(estimate < 0, "#013766", "#ac0e28"), "#D2D2D2FF")) %>%
  mutate("sig" = as.numeric(fdr < 0.05)) %>%
  drop_na("name")

annot_regions <- makeGRangesFromDataFrame(regions_gviz[-944, ], seqinfo = info)

dataTrack_regions <- DataTrack(range = annot_regions,
                          data = df_regions$estimate,
                          name = "        Region",
                          type = "heatmap",
                          ylim = c(-0.04,0.04),
                          gradient = c('#1683D0','white','#d60000'),
                          background.title='#e76f3d')
dataTrack_regions_sig <- DataTrack(range = annot_regions,
                          data = df_regions$sig,
                          name = "        Region sig",
                           type = "heatmap",
                           gradient = c('white','black'),background.title='black',
                           background.title ='#e76f3d',
                           ylim = c(0,1))

dataTrack_regions_k <- DataTrack(range = annot_regions,
                          data = df_regions_k$estimate,
                          name = "        Region",
                          type = "heatmap",
                          #ylim = c(-0.04,0.04),
                          gradient = c('#1683D0','white','#d60000'),
                          background.title='#1683D0')
dataTrack_regions_sig_k <- DataTrack(range = annot_regions,
                          data = df_regions_k$sig,
                          name = "        Region sig",
                           type = "heatmap",
                           gradient = c('white','black'),background.title='black',
                           background.title ='#1683D0',
                           ylim = c(0,1))
# ============================================================================
# windows
# ============================================================================

windows <- bin_coords
windows_gviz <- windows %>%
  dplyr::select(chr, start, end, name) %>%
  distinct() %>%
  full_join(VPR_window_results_plot, by = c("name" = "bin")) %>%
  replace_na(list(pvalue = 1, fdr = 1, effect = "fixed", estimate = 0))
windows_gviz_k <- windows %>%
  dplyr::select(chr, start, end, name) %>%
  distinct() %>%
  full_join(KRAB_window_results_plot, by = c("name" = "bin")) %>%
  replace_na(list(pvalue = 1, fdr = 1, effect = "fixed", estimate = 0))

scores <- -log10(windows_gviz$fdr)[-944]
significant <- as.numeric(windows_gviz$fdr < 0.05)[-944]
coeffs <- windows_gviz$estimate[-944]

df_windows <- windows_gviz %>%
  mutate("score" = -log10(fdr)) %>%
  mutate("color" = ifelse(fdr < 0.05, ifelse(estimate < 0, "#013766", "#ac0e28"), "#D2D2D2FF")) %>%
  mutate("sig" = as.numeric(fdr < 0.05)) %>%
  drop_na("name")

scores_k <- -log10(windows_gviz_k$fdr)[-944]
significant_k <- as.numeric(windows_gviz_k$fdr < 0.05)[-944]
coeffs_k <- windows_gviz_k$estimate[-944]

df_windows_k <- windows_gviz_k %>%
  mutate("score" = -log10(fdr)) %>%
  mutate("color" = ifelse(fdr < 0.05, ifelse(estimate < 0, "#013766", "#ac0e28"), "#D2D2D2FF")) %>%
  mutate("sig" = as.numeric(fdr < 0.05)) %>%
  drop_na("name")

annot_windows <- makeGRangesFromDataFrame(windows_gviz, seqinfo = info)

dataTrack_windows <- DataTrack(range = annot_windows,
                          data = df_windows$estimate,
                          name = "        Windows",
                          type = "heatmap",
                          ylim = c(-0.04,0.04),
                          gradient = c('#1683D0','white','#d60000'),
                          background.title='#e76f3d')
dataTrack_windows_sig <- DataTrack(range = annot_windows,
                          data = df_windows$sig,
                          name = "        Windows sig",
                           type = "heatmap",
                           gradient = c('white','black'),background.title='black',
                           background.title ='#e76f3d',
                           ylim = c(0,1))

dataTrack_windows_k <- DataTrack(range = annot_windows,
                          data = df_windows_k$estimate,
                          name = "        Windows",
                          type = "heatmap",
                          #ylim = c(-0.04,0.04),
                          gradient = c('#1683D0','white','#d60000'),
                          background.title='#1683D0')
dataTrack_windows_sig_k <- DataTrack(range = annot_windows,
                          data = df_windows_k$sig,
                          name = "        Windows sig",
                           type = "heatmap",
                           gradient = c('white','black'),background.title='black',
                           background.title ='#1683D0',
                           ylim = c(0,1))
# ============================================================================
# single_bins
# ============================================================================
bins <- bin_coords
bins_gviz <- VPR_bin_results_plot %>%
  replace_na(list(pvalue = 1, adj_pvalue = 1, effect = "fixed", estimate = 0))
bins_gviz_k <- KRAB_bin_results_plot %>%
  replace_na(list(pvalue = 1, adj_pvalue = 1, effect = "fixed", estimate = 0))
scores <- -log10(bins_gviz$adj_pvalue)
scores_k <- -log10(bins_gviz_k$adj_pvalue)
significant <- as.numeric(bins_gviz$adj_pvalue < 0.05)
significant_k <- as.numeric(bins_gviz_k$adj_pvalue < 0.05)
coeffs <- bins_gviz$estimate
coeffs_k <- bins_gviz_k$estimate

df_bins <- bins_gviz %>%
  mutate("score" = -log10(adj_pvalue)) %>%
  mutate("color" = ifelse(adj_pvalue < 0.05, ifelse(estimate < 0, "#013766", "#ac0e28"), "#D2D2D2FF")) %>%
  mutate("sig" = as.numeric(adj_pvalue < 0.05)) %>%
  drop_na("Gene")
df_bins_k <- bins_gviz_k %>%
  mutate("score" = -log10(adj_pvalue)) %>%
  mutate("color" = ifelse(adj_pvalue < 0.05, ifelse(estimate < 0, "#013766", "#ac0e28"), "#D2D2D2FF")) %>%
  mutate("sig" = as.numeric(adj_pvalue < 0.05)) %>%
  drop_na("Gene")
bin_coords_bins <- bin_coords[bin_coords$name%in%VPR_bin_results_plot$Gene,]
annot_bins <- makeGRangesFromDataFrame(bin_coords_bins, seqinfo = info)

dataTrack_bins <- DataTrack(range = annot_bins,
                          data = df_bins$estimate,
                          name = "        bins",
                          type = "heatmap",
                          ylim = c(-0.1,0.1),
                          gradient = c('#1683D0','white','#d60000'),
                          background.title='#e76f3d')
dataTrack_bins_sig <- DataTrack(range = annot_bins,
                          data = df_bins$sig,
                          name = "        bins sig",
                           type = "heatmap",
                           gradient = c('white','black'),background.title='black',
                           background.title ='#e76f3d',
                           ylim = c(0,1))

dataTrack_bins_k <- DataTrack(range = annot_bins,
                          data = df_bins_k$estimate,
                          name = "        bins",
                          type = "heatmap",
                          #ylim = c(-0.1,0.1),
                          gradient = c('#1683D0','white','#d60000'),
                          background.title='#1683D0')
dataTrack_bins_sig_k <- DataTrack(range = annot_bins,
                          data = df_bins_k$sig,
                          name = "        bins sig",
                           type = "heatmap",
                           gradient = c('white','black'),background.title='black',
                           background.title ='#1683D0',
                           ylim = c(0,1))
# ============================================================================
# Actual plots
# ============================================================================
# plotTracks(list(gtrack,dataTrack,geneTrack),transcriptAnnotation = "symbol",from = 40749815, to = 42748835, collapseTranscripts="longest", col.line = NULL)
pdf("VPR_KRAB_heatmaps_separeate.pdf", width = 10, height = 6)
# Whole TAD
#plotTracks(list(dataTrack_regions, dataTrack_regions_sig, dataTrack_regions_k, dataTrack_regions_sig_k, dataTrack_windows, dataTrack_windows_sig, dataTrack_windows_k, dataTrack_windows_sig_k, dataTrack_bins, dataTrack_bins_sig, dataTrack_bins_k, dataTrack_bins_sig_k, ATACtrack, H3K27ACtrack, PHOX2Btrack, HAND2track, gtrack, superenhancer,geneTrack), from = 40749815, to = 42748835, sizes = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2), groupAnnotation = "group", Boeva = "#7C0A02", legend = TRUE, cex.title = 0.8, rot.title = 1, title.width = 2)
plotTracks(list(dataTrack_regions, dataTrack_regions_sig,dataTrack_windows, dataTrack_windows_sig,dataTrack_bins, dataTrack_bins_sig, ATACtrack, H3K27ACtrack, PHOX2Btrack, HAND2track, gtrack, superenhancer, hl_phox), from = 40749815, to = 42748835, sizes = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,2), groupAnnotation = "group", Boeva = "#7C0A02", transcriptAnnotation = "symbol",legend = TRUE, cex.title = 0.5, rot.title = 1,collapseTranscripts=TRUE, title.width = 2)
plotTracks(list(dataTrack_regions_k, dataTrack_regions_sig_k,dataTrack_windows_k, dataTrack_windows_sig_k,dataTrack_bins_k, dataTrack_bins_sig_k, ATACtrack, H3K27ACtrack, PHOX2Btrack, HAND2track, gtrack, superenhancer, hl_phox), from = 40749815, to = 42748835, sizes = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), groupAnnotation = "group", Boeva = "#7C0A02", legend = TRUE, cex.title = 0.8, rot.title = 1,transcriptAnnotation = "symbol", collapseTranscripts=TRUE, title.width = 2)
dev.off()
# Promoter
pdf("VPR_KRAB_heatmaps_promoter_separate.pdf", width = 10, height = 5)
# plotTracks(list(dataTrack_regions, dataTrack_regions_sig, dataTrack_regions_k, dataTrack_regions_sig_k, dataTrack_windows, dataTrack_windows_sig, dataTrack_windows_k, dataTrack_windows_sig_k, dataTrack_bins, dataTrack_bins_sig, dataTrack_bins_k, dataTrack_bins_sig_k, ATACtrack, H3K27ACtrack, PHOX2Btrack, HAND2track, gtrack, superenhancer,geneTrack), from = 41734346, to = 41763745, sizes = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2), groupAnnotation = "group", Boeva = "#7C0A02", legend = TRUE, cex.title = 0.8, rot.title = 1, title.width = 2)
# plotTracks(list(dataTrack_regions, dataTrack_regions_sig, dataTrack_windows, dataTrack_windows_sig, dataTrack_bins, dataTrack_bins_sig, ATACtrack, H3K27ACtrack, PHOX2Btrack, HAND2track, gtrack, superenhancer,geneTrack), from = 41723346, to = 41774745, sizes = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 4), groupAnnotation = "group", Boeva = "#7C0A02", legend = TRUE, cex.title = 0.8, rot.title = 1, title.width = 2)
plotTracks(list(dataTrack_regions, dataTrack_regions_sig,dataTrack_windows, dataTrack_windows_sig,dataTrack_bins, dataTrack_bins_sig, ATACtrack, H3K27ACtrack, PHOX2Btrack, HAND2track, gtrack, superenhancer, hl_phox),  from = 41723346, to = 41774745, sizes = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,2), groupAnnotation = "group", Boeva = "#7C0A02", transcriptAnnotation = "symbol",legend = TRUE, cex.title = 0.5, rot.title = 1,collapseTranscripts=TRUE, title.width = 2)
plotTracks(list(dataTrack_regions_k, dataTrack_regions_sig_k,dataTrack_windows_k, dataTrack_windows_sig_k,dataTrack_bins_k, dataTrack_bins_sig_k, ATACtrack, H3K27ACtrack, PHOX2Btrack, HAND2track, gtrack, superenhancer, hl_phox), from = 41723346, to = 41774745, sizes = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2), groupAnnotation = "group", Boeva = "#7C0A02", legend = TRUE, cex.title = 0.8, rot.title = 1,transcriptAnnotation = "symbol", collapseTranscripts=TRUE, title.width = 2)
dev.off()
