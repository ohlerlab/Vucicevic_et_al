# ============================================================================
# Track plots
# ============================================================================
# This script generates track plots for visualizing genomic data.
# It uses various R libraries including tidyverse, rtracklayer, AnnotationHub, Gviz, Rsamtools, and biomaRt.
# The script includes the following sections:
# 
# 1. Loading required libraries
# 2. Defining color schemes
# 3. Setting working directory
# 4. Loading results from CSV files
# 5. Importing annotation tracks from BigWig files
# 6. Creating DataTrack objects for visualization
# 7. Defining fixed parameters for genomic visualization
# 8. Loading and processing regions data
# 9. Creating DataTrack objects for VPR and KRAB constructs
# 10. Plotting the tracks and saving the plots as PDF files
# 
# The script generates two main plots:
# 1. A plot showing the distance score for both VPR and KRAB constructs across a specified genomic region.
# 2. A plot showing the promoter region with bins for both VPR and KRAB constructs.
# 
# The plots include various annotation tracks such as ATAC-seq, H3K27ac, PHOX2B, HAND2, and superenhancers.
# The final plots are saved as PDF files in the specified working directory.

## Report Plotting Test
library(tidyverse)
library(rtracklayer)
library(AnnotationHub)
library(Gviz)
library(Rsamtools)
library(biomaRt)

# ==== colors =====
cold_colors <- c(
  "#003399",
  "#00528A",
  "#007877",
  "#009966"
)
my_colors_tracks <- c("#013766", "white", "#ac0e28")
names(my_colors_tracks) <- c("down", "nc", "up")

# ==== parameters =====
setwd("~/Documents/Phox2b/PAPER/visualization/")
# ============================================================================
# Loading results
# ============================================================================
# ==== Guides =====
VPR_guide_results_plot <- read_csv("VPR_guide_results_plot.csv")
KRAB_guide_results_plot <- read_csv("KRAB_guide_results_plot.csv")
# ==== Bins =====
VPR_window_results_plot <- read_csv("VPR_window_results_plot.csv")
KRAB_window_results_plot <- read_csv("KRAB_window_results_plot.csv")
# ==== Regions =====
VPR_results_plot <- read_csv("VPR_regions_results_plot.csv")
KRAB_results_plot <- read_csv("KRAB_regions_results_plot.csv")

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
geneTrack <- BiomartGeneRegionTrack(chromosome = chr, genome = gen, biomart = biomart, name = "GeneTrack", background.title = "black")

# ==== regions =====
regions <- readr::read_delim("~/Documents/Phox2b/merged_windows/merged_windowed_500.bed", delim = "\t", col_names = c("Chr", "start", "end", "bins"))
regions <- dplyr::mutate(regions, name = paste("Region", seq(1, nrow(regions), 1), sep = "_")) %>%
  dplyr::mutate(BINS = str_split(bins, pattern = ",")) %>%
  tidyr::unnest(cols = c(BINS))
regions_gviz <- regions %>%
  dplyr::select(Chr, start, end, name) %>%
  dplyr::distinct() %>%
  dplyr::full_join(VPR_results_plot, by = c("name" = "Gene")) %>%
  tidyr::replace_na(list(pvalue = 1, fdr = 1, effect = "fixed", estimate = 0))
regions_gviz_k <- regions %>%
  dplyr::select(Chr, start, end, name) %>%
  distinct() %>%
  full_join(KRAB_results_plot, by = c("name" = "Gene")) %>%
  tidyr::replace_na(list(pvalue = 1, fdr = 1, effect = "fixed", estimate = 0))




# regions
scores <- -log10(regions_gviz$fdr)[-944]
significant <- as.numeric(regions_gviz$fdr < 0.05)[-944]
coeffs <- regions_gviz$estimate[-944]

# Distance plot VPR
df <- regions_gviz %>%
  mutate("score" = -log10(fdr)) %>%
  mutate("color" = ifelse(fdr < 0.05, ifelse(estimate < 0, "#1683D0", "#d60000"), "#D2D2D2FF")) %>%
  mutate("sig" = as.numeric(fdr < 0.05)) %>%
  drop_na("name")

annot <- makeGRangesFromDataFrame(regions_gviz[-944, ], seqinfo = info)

dataTrack_ns <- DataTrack(
  range = annot,
  data = df$score * sign(df$estimate),
  name = "        VPR",
  type = "hist",
  ylim = c(-60, 10),
  col.histogram = "#D2D2D2FF",
  fill.histogram = "#D2D2D2FF"
)

dataTrack_down <- DataTrack(
  range = annot[df$color == "#1683D0"],
  data = -df$score[df$color == "#1683D0"],
  name = "        VPR",
  ylim = c(-60, 10),
  type = "hist",
  col.histogram = "#1683D0",
  fill.histogram = "#1683D0"
)

dataTrack_up <- DataTrack(
  range = annot[df$color == "#d60000"],
  data = df$score[df$color == "#d60000"],
  name = "        VPR",
  ylim = c(-60, 10),
  type = "hist",
  col.histogram = "#d60000",
  fill.histogram = "#d60000"
)

# dataTrack <- OverlayTrack(trackList = list(dataTrack_down,dataTrack_up),background.title='#e76f3d')
dataTrack <- OverlayTrack(trackList = list(dataTrack_ns, dataTrack_down, dataTrack_up), background.title = "#e76f3d")

# KRAB
# Distance plot KRAB
df_k <- regions_gviz_k %>%
  mutate("score" = -log10(fdr)) %>%
  mutate("color" = ifelse(fdr < 0.05, ifelse(estimate < 0, "#1683D0", "#d60000"), "#D2D2D2FF")) %>%
  drop_na("name")

dataTrack_ns_k <- DataTrack(
  range = annot,
  data = df_k$score * sign(df$estimate),
  name = "        KRAB",
  type = "hist",
  ylim = c(-60, 10),
  col.histogram = "#D2D2D2FF",
  fill.histogram = "#D2D2D2FF"
)

dataTrack_down_k <- DataTrack(
  range = annot[df_k$color == "#1683D0"],
  data = -df_k$score[df_k$color == "#1683D0"],
  name = "        KRAB",
  ylim = c(-60, 10),
  type = "hist",
  col.histogram = "#1683D0",
  fill.histogram = "#1683D0"
)

dataTrack_up_k <- DataTrack(
  range = annot[df_k$color == "#d60000"],
  data = df_k$score[df_k$color == "#d60000"],
  name = "        KRAB",
  ylim = c(-60, 10),
  type = "hist",
  col.histogram = "#d60000",
  fill.histogram = "#d60000"
)

# dataTrack_k <- OverlayTrack(trackList = list(dataTrack_down_k,dataTrack_up_k),background.title='#1683D0')
dataTrack_k <- OverlayTrack(trackList = list(dataTrack_ns_k, dataTrack_down_k, dataTrack_up_k), background.title = "#1683D0")

# mini scale
dataTrack_ns_k_s <- DataTrack(
  range = annot,
  data = df_k$score * df$estimate,
  name = "        KRAB",
  ylim = c(0, 3),
  type = "hist",
  col.histogram = "#d60000",
  fill.histogram = "#d60000"
)

dataTrack_down_k_s <- DataTrack(
  range = annot[df_k$color == "#1683D0"],
  data = df_k$score[df_k$color == "#1683D0"],
  name = "        KRAB",
  ylim = c(0, 3),
  type = "hist",
  col.histogram = "#1683D0",
  fill.histogram = "#1683D0"
)

dataTrack_up_k_s <- DataTrack(
  range = annot[df_k$color == "#d60000"],
  data = df_k$score[df_k$color == "#d60000"],
  name = "        KRAB",
  ylim = c(0, 3),
  type = "hist",
  col.histogram = "#d60000",
  fill.histogram = "#d60000"
)

dataTrack_k_s <- OverlayTrack(trackList = list(dataTrack_ns_k, dataTrack_down_k, dataTrack_up_k), background.title = "#168fd0")
# additional data
geneTrack <- BiomartGeneRegionTrack(
  chromosome = chr,
  genome = gen,
  biomart = biomart,
  name = "GeneTrack",
  background.title = "black",
  start = 40749815, end = 42748835,
  col.line = NULL,
  col = NULL,
  transcriptAnnotation = "symbol",
  collapseTranscripts = "longest"
)

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

# plotTracks(list(gtrack,dataTrack,geneTrack),transcriptAnnotation = "symbol",from = 40749815, to = 42748835, collapseTranscripts="longest", col.line = NULL)
pdf("both_constructs_distance_score.pdf", width = 10, height = 5)
# Whole TAD
# plotTracks(list(dataTrack, ATACtrack, H3K27ACtrack, PHOX2Btrack, HAND2track, gtrack, superenhancer, hl_phox), from = 40749815, to = 42748835, sizes = c(4, 1, 1, 1, 1, 1, 1, 2), groupAnnotation = "group", Boeva = "#7C0A02", legend = TRUE, cex.title = 0.8, rot.title = 1, title.width = 2)
# plotTracks(list(dataTrack_k, ATACtrack, H3K27ACtrack, PHOX2Btrack, HAND2track, gtrack, superenhancer, hl_phox), from = 40749815, to = 42748835, sizes = c(4, 1, 1, 1, 1, 1, 1, 2), groupAnnotation = "group", Boeva = "#7C0A02", legend = TRUE, cex.title = 0.8, rot.title = 1, title.width = 2)
plotTracks(list(dataTrack, dataTrack_k, ATACtrack, H3K27ACtrack, PHOX2Btrack, HAND2track, gtrack, superenhancer, hl_phox), from = 40749815, to = 42748835, sizes = c(4, 4, 1, 1, 1, 1, 1, 1, 2), groupAnnotation = "group", Boeva = "#7C0A02", legend = TRUE, cex.title = 0.8, rot.title = 1, title.width = 2)
dev.off()


# Promoter region with bins
df <- VPR_bin_results %>%
  mutate("score" = -log10(fdr)) %>%
  mutate("color" = ifelse(fdr < 0.05, ifelse(coeff < 0, "#1683D0", "#d60000"), "#D2D2D2FF")) %>%
  mutate("sig" = as.numeric(fdr < 0.05)) %>%
  filter(!is.na(fdr))
annot <- makeGRangesFromDataFrame(df, seqinfo = info, start.field = "start.bin", end.field = "end.bin", seqnames.field = "Chr.bin")

# ==== VPR =====
dataTrack_ns <- DataTrack(
  range = annot,
  data = df$score * sign(df$coeff),
  name = "        VPR",
  type = "hist",
  ylim = c(-8, 8),
  col.histogram = "#D2D2D2FF",
  fill.histogram = "#D2D2D2FF"
)

dataTrack_down <- DataTrack(
  range = annot[df$color == "#1683D0"],
  data = -df$score[df$color == "#1683D0"],
  name = "        VPR",
  ylim = c(-8, 8),
  type = "hist",
  col.histogram = "#1683D0",
  fill.histogram = "#1683D0"
)

dataTrack_up <- DataTrack(
  range = annot[df$color == "#d60000"],
  data = df$score[df$color == "#d60000"],
  name = "        VPR",
  ylim = c(-8, 8),
  type = "hist",
  col.histogram = "#d60000",
  fill.histogram = "#d60000"
)

# dataTrack <- OverlayTrack(trackList = list( dataTrack_down,dataTrack_up),background.title='#e76f3d')
dataTrack <- OverlayTrack(trackList = list(dataTrack_ns, dataTrack_down, dataTrack_up), background.title = "#e76f3d")

# KRAB
df_k <- KRAB_bin_results %>%
  mutate("score" = -log10(fdr)) %>%
  mutate("color" = ifelse(fdr < 0.05, ifelse(coeff < 0, "#1683D0", "#d60000"), "#D2D2D2FF")) %>%
  mutate("sig" = as.numeric(fdr < 0.05))
dataTrack_k <- DataTrack(
  range = annot,
  data = df_k$coeff,
  name = "    KRAB",
  ylim = c(-0.04, 0.04),
  type = "heatmap",
  gradient = c("#1683D0", "white", "#d60000"),
  background.title = "#1683D0"
)
dataTrack2_k <- DataTrack(
  range = annot,
  data = df_k$score,
  name = "    KRAB",
  ylim = c(0, 8),
  type = "heatmap",
  gradient = c("white", "black"),
  background.title = "#1683D0"
)
# dataTrack_k <- OverlayTrack(trackList = list(dataTrack_ns_k,dataTrack_down_k,dataTrack_up_k),background.title='#1683D0')
pdf("both_constructs_distance_score_promoter.pdf", width = 10, height = 5)
# plotTracks(list(dataTrack, ATACtrack, H3K27ACtrack, PHOX2Btrack, HAND2track, gtrack, superenhancer, hl_phox), from = 41734346, to = 41763745, sizes = c(4, 1, 1, 1, 1, 1, 1, 2), groupAnnotation = "group", Boeva = "#7C0A02", legend = TRUE, cex.title = 0.8, rot.title = 1, title.width = 2)
# plotTracks(list(dataTrack_k, ATACtrack, H3K27ACtrack, PHOX2Btrack, HAND2track, gtrack, superenhancer, hl_phox), from = 41734346, to = 41763745, sizes = c(4, 1, 1, 1, 1, 1, 1, 2), groupAnnotation = "group", Boeva = "#7C0A02", legend = TRUE, cex.title = 0.8, rot.title = 1, title.width = 2)
plotTracks(list(dataTrack, dataTrack_k, ATACtrack, H3K27ACtrack, PHOX2Btrack, HAND2track, gtrack, superenhancer, hl_phox), from = 41734346, to = 41763745, sizes = c(4, 4, 1, 1, 1, 1, 1, 1, 2), groupAnnotation = "group", Boeva = "#7C0A02", legend = TRUE, cex.title = 0.8, rot.title = 1, title.width = 2)
dev.off()
