# Visualization Scripts for Phox2B project

# Author: Lorena Sofia Lopez Zepeda``
# Version: 2021-11-21

# Packages
library(tidyverse)

# Parameters
setwd("~/Documents/Phox2b/PAPER/visualization/")
# ============================================================================
# Loading results
# ============================================================================
# ==== Guides =====
VPR_guide_results_plot <- read_csv("VPR_guide_results_plot.csv")
KRAB_guide_results_plot <- read_csv("KRAB_guide_results_plot.csv")
# ==== Windows =====
VPR_window_results_plot <- read_csv("VPR_window_results_plot.csv")
KRAB_window_results_plot <- read_csv("KRAB_window_results_plot.csv")
# ==== Bins =====
VPR_bin_results_plot <- read_csv("CT_Results_input_VPR.txt")
KRAB_bin_results_plot <- read_csv("CT_Results_input_KRAB.txt")
# ==== Regions =====
VPR_results_plot <- read_csv("VPR_regions_results_plot.csv")
KRAB_results_plot <- read_csv("KRAB_regions_results_plot.csv")
# ============================================================================
# Plotting functions and values
# ============================================================================
# establishing theme and colors for the plots
mycolors <- c("#6B8993FF", "#E7695DFF", "#D2D2D2FF")
mycolorstable <- c("#6B8993FF", "#E7695DFF", "#969696FF")
names(mycolors) <- c("Depleted", "Enriched", "NoChange")
names(mycolorstable) <- c("Depleted", "Enriched", "NoChange")
theme_dusa <- function() {
  ggplot2::theme(
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent"),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent")
  )
}
# ==== phox_volcano_plot =====
phox_volcano_plot <- function(df) {
  ggplot(df, aes(x = estimate, y = -log10(fdr), col = Result)) +
    geom_point(alpha = 0.70) +
    theme_dusa() +
    scale_color_manual(values = mycolors)
}
phox_volcano_plot_guides <- function(df) {
  ggplot(df, aes(x = logfold, y = -log10(fdr), col = Result)) +
    geom_jitter(alpha = 0.5, width = 0.25, height = 0.25) +
    theme_dusa() +
    scale_color_manual(values = mycolors)
}
# ==== tabulation plots =====
tab_plot <- function(df) {
  ggplot(df, aes(VPR, KRAB)) +
    geom_point(aes(size = Freq), color = "grey") +
    theme_classic() +
    scale_size_continuous(range = c(10, 30)) +
    geom_text(aes(label = Freq)) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18)
    )
}
# ============================================================================
# Regions plot
# ============================================================================

pdf('region_images.pdf',width = 6, height = 4)
# ==== VPR =====
v_r <- phox_volcano_plot(VPR_results_plot)
v_r +
  xlim(-0.075, 0.075) +
  ylim(0, 60) +
  ggtitle("VPR regions")

# ==== KRAB =====
k_r <- phox_volcano_plot(KRAB_results_plot)
k_r +
  xlim(-0.075, 0.075) +
  ylim(0, 60) +
  ggtitle("KRAB regions")

# ==== Both =====
full_regions_results_plot <- rbind(VPR_results_plot, KRAB_results_plot) %>%
  drop_na(Gene)
b_r <- phox_volcano_plot(full_regions_results_plot)
b_r +
  xlim(-0.075, 0.075) +
  ylim(0, 60) +
  scale_color_manual(values = mycolors) +
  facet_wrap(~Construct) +
  ggtitle("Region Results by Construct")

# ==== Cross tabulation =====
library(GGally)
ggtable(full_regions_results_plot, "Result", "Construct", title = "Regions results by construct")
ggtable(full_regions_results_plot, "Result", "Construct", cells = "row.prop", title = "Regions results by construct (proportion)", mapping = ggplot2::aes(color = Result, )) +
  scale_color_manual(values = mycolorstable)
dev.off()

# ============================================================================
# Bins
# ============================================================================

pdf('window_images.pdf',width = 6, height = 4)

# ==== VPR =====
v_w <- phox_volcano_plot(VPR_window_results_plot)
v_w +
  xlim(-0.075, 0.075) +
  ylim(0, 8) +
  ggtitle("VPR bins")

# ==== KRAB =====
k_w <- phox_volcano_plot(KRAB_window_results_plot)
k_w +
  xlim(-0.075, 0.075) +
  ylim(0, 8) +
  ggtitle("KRAB bins")

# ==== Both =====
full_window_results_plot <- rbind(VPR_window_results_plot, KRAB_window_results_plot)
b_w <- phox_volcano_plot(full_window_results_plot)
b_w+
xlim(-0.075, 0.075) +
  ylim(0, 8) +
  facet_wrap(~Construct) +
  ggtitle("Window Results by Construct")

# ==== tabulation =====
ggtable(full_window_results_plot, "Result", "Construct", title = "Window results by construct")
ggtable(full_window_results_plot, "Result", "Construct", cells = "row.prop", title = "Window results by construct (proportion)", mapping = ggplot2::aes(color = Result, )) +
  scale_color_manual(values = mycolorstable)
dev.off()

# ==== Bin/construct table =====
VPR_hits <- VPR_window_results_plot %>% filter(Result != "NoChange")
KRAB_hits <- KRAB_window_results_plot %>% filter(Result != "NoChange")
shared_windows <- inner_join(VPR_hits, KRAB_hits, by = "bin", suffix = c(".VPR", ".KRAB")) %>%
  dplyr::select(bin, Result.VPR, Result.KRAB)

# ==== Shared bins in constructs =====
shared_windows_table <- shared_windows %>%
  dplyr::mutate("VPR_results" = ifelse(Result.VPR == "Depleted", -1, 1)) %>%
  dplyr::mutate("KRAB_results" = ifelse(Result.KRAB == "Depleted", -1, 1)) %>%
  mutate("Similar" = (sign(VPR_results) == sign(KRAB_results))) %>%
  dplyr::select(bin, VPR_results, KRAB_results, Similar)

# ==== Shared bins plot =====
shared_bins <- as.data.frame(table(shared_windows$Result.VPR, shared_windows$Result.KRAB))
colnames(shared_bins) <- c("VPR", "KRAB", "Freq")
g <- tab_plot(shared_bins)
pdf("significant_bins_both.pdf")
g
dev.off()

# ============================================================================
# Guides
# ============================================================================

pdf("guide_images.pdf", width = 6, height = 4)

# ==== VPR =====
v_g <- phox_volcano_plot_guides(VPR_guide_results_plot)
v_g +
  xlim(-7.5, 7.5) +
  ylim(0, 20) +
  ggtitle("VPR guides (FDR)")

# ==== KRAB =====
k_g <- phox_volcano_plot_guides(KRAB_guide_results_plot)
k_g +
  xlim(-7.5, 7.5) +
  ylim(0, 20) +
  ggtitle("KRAB guides (FDR)")

# ==== Both =====
VPR_guide_results_plot$Construct <- rep("VPR", nrow(VPR_guide_results_plot))
KRAB_guide_results_plot$Construct <- rep("KRAB", nrow(KRAB_guide_results_plot))
full_guide_results_plot <- rbind(VPR_guide_results_plot, KRAB_guide_results_plot)
b_g <- phox_volcano_plot_guides(full_guide_results_plot)
b_g +
  xlim(-7.5, 7.5) +
  ylim(0, 20) +
  facet_wrap(~Construct) +
  ggtitle("Guide Results by Construct")

# ==== tabulation =====
ggtable(full_guide_results_plot, "Result", "Construct", title = "Guide results by construct")
ggtable(full_guide_results_plot, "Result", "Construct", cells = "row.prop", title = "Guide results by construct (proportion)", mapping = ggplot2::aes(color = Result, )) +
  scale_color_manual(values = mycolorstable)

dev.off()


# GUIDE/construct table
VPR_hits <- VPR_guide_results_plot %>% filter(Result != "NoChange")
KRAB_hits <- KRAB_guide_results_plot %>% filter(Result != "NoChange")
shared_windows <- inner_join(VPR_hits, KRAB_hits, by = "GUIDE", suffix = c(".VPR", ".KRAB")) %>%
  dplyr::select(GUIDE, Result.VPR, Result.KRAB)

# Shared GUIDE in constructs
shared_windows_table <- shared_windows %>%
  dplyr::mutate("VPR_results" = ifelse(Result.VPR == "Depleted", -1, 1)) %>%
  dplyr::mutate("KRAB_results" = ifelse(Result.KRAB == "Depleted", -1, 1)) %>%
  mutate("Similar" = (sign(VPR_results) == sign(KRAB_results))) %>%
  dplyr::select(GUIDE, VPR_results, KRAB_results, Similar)

# ==== tabulation plot =====
shared_guides <- as.data.frame(table(shared_windows$Result.VPR, shared_windows$Result.KRAB))
colnames(test) <- c("VPR", "KRAB", "Freq")
g <- tab_plot(shared_guides)
pdf("significant_GUIDE.pdf")
g
dev.off()
