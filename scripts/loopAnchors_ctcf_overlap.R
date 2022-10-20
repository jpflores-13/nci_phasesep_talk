## Make pie charts side by side showing % of loop anchors that overlap a CTCF peak 
## (control on left, gained on right)

##What percentage of loop anchors (from all loops identified in untreated cells) overlap a CTCF peak?
##What percentage of loop anchors (from gained loops) overlap a CTCF peak?

library(tidyverse)
library(mariner)
library(plyranges)
library(InteractionSet)
library(purrr)
library(plotgardener)
library(RColorBrewer)
library(nullranges)

## Read in and prune narrowPeak files
peakFiles <- list.files(path = "data/raw/chip", full.names = TRUE)
chipPeaks <- 
  read_narrowpeaks(peakFiles) |>
  keepStandardChromosomes(pruning.mode = 'coarse') |> 
  select(peak)

## Read in loops 
all_loops <- readRDS("data/raw/hic/loops/YAPP_hic_diff_loopCounts_noDroso.rds") |> 
  reduceRegions() 

## create metadata columns for contact frequency & size 
mcols(all_loops)$loop_size <- pairdist(all_loops)
mcols(all_loops)$loop_type <- case_when(
  mcols(all_loops)$padj < 0.05 & mcols(all_loops)$log2FoldChange > 1 ~ "gained",
  # mcols(all_loops)$padj < 0.1 & mcols(all_loops)$log2FoldChange > 0 ~ "gained",
  mcols(all_loops)$padj < 0.05 & mcols(all_loops)$log2FoldChange < -1 ~ "lost",
  mcols(all_loops)$padj > 0.05 ~ "static",
  is.character("NA") ~ "other")

## contact frequency
mcols(all_loops)$sorb_contacts <- mcols(all_loops)$YAPP_HEK_sorbitol_4_2_inter_30.hic +
  mcols(all_loops)$YAPP_HEK_sorbitol_5_2_inter_30.hic + mcols(all_loops)$YAPP_HEK_sorbitol_6_2_inter_30.hic

## use matchRanges to select a null set of control sample loops that is matched for size & contact frequency
nullSet <- matchRanges(focal = all_loops[mcols(all_loops)$loop_type == "gained"],
                       pool = all_loops[!mcols(all_loops)$loop_type == "gained"],
                       covar = ~ loop_size + sorb_contacts, 
                       method = 'stratified',
                       replace = FALSE) |> 
  reduceRegions() |> 
  regions()

gained_loops <- all_loops |> 
  subset(padj < 0.05 & log2FoldChange > 1) |> 
  reduceRegions() |> 
  regions()

all_loops <- regions(all_loops)  

## returns only the ranges in the first object that have overlaps with any ranges in the second object.

chip_gained <- length(subsetByOverlaps(gained_loops, chipPeaks)) 

chip_nullSet <- length(subsetByOverlaps(nullSet, chipPeaks)) 

chip_control <- length(subsetByOverlaps(all_loops,chipPeaks))

## Calculate percentages of the total number of loops and round
gainedP <- signif(chip_gained/length(gained_loops)*100, 3)
nullSetP <- signif(chip_nullSet/length(nullSet)*100, 3)
controlP <- signif(chip_control/length(all_loops)*100, 3)

## Create data frame to plot results for gained YAPP loop %
df_yapp <- data.frame(
  group = c("overlaps_with_ctcfPeaks", "all_gained"),
  number = c(chip_gained, length(gained_loops)),
  percent = c(gainedP,(100 - gainedP))
)

## Plot results
library(ggplot2)
pdf("plots/yapp_ctcf_piechart.pdf", width = 4, height = 4)
ggplot(data = df_yapp, aes(x = 1, y = percent, fill = group))+
  geom_col(col = "white")+
  coord_polar("y") +
  scale_fill_manual(values = c("#2171b5", "#bdd7e7"))+
  labs(title = "% of loop anchors overlapping a CTCF peak") +
  # theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = -5, face = "bold"),
    legend.position = "none", text = element_text(face = "bold")
  ) +
  annotate(geom = "text",
           x = c(1.0, 1.0),
           y = c(20, 50),
           label = paste0(df_yapp$group, " (", df$df_yapp, "%)"),
           col = c("white", "white"),
           size = 3.5)
dev.off()

# ## Create data frame to plot results for nullSet loop %
# df_nullSet  <- data.frame(
#   group = c("YAPP","Matched", "CTCF"),
#   number = c(chip_gained, chip_nullSet, chip_control),
#   percent = c(gainedP, nullSetP, (controlP-(gainedP+nullSetP))))
# 
# ## Plot results
# library(ggplot2)
# pdf("plots/yapp_matched_ctcf_piechart.pdf", width = 4, height = 4)
# ggplot(data = df_nullSet, aes(x = 1, y = percent, fill = group))+
#   geom_col(col = "white")+
#   coord_polar("y") +
#   scale_fill_manual(values = c("#2171b5", "#bdd7e7", "#6baed6"))+
#   labs(title = "Loop anchors shared with CTCF") +
#   theme_void() +
#   theme(
#     plot.title = element_text(hjust = 0.5, vjust = -5, face = "bold"),
#     legend.position = "none", text = element_text(face = "bold")
#   ) +
#   annotate(geom = "text",
#            x = c(1, 1, 1.0),
#            y = c(20, 80, 60),
#            label = paste0(df_nullSet$group, " (", df_nullSet$percent, "%)"),
#            col = c("black", "black", "white"),
#            size = c(3, 3, 4.5))
# dev.off()

## Create data frame to plot results for gained YAPP loop %
df_control <- data.frame(
  group = c("overlaps_with_ctcfPeaks", "control"),
  number = c(chip_control, length(all_loops)),
  percent = c(controlP,(100 - controlP))
)

## Plot results
library(ggplot2)
pdf("plots/yapp_ctcf_piechart.pdf", width = 4, height = 4)
ggplot(data = df_control, aes(x = 1, y = percent, fill = group))+
  geom_col(col = "white")+
  coord_polar("y") +
  scale_fill_manual(values = c("#2171b5", "#bdd7e7"))+
  labs(title = "% of loop anchors overlapping a CTCF peak") +
  # theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, vjust = -5, face = "bold"),
    legend.position = "none", text = element_text(face = "bold")
  ) +
  annotate(geom = "text",
           x = c(1.0, 1.0),
           y = c(20, 50),
           label = paste0(df_control$group, " (", df$df_control, "%)"),
           col = c("white", "white"),
           size = 3.5)
dev.off()