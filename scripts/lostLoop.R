
## Lost differential loop visualizations (HiC rectangles)
## NCI Phase Separation Talk

# Load data ---------------------------------------------------------------
library(plotgardener)
library(hictoolsr)
library(InteractionSet)
library(tidyverse)
library(glue)
library(dbscan)
library(data.table)
library(org.Hs.eg.db)

diff_loopCounts <- readRDS("data/raw/hic/loops/YAPP_hic_diff_loopCounts.rds")

# Setting lost loops -------------------------------------------

lost_adj <- subset(diff_loopCounts, padj < 0.1 & log2FoldChange < 0)

## filter for the best lost loops
bestLost <- head(lost_adj[order(lost_adj$log2FoldChange, decreasing = T)],100)
bestLost <- head(lost_adj[order(lost_adj$pvalue, decreasing = F)], 100)

# # all lost loops
# loopRegions <-
#   GRanges(seqnames = as.character(seqnames(anchors(x = lost, "first"))),
#           ranges = IRanges(start = start(anchors(lost, 'first')),
#                            end = end(anchors(lost, 'second'))))

# top 50 lost loops
# bestLost <- swapAnchors(bestLost)

loopRegions_lost <- 
  GRanges(seqnames = as.character(seqnames(anchors(x = bestLost, "first"))),
          ranges = IRanges(start = start(anchors(bestLost, 'first')),
                           end = end(anchors(bestLost, 'second'))))

## Expand regions by buffer
buffer <- 200e3
loopRegions_lost <- loopRegions_lost + buffer
loopRegions_lost <- as.data.frame(loopRegions_lost)

## Use tidyverse to remove `chr`in seqnames
loopRegions_lost <- loopRegions_lost |>
  mutate(seqnames = str_remove(seqnames, "chr"))

# load loop lists --------------------------------------------------------

# cont_loops <- readRDS("data/processed/hic/cont_bothDroso_loops.rds")
# 
# sorb_loops <- readRDS("data/processed/hic/sorb_bothDroso_loops.rds")
# 
# omega_loops <- readRDS("data/processed/hic/omega_bothDroso_loops.rds")

# Create Survey Plots -----------------------------------------------------

##make pdf
pdf(file = "plots/lostLoop.pdf",
    width = 5.75,
    height = 6)

## Loop through each region

## Define parameters
p <- pgParams(assembly = "hg38",
              resolution = 10e3,
              chrom = paste0("chr",loopRegions_lost$seqnames[16]),
              chromstart = loopRegions_lost$start[16],
              chromend = loopRegions_lost$end[16],
              zrange = c(0,100),
              norm = "SCALE",
              x = 0.25,
              width = 5,
              length = 5,
              height = 2,
              fill = "#37a7db",
              linecolor = "#37a7db")


# Begin Visualization -----------------------------------------------------
## Make page
pageCreate(width = 5.75, height = 6,
           xgrid = 0, ygrid = 0, showGuides = T)

## Plot middle Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls 

control <- plotHicRectangle(data = "data/raw/hic/hicFiles/YAPP/YAPP_HEK_control_inter_30.hic",
                            params = p,
                            y = 0.5)

# annoHeatmapLegend(control, orientation = "v",
#                   fontsize = 8,
#                   fontcolor = "black",
#                   digits = 2,
#                   x = 5.5,
#                   y = 0.5,
#                   width = 0.1,
#                   height = 1.5,
#                   just = c("left", "top"),
#                   default.units = "inches")
# 
# annoPixels(control,
#            data = cont_loops,
#            shift = 0.5,
#            type = "arrow",
#            col = "#005AB5")
# 
## Plot bottom Hi-C rectangle + SIP `-isDroso = TRUE` & `-isDroso = TRUE` calls

sorb <- 
  plotHicRectangle(data = "data/raw/hic/hicFiles/YAPP/YAPP_HEK_sorbitol_inter_30.hic", 
                   params = p,
                   y = 2.6)
# 
# annoHeatmapLegend(sorb, orientation = "v",
#                   fontsize = 8,
#                   fontcolor = "black",
#                   digits = 2,
#                   x = 5.5,
#                   y = 2.6,
#                   width = 0.1,
#                   height = 1.5,
#                   just = c("left", "top"),
#                   default.units = "inches")
# 
# annoPixels(sorb, data = sorb_loops,
#            shift = 0.5,
#            type = "arrow",
#            col = "#DC3220")
# 
## Plot genes
plotGenes(param = p,
          chrom = p$chrom,
          x = 0.25,
          y = 4.7,
          height = 0.5)

## Plot genome label
plotGenomeLabel(params = p,
                x = 0.25,
                y = 5.3)

# Annotate Hi-C rectangles by treatment ------------------------------------

plotText(label = "untreated",
         x = 0.25,
         y = 0.5,
         just = c("top", "left"))

plotText(label = "+ sorbitol",
         x = 0.25,
         y = 2.6,
         just = c("top", "left"))

dev.off()