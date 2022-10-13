
## Gained differential loop visualizations (HiC rectangles)
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

# Setting static and lost loops -------------------------------------------

## gained loops with a p-adj. value of < 0.1 and a (+) log2FC
gained_adj <- subset(diff_loopCounts, padj < 0.1 & log2FoldChange > 0)

## filter for the 100 best gained loops
bestGained <- head(gained_adj[order(gained_adj$log2FoldChange, decreasing = T)],100)
bestGained <- head(gained_adj[order(gained_adj$padj, decreasing = F)], 100)

# top 100 gained loops
## If the below function doesn't work, might need to use swapAchors()

loopRegions_gained <- 
  GRanges(seqnames = as.character(seqnames(anchors(x = bestGained, "first"))),
          ranges = IRanges(start = start(anchors(bestGained, "first")),
                           end = end(anchors(bestGained, "second"))))

## Expand regions by buffer
buffer <- 200e3
loopRegions_gained <- loopRegions_gained + buffer
loopRegions_gained <- as.data.frame(loopRegions_gained)

## Use tidyverse to remove `chr`in seqnames
loopRegions_gained <- loopRegions_gained |>
  mutate(seqnames = str_remove(seqnames, "chr"))

# load loop lists --------------------------------------------------------

# cont_loops <- readRDS("data/processed/hic/cont_bothDroso_loops.rds")
# 
# sorb_loops <- readRDS("data/processed/hic/sorb_bothDroso_loops.rds")
# 
# omega_loops <- readRDS("data/processed/hic/omega_bothDroso_loops.rds")

# Create Survey Plots -----------------------------------------------------

##make pdf
pdf(file = "plots/gainedLoop.pdf",
    width = 5.75,
    height = 6)

## Loop through each region
  
  ## Define parameters
  p <- pgParams(assembly = "hg38",
                resolution = 10e3,
                chrom = paste0("chr",loopRegions_gained$seqnames[22]),
                chromstart = loopRegions_gained$start[22],
                chromend = loopRegions_gained$end[22],
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