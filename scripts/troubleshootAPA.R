## Generate APA matrices to look for outlier
## interactions.

## Code to install the development
## version of mariner:
if (!require("mariner", quietly = TRUE)) {
  remotes::install_github("EricSDavis/mariner@dev", force = TRUE)
}

## Install hictoolsr but exclude conflicting functions
if (!require("hictoolsr", quietly = TRUE,
             exclude=c("as_ginteractions",
                       "makeGInteractionsFromDataFrame"))) {
  remotes::install_github("EricSDavis/hictoolsr@1.1.2", force = TRUE)
}


## Load required libraries
library(InteractionSet)
library(mariner)
library(plotgardener)
library(RColorBrewer)
library(purrr)

## Load utility functions
source("scripts/utils/findBadLoops.R")
source("scripts/utils/layoutFunctions.R")
source("scripts/utils/plotApas.R")

## Define hic files
hicFiles <- list.files("data/raw/hic/hicFiles",
                       full.names = TRUE)

## Load differential loop calls
allLoops <- readRDS("data/raw/hic/loops/YAPP_hic_diff_loopCounts.rds")

## Separate into gained/lost loops
gainedLoops <- 
  allLoops |> 
  subset(padj < 0.1 & log2FoldChange > 0) |>
  reduceRegions()

lostLoops <-
  allLoops |> 
  subset(padj < 0.1 & log2FoldChange < 0) |>
  reduceRegions()

## Define variables
res <- 10e03
buffer <- 10
norm <- "SCALE"

## Filter out short loops
gainedLoops <- filterBedpe(gainedLoops, res = res, buffer = buffer)
lostLoops <- filterBedpe(lostLoops, res=res, buffer = buffer)

## Expand pixels to matrices and extract
gainedIset <- 
  pixelsToMatrices(gainedLoops, buffer) |>
  pullHicMatrices(binSize = res,
                  files = hicFiles,
                  norm = norm)
lostIset <- 
  pixelsToMatrices(lostLoops, buffer) |>
  pullHicMatrices(binSize = res,
                  files = hicFiles,
                  norm = norm)

## Aggregate matrices (before cleaning)
aggGainedCont <- 
  assay(gainedIset) |>
  aperm(c(3,4,1,2)) |>
  {\(x) apply(x[,,,1], c(1,2), sum)}()

aggLostCont <- 
  assay(lostIset) |>
  aperm(c(3,4,1,2)) |>
  {\(x) apply(x[,,,1], c(1,2), sum)}()

aggGainedSorb <- 
  assay(gainedIset) |>
  aperm(c(3,4,1,2)) |>
  {\(x) apply(x[,,,2], c(1,2), sum)}()

aggLostSorb <- 
  assay(lostIset) |>
  aperm(c(3,4,1,2)) |>
  {\(x) apply(x[,,,2], c(1,2), sum)}()


## Remove loops with missing regions -------------------------------------------

## Indices of loops with missing values
## for gained loops
gainedBadIndices <- 
  gainedIset |>
  assay() |>
  aperm(c(3,4,1,2)) |>
  {\(x) {
    cont <- findBadLoops(x[,,,1])
    sorb <- findBadLoops(x[,,,2])
    sort(unique(c(cont, sorb)))
  }}() 

## Indices of loops with missing values
## for lost loops
lostBadIndices <- 
  lostIset |>
  assay() |>
  aperm(c(3,4,1,2)) |>
  {\(x) {
    cont <- findBadLoops(x[,,,1])
    sorb <- findBadLoops(x[,,,2])
    sort(unique(c(cont, sorb)))
  }}() 

## Define which loops to accept by removing
## bad indices
gainedAccept <- which(!1:length(gainedIset) %in% gainedBadIndices)
lostAccept <- which(!1:length(lostIset) %in% lostBadIndices)


## Aggregate matrices (after cleaning)
aggGainedContClean <- 
  assay(gainedIset[gainedAccept]) |>
  aperm(c(3,4,1,2)) |>
  {\(x) apply(x[,,,1], c(1,2), sum)}()

aggLostContClean <- 
  assay(lostIset[lostAccept]) |>
  aperm(c(3,4,1,2)) |>
  {\(x) apply(x[,,,1], c(1,2), sum)}()

aggGainedSorbClean <- 
  assay(gainedIset[gainedAccept]) |>
  aperm(c(3,4,1,2)) |>
  {\(x) apply(x[,,,2], c(1,2), sum)}()

aggLostSorbClean <- 
  assay(lostIset[lostAccept]) |>
  aperm(c(3,4,1,2)) |>
  {\(x) apply(x[,,,2], c(1,2), sum)}()


## Visualize APA with plotgardener ---------------------------------------------

## Before
pdf("plots/troubleshootAPA_before.pdf", width = 4, height = 4)
plotApas(apas = list(aggGainedCont, aggGainedSorb,
                     aggLostCont, aggLostSorb))
dev.off()

## After
pdf("plots/troubleshootAPA_after.pdf", width = 4, height = 4)
plotApas(apas = list(aggGainedContClean, aggGainedSorbClean,
                     aggLostContClean, aggLostSorbClean))
dev.off()



