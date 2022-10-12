#' Function to plot 4x4 grid of APAs
#' in a super specific order
#' 
#' Values are normalized to the maximum
#' center pixel of each column.
#' 
#' @param apas List of apa matrices plotted in column order
#' 
#' @returns 4x4 grid of apa plots
#' @noRd
plotApas <- function(apas = list(aggGainedCont, aggGainedSorb,
                                 aggLostCont, aggLostSorb)) {
  
  ## Create page
  pageCreate(width = 4, height = 4, showGuides = FALSE)
  
  ## Define common params
  p <- pgParams(
    width = 1.5,
    height = 1.5,
    just = c('left', 'top'),
    palette = colorRampPalette(brewer.pal(9, 'YlGnBu')),
    space = 0.1
  )
  
  ## Define positions for placing plots
  rowPos <- layoutRow(y=0.5, h=p$height, s=p$space, n=2)
  colPos <- layoutCol(x=0.5, w=p$width, s=p$space, n=2)
  
  ## Normalize to center pixels of columns
  centerPixels <- lapply(apas, \(x) x[buffer+1, buffer+1]) |> unlist()
  zCol1 <- c(0, max(centerPixels[1:2]))
  zCol2 <- c(0, max(centerPixels[3:4]))
  zranges <- list(zCol1, zCol1, zCol2, zCol2)
  
  ## Plot APAs & annotate
  plots <- 
    pmap(.l = list(apas,
                   rep(rowPos, each=2),
                   rep(colPos, 2),
                   zranges),
         .f = \(apa, x, y, z) {
           plotApa(apa = apa,
                   x = x,
                   y = y,
                   zrange = z,
                   params=p)
         })
  
  ## Annotate rows
  plotText(label = c("Control", "Sorb"),
           x = colPos[1] - p$space/2,
           y = rowPos + p$height/2,
           rot = 90,
           just = c("center", "bottom"))
  
  ## Annotate columns
  plotText(label = c("Gained Loops", "Lost Loops"),
           x = colPos + p$width/2,
           y = rowPos[1] - p$space/2,
           just = c("center", "bottom"))
  
  ## Annotate heatmap legends
  annoHeatmapLegend(plot = plots[[2]],
                    orientation = 'h',
                    x = colPos[1],
                    y = rowPos[2] + p$height + p$space,
                    fontcolor = 'black',
                    width = p$width,
                    height = p$space)
  
  annoHeatmapLegend(plot = plots[[4]],
                    orientation = 'h',
                    x = colPos[2],
                    y = rowPos[2] + p$height + p$space,
                    fontcolor = 'black',
                    width = p$width,
                    height = p$space)
  
  
}