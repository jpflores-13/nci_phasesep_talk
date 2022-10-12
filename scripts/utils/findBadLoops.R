## Identify mats with all 0s in
## any rows or columns
findBadLoops <- function(a, thresh = 1) {
  badRows <- 
    apply(a, 3, rowSums) |>
    apply(2, min) |>
    {\(x) which(x < thresh)}()
  
  badCols <- 
    apply(a, 3, colSums) |>
    apply(2, min) |>
    {\(x) which(x < thresh)}()
  
  badIndices <-sort(c(badRows, badCols))
  return(badIndices)
}