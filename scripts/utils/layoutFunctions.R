## Define functions for generating row/col positions
layoutRow <- \(y, h, s, n) y + ((h+s) * seq(0,(n-1)))
layoutCol <- \(x, w, s, n) layoutRow(y = x, h = w, s = s, n = n)