plot.PPtree <- 
function (x, ...) {
  
  # Get the hclust object and plot it
  hTree <- x$hTree
  plot(hTree, ...)
  
}
