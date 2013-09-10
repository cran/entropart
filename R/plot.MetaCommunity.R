plot.MetaCommunity <- 
function (x, ...) {
  barplot(cbind(x$Psi,
                rep(0, x$Nspecies),
                x$Ps              
          ),
          beside = FALSE,
          width = c(x$Wi, .5, 1),
          names.arg = c(names(x$Wi), "", "Metacommunity"),
          ylab = "Species frequencies",
          ...
  )
}