DivProfile <-
function(q.seq = seq(0, 2, .1), MC, Biased = TRUE, Correction = "Best", Tree = NULL, Normalize = TRUE, Z = NULL, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # Preprocess the tree
  ppTree <- Preprocess.Tree(Tree)
  if (Normalize) {
    Height <- 1
  } else {
    Height <- ppTree$Height
  }  

  # Calculate diversity profile. Parallelize.
  Diversity.seq <- simplify2array(mclapply(q.seq, DivPart, MC = MC, Biased = Biased, Correction = Correction, Tree = ppTree, Normalize = Normalize, Z=Z, CheckArguments =FALSE))

  # Rearrange complex structures
  Dalpha <- unlist(Diversity.seq["CommunityAlphaDiversities", ])
  arrDalpha <- simplify2array(tapply(Dalpha, names(Dalpha), c))
  row.names(arrDalpha) <- q.seq
  Ealpha <- unlist(Diversity.seq["CommunityAlphaEntropies", ])
  arrEalpha <- simplify2array(tapply(Ealpha, names(Ealpha), c))
  row.names(arrEalpha) <- q.seq
  # Prepare a matrix of results
  DivProfile <- list(MetaCommunity = unlist(Diversity.seq["MetaCommunity", 1], use.names=FALSE),
                     Order = unlist(Diversity.seq["Order", ], use.names=FALSE), 
                     Biased = unlist(Diversity.seq["Biased", 1], use.names=FALSE), 
                     Correction = unlist(Diversity.seq["Correction", 1], use.names=FALSE),
                     Normalized = unlist(Diversity.seq["Normalized", 1], use.names=FALSE),
                     CommunityAlphaDiversities = arrDalpha, 
                     CommunityAlphaEntropies = arrEalpha, 
                     TotalAlphaDiversity = unlist(Diversity.seq["TotalAlphaDiversity", ]), 
                     TotalBetaDiversity = unlist(Diversity.seq["TotalBetaDiversity", ]), 
                     GammaDiversity = unlist(Diversity.seq["GammaDiversity", ]), 
                     TotalAlphaEntropy =  unlist(Diversity.seq["TotalAlphaEntropy", ]), 
                     TotalBetaEntropy = unlist(Diversity.seq["TotalBetaEntropy", ]), 
                     GammaEntropy =  unlist(Diversity.seq["GammaEntropy", ])
                    )
  if(!is.null(Tree))
    DivProfile$Tree <- deparse(substitute(Tree)) 
  if(is.null(Z)) {
    DivProfile$Method <- "HCDT"
  } else {
    DivProfile$Method <- "Similarity-based"
    DivProfile$Z <- deparse(substitute(Z))  
  }
  class(DivProfile) <- "DivProfile"
    
  return (DivProfile)
}


is.DivProfile <-
function (x) 
{
  inherits(x, "DivProfile")
}


plot.DivProfile <- 
function (x, ..., main = NULL, xlab = "Order of Diversity", ylab = NULL, Which = "All") 
{
  # Save graphical parameters
  if (Which == "All") {
    op <- graphics::par(no.readonly = TRUE)
    graphics::par(mfrow=c(2, 2))    
  }
  if (Which == "All" | (Which == "Alpha" & is.null(main))) main <- "Total Alpha Diversity"
  if (Which == "All" | (Which == "Alpha" & is.null(ylab))) ylab <- expression(paste(alpha, " diversity"))
  if (Which == "All" | Which == "Alpha") {
    graphics::plot(y=x$TotalAlphaDiversity, x=x$Order, type="l", main=main, xlab=xlab, ylab=ylab, ...)
  }
  if (Which == "All" | (Which == "Communities" & is.null(main))) main <- "Alpha Diversity of Communities"
  if (Which == "All" | (Which == "Communities" & is.null(ylab))) ylab <- expression(paste(alpha, " diversity"))
  if (Which == "All" | Which == "Communities") {
    graphics::plot(x$CommunityAlphaDiversities[, 1] ~ x$Order, type="n", xlim=c(min(x$Order), max(x$Order)), ylim=c(min(x$CommunityAlphaDiversities), max(x$CommunityAlphaDiversities)), main=main, xlab=xlab, ylab=ylab, ...)
    for (Community in (1:ncol(x$CommunityAlphaDiversities))) {
      graphics::lines(x=x$Order, y=x$CommunityAlphaDiversities[, Community], lty=Community)
    }  
    if (Which == "Communities") {
      graphics::legend("topright", colnames(x$CommunityAlphaDiversities), lty=1:ncol(x$CommunityAlphaDiversities), inset=0.01)
    }
  }
  if (Which == "All" | (Which == "Beta" & is.null(main))) main <- "Beta Diversity"
  if (Which == "All" | (Which == "Beta" & is.null(ylab))) ylab <- expression(paste(beta, " diversity"))
  if (Which == "All" | Which == "Beta") {
    graphics::plot(y=x$TotalBetaDiversity, x=x$Order, type="l", main=main, xlab=xlab, ylab=ylab, ...)
  }
  if (Which == "All" | (Which == "Gamma" & is.null(main))) main <- "Gamma Diversity"
  if (Which == "All" | (Which == "Gamma" & is.null(ylab))) ylab <- expression(paste(gamma, " diversity"))
  if (Which == "All" | Which == "Gamma") {
    graphics::plot(y=x$GammaDiversity, x=x$Order, type="l", main=main, xlab=xlab, ylab=ylab, ...)
  }
  # Restore parameters
  if (Which == "All") {
    graphics::par(op)
  }
}


summary.DivProfile <-
function(object, ...) 
{
  cat("Diversity profile of MetaCommunity", object$MetaCommunity, fill=TRUE)
  if (!object$Biased)  
    cat(" with correction:", object$Correction)
  cat("\n")
  
  if (!is.null(object$Tree)) {
    cat("Phylogenetic or functional diversity was calculated according to the tree", object$Tree, "\n", fill=TRUE)
    cat("Diversity is", ifelse(object$Normalized, "normalized", "not normalized"), "\n", fill=TRUE)
  }
  
  cat("Diversity against its order:\n")
  Values <- cbind(object$Order, object$TotalAlphaDiversity, object$TotalBetaDiversity, object$GammaDiversity)
  colnames(Values) <- c("Order", "Alpha Diversity", "Beta Diversity", "Gamma Diversity")
  print(Values)
  
  return(invisible(NULL))
}