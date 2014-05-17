HqzBeta <-
function(Ps, Pexp = NULL, q = 1, Z = diag(length(Ps)), CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (length(Ps) != length(Pexp)) {
    stop("Ps and Pexp should have the same length.")
  }  
  # If names are missing, the probability vector and the similarity vector are assumed to be in the same order
  if (is.null(colnames(Z)) | is.null(names(Ps))) {
    if (ncol(Z) != length(Ps))
      # The matrix is square (this has been checked)
      stop("The matrix dimension must equal the probability vector length.")    
    # Eliminate zeros
    Z <- Z[Pexp != 0, Pexp != 0]
    Ps <- Ps[Pexp != 0]
    Pexp <- Pexp[Pexp != 0]
  } else { # Matrix and Ps are be named.
    # Reorder Ps and Pexp
    if (!setequal(names(Ps), names(Pexp)))
      stop("Ps and Pexp should have the names.")
    Pexp <- Pexp[names(Ps)]
    # Eliminate zeros
    Ps <- Ps[Pexp != 0]
    Pexp <- Pexp[Pexp != 0]
    if (length(setdiff(names(Ps), colnames(Z))) != 0)
      # The matrix is square (this has been checked)
      stop("Some species are missing in the similarity matrix.")    
    Z <- Z[names(Ps), names(Ps)]
  }
  
  # Calculate (Zp)
  Zps <- Z %*% Ps
  Zpexp <- Z %*% Pexp
  
  dataBeta <- Ps * (lnq(1/Zpexp, q)-lnq(1/Zps, q))

  return (sum(dataBeta))
}
