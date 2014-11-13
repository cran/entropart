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
    if (ncol(as.matrix(Z)) != length(Ps))  # as.matrix(Z) in case it has been reduced to a single value because of zeros
      # The matrix is square (this has been checked)
      stop("The matrix dimension must equal the probability vector length.")    
    # Eliminate zeros
    Z <- as.matrix(Z)[Ps != 0, Ps != 0]
    Pexp <- Pexp[Ps != 0]
    Ps <- Ps[Ps != 0]
  } else { # Matrix and Ps are be named.
    # Reorder Ps and Pexp
    if (!setequal(names(Ps), names(Pexp)))
      stop("Ps and Pexp should have the names.")
    Pexp <- Pexp[names(Ps)]
    # Eliminate zeros
    Pexp <- Pexp[Ps != 0]
    Ps <- Ps[Ps != 0]
    if (length(setdiff(names(Ps), colnames(Z))) != 0)
      # The matrix is square (this has been checked)
      stop("Some species are missing in the similarity matrix.")    
    Z <- as.matrix(Z)[names(Ps), names(Ps)]
  }
  
  # Calculate (Zp)
  Zps <- Z %*% Ps
  Zpexp <- Z %*% Pexp
  
  dataBeta <- Ps * (lnq(1/Zpexp, q)-lnq(1/Zps, q))

  return (sum(dataBeta))
}
