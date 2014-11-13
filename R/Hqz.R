Hqz <- function(Ps, q = 1, Z = diag(length(Ps)), CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()  
  # If names are missing, the probability vector and the similarity vector are assumed to be in the same order
  if (is.null(colnames(Z)) | is.null(names(Ps))) {
    if (ncol(as.matrix(Z)) != length(Ps)) # as.matrix(Z) in case it has been reduced to a single value because of zeros
      # The matrix is square (this has been checked)
      stop("The matrix dimension must equal the probability vector length.")    
    # Eliminate zeros
    Z <- as.matrix(Z)[Ps != 0, Ps != 0]
    Ps <- Ps[Ps != 0]
  } else { # Matrix and Ps are be named.  
    # Eliminate zeros
    Ps <- Ps[Ps != 0]
    if (length(setdiff(names(Ps), colnames(Z))) != 0)
      # The matrix is square (this has been checked)
      stop("Some species are missing in the similarity matrix.")    
    Z <- as.matrix(Z)[names(Ps), names(Ps)]
  }
  
  # Calculate (Zp)
  Zp <- Z %*% Ps
  if (q == 1) {
    # Limit value
    Entropy <- -Ps %*% log(Zp)
  } else {
    # Calculate (Zp)^(q-1)
    Zpqm1 <- Zp^(q-1)
    # Calculate Hqz
    Entropy <- (1-(Ps %*% Zpqm1))/(q-1)
  }
  # Return the value of entropy, as a number rather than a 1x1 matrix
  return (as.numeric(Entropy))
}
