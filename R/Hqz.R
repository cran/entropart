Hqz <- function(Ps, q, Z, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()  
  if (dim(Z)[1] != length(Ps))
    # The matrix is square (this has been checked)
    stop("The matrix dimension must equal the probability vector length.")
  # Matrix and Ps may be named. If the matrix has names, they are identical in rows and columns (this has been checked)
  if (!is.null(colnames(Z)) & !is.null(names(Ps))) {
    # If the matrix and the probability vector are named, names must be identical
    if (setequal(colnames(Z), names(Ps))) {
      # Same names, reorder Ps to fit Z
      Ps <- Ps[colnames(Z)]
    } else {
      stop("Names differ between the matrix and the probability vector.")
    }
  }
    
  # Calculate (Zp) first to test it later
  Zp <- Z %*% Ps
  if (q == 1) {
    # Force 0log0=0 (replace by 0log1)
    Zp[Zp == 0] <- 1
    # Limit value
    Entropy <- -Ps %*% log(Zp)
  } else {
    # Calculate (Zp)^(q-1)
    Zpqm1 <- Zp^(q-1)
    # Force 0^(q-1)=0
    Zpqm1[Zp == 0] <- 0
    # Calculate Hqz
    Entropy <- (1-(Ps %*% Zpqm1)^(q-1))/(q-1)
  }
  # Return the value of diversity, as a number rather than a 1x1 matrix
  return (as.numeric(Entropy))
}
