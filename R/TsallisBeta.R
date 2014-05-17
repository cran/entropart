TsallisBeta <-
function(Ps, Pexp = NULL, q = 1, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (length(Ps) != length(Pexp)) {
    stop("Ps and Pexp should have the same length.")
  }  

  dataBeta <- Ps^q * lnq(Ps/Pexp, q)
  dataBeta[Ps == 0] <- 0

  return (sum(dataBeta))
}
