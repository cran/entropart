bcSimpsonBeta <-
function(Ns, Nexp, Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (bcTsallisBeta(Ns, Nexp, 2, Correction, CheckArguments=FALSE))
}
