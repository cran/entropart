bcShannonBeta <-
function(Ns, Nexp, Correction="Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (bcTsallisBeta(Ns, Nexp, 1, Correction, CheckArguments=FALSE))
}
