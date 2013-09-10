ShannonBeta <-
function(Ps, Pexp, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (TsallisBeta(Ps, Pexp, 1, CheckArguments=FALSE))
}
