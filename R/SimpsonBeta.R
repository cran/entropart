SimpsonBeta <-
function(Ps, Pexp, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (TsallisBeta(Ps, Pexp, 2, CheckArguments=FALSE))
}
