Richness <-
function(Ns, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (sum(Ns > 0))
}
