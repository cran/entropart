Diversity <-
function(Ps, q, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  Entropy <- Tsallis(Ps, q, CheckArguments=FALSE)
  
  return (expq(Entropy, q))
}
