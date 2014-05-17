Diversity <-
function(Ps, q = 1, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  Entropy <- Tsallis(Ps, q, CheckArguments=FALSE)
  
  return (expq(Entropy, q))
}
