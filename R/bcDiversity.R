bcDiversity <-
function(Ns, q, Correction="Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  Entropy <- bcTsallis(Ns, q, Correction, CheckArguments=FALSE)
  
  return (expq(Entropy, q))
}
