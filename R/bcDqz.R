bcDqz <-
function (Ns, q = 1, Z = diag(length(Ns)), Correction = "Best", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  Entropy <- bcHqz(Ns, q, Z, Correction, CheckArguments=FALSE)
  
  return (expq(Entropy, q))
}
