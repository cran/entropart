Tsallis <-
function(Ps, q, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  dataTsallis <- -Ps^q * lnq(Ps, q)
  # Eliminate unobserved species
  dataTsallis[Ps == 0] <- 0
  
  return (sum(dataTsallis))
}
