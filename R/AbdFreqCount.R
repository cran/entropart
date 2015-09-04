AbdFreqCount <- 
function (Ns, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  
  NsInt <- as.integer(Ns)
  if (any(NsInt != Ns)) stop ("The abundance frequency count requires integer abundances")
  
  # Eliminate 0
  Ns <- Ns[Ns > 0]
  
  DistNs <- tapply(Ns, Ns, length)
  afc <- matrix(c(as.integer(names(DistNs)), DistNs), ncol = 2)
  colnames(afc) <- c("Abundance", "NbSpecies")
  return(afc)
}