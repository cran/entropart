bcSimpson <-
function(Ns, Correction="Lande", CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  if (Correction == "Lande" | Correction == "Best") {
    Nrecords <- sum(Ns)
    return(Nrecords/(Nrecords-1)*Tsallis(Ns/Nrecords, 2, CheckArguments=FALSE))  
  } else {
    return(bcTsallis(Ns, 2, Correction, CheckArguments=FALSE)) 
  }
}
