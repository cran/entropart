bcRao <-
  function(Ns, Tree, Correction="Lande", CheckArguments = TRUE)
  {
    if (CheckArguments)
      CheckentropartArguments()
    
    if (Correction == "Lande"  | Correction == "Best") {
      Nrecords <- sum(Ns)
      return(Nrecords/(Nrecords-1)*AllenH(Ns/sum(Ns), 2, Tree, Normalize=FALSE, CheckArguments=FALSE))  
    } else {
      return(bcPhyloEntropy(Ns, q = 2, Tree = Tree, Normalize = FALSE, Correction = Correction, CheckArguments = FALSE)) 
    }
  }
