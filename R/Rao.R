Rao <-
function(Ps, Tree, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return(AllenH(Ps, 2, Tree, Normalize=FALSE, CheckArguments=FALSE))
}
