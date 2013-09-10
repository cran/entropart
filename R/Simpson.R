Simpson <-
function(Ps, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (Tsallis(Ps, 2, CheckArguments=FALSE))
}
