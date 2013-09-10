Shannon <-
function(Ps, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (Tsallis(Ps, 1, CheckArguments=FALSE))
}
