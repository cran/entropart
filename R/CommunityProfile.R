CommunityProfile <-
function(FUN, NorP, q.seq = seq(0, 2, 0.1), ..., CheckArguments = TRUE) 
{
  if (CheckArguments) {
    CheckentropartArguments()
  }

  # Phylogenetic diversity profile
  return (list(x=q.seq,
               y=sapply(q.seq, function(q) FUN(NorP, q, ..., CheckArguments = FALSE))
              )
         )
  
}