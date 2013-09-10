expq <-
function(x, q)
{
  if (q == 1) {
    return (exp(x))
  } else {
    Exponential <- (x*(1-q)+1)^(1/(1-q))
    if (q > 1) {
      Exponential[x > 1/(q-1)] <- NA
    }
    return (Exponential)
  }
}
