lnq <-
function(x, q)
{
  if (q == 1) {
    return (log(x))
  } else {
    Log <- (x^(1-q)-1)/(1-q)
    Log[x < 0] <- NA
    return (Log)
  }
}
