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


lnq.CommunityProfile <-
function(Profile)
{
  if (!is.CommunityProfile(Profile))
    stop("Profile must be a CommunityProfile")
  
  CP <- Profile
  CP$y <- sapply(1:length(CP$x), function(i) lnq(CP$y[i], CP$x[i]))
  if (!is.null(CP$low))
    CP$low <- sapply(1:length(CP$x), function(i) lnq(CP$low[i], CP$x[i]))
  if (!is.null(CP$hi))
    CP$hi <- sapply(1:length(CP$x), function(i) lnq(CP$hi[i], CP$x[i]))
  
  return (CP)
}
