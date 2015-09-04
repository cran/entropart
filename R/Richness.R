Richness <-
function(NorP, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE, CheckArguments = TRUE, Ps = NULL, Ns = NULL) 
{
  UseMethod("Richness")
}


Richness.ProbaVector <-
function(NorP, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE, CheckArguments = TRUE, Ps = NULL, Ns = NULL)  
{
  if (CheckArguments)
    CheckentropartArguments()
  
  return (sum(NorP > 0))
}


Richness.AbdVector <-
function(NorP, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE, CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  return(bcRichness(Ns=NorP, Correction=Correction, Alpha=Alpha, JackOver=JackOver, CheckArguments=CheckArguments))
}


Richness.integer <-
function(NorP, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE, CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ns)) {
      NorP <- Ns
    } else {
      stop("An argument NorP or Ns must be provided.")
    }
  }
  return(bcRichness(Ns=NorP, Correction=Correction, Alpha=Alpha, JackOver=JackOver, CheckArguments=CheckArguments))
}


Richness.numeric <-
function(NorP, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE, CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  if (missing(NorP)){
    if (!missing(Ps)) {
      NorP <- Ps
    } else {
      if (!missing(Ns)) {
        NorP <- Ns
      } else {
        stop("An argument NorP or Ps or Ns must be provided.")
      }
    }
  }
  
  if (abs(sum(NorP) - 1) < 3*.Machine$double.eps) {
    # Probabilities sum to 1, allowing rounding error
    return(Richness.ProbaVector(NorP, CheckArguments=CheckArguments))
  } else {
    # Abundances
    return(Richness.AbdVector(NorP, Correction=Correction, Alpha=Alpha, JackOver=JackOver, CheckArguments=CheckArguments))
  }
}


bcRichness <- 
function(Ns, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE, CheckArguments = TRUE)
{
  if (CheckArguments)
    CheckentropartArguments()
  
  # No correction
  if (Correction == "None") {
    return(Richness.ProbaVector(Ns/sum(Ns), CheckArguments = FALSE))
  }
  
  # Eliminate 0 and calculate basic statistics
  Ns <- Ns[Ns > 0]
  N <- sum(Ns)
  S <- length(Ns)
  AFC <- AbdFreqCount(Ns)
  S1 <- AFC[AFC[, 1] == 1][2]
  S2 <- AFC[AFC[, 1] == 2][2]
  
  # No correction
  if (Correction == "None") {
    return(S)
  }
  
  # Chao1
  if ((Correction == "Chao1") | (Correction == "iChao1")) {
    if (is.na(S1)) {
      S1 <- 0
    }
    if (is.na(S2)) {
      S0 <- (N-1)/N*S1*(S1-1)/2
    } else {
      S0 <- (N-1)/N*S1*S1/2/S2
    }
  }
  if (Correction == "Chao1") {
    return(S + S0)
  }
  if (Correction == "iChao1") {
    S3 <- AFC[AFC[, 1] == 3][2]
    S4 <- AFC[AFC[, 1] == 4][2]
    if (is.na(S3)) {
      S3 <- 0
    }
    if (is.na(S4)) {
      S4 <- 1
    }
    iS0 <- S3/4/S4*max(S1-S2*S3/2/S4, 0)
    return(S + S0 + iS0)
  }
  
  if (Correction == "Jackknife") {
    # Adapted from jackknife in SPECIES
    # Max possible order
    k <- min(nrow(AFC) - 1, 10)
    # Max number of individuals
    m <- max(AFC[, 1])
    # Complete the abundance frequency count for all counts between 1 and m
    ntemp <- cbind(c(1:m), rep(0, m))
    ntemp[AFC[, 1], 2] <- AFC[, 2]
    AFC <- ntemp
    # Prepare a matrix with k+1 rows and 5 columns
    gene <- matrix(0, nrow=k + 1, ncol=5)
    gene[1, 1] <- S
    for (i in 1:k) {
      gene[i + 1, 1] <- S
      gene[i + 1, 4] <- S
      for (j in 1:i) {
        gene[i + 1, 1] <- gene[i + 1, 1] + (-1)^(j + 1) * 2^i * stats::dbinom(j, i, 0.5) * AFC[j, 2]
        gene[i + 1, 4] <- gene[i + 1, 4] + (-1)^(j + 1) * 2^i * stats::dbinom(j, i, 0.5) * AFC[j, 2] * prod(1:j)
      }
      gene[i + 1, 2] <- -gene[i + 1, 1]
      for (j in 1:i) {
        gene[i + 1, 2] <- gene[i + 1, 2] + ((-1)^(j + 1) * 2^i * stats::dbinom(j, i, 0.5) + 1)^2 * AFC[j, 2]
      }
      gene[i + 1, 2] <- gene[i + 1, 2] + sum(AFC[(i + 1):nrow(AFC), 2])
      gene[i + 1, 2] <- sqrt(gene[i + 1, 2])
    }
    if (k > 1) {
      for (i in 2:k) {
        gene[i, 3] <- -(gene[i + 1, 1] - gene[i, 1])^2/(S - 1)
        for (j in 1:(i - 1)) {
          gene[i, 3] <- gene[i, 3] + ((-1)^(j + 1) * 2^(i) * stats::dbinom(j, i, 0.5) - (-1)^(j + 1) * 2^(i - 1) * stats::dbinom(j, i - 1, 0.5))^2 * AFC[j, 2] * S/(S - 1)
        }
        gene[i, 3] <- gene[i, 3] + AFC[i, 2] * S/(S - 1)
        gene[i, 3] <- sqrt(gene[i, 3])
        gene[i, 5] <- (gene[i + 1, 1] - gene[i, 1])/gene[i, 3]
      }
    }
    # Threshold for Burnham and Overton's test
    coe <- qnorm(1 - Alpha/2, 0, 1)
    # Which orders pass the test?
    x <- (gene[2:(k + 1), 5] < coe)
    if (sum(x, na.rm=TRUE) == 0) {
      # If none, keep the max value of k
      jackest <- gene[k + 1, 1]
      # Estimated standard error
      sej <- gene[k + 1, 2]
    }
    else {
      # Else, keep the smallest value (+1 because jack1 is in line 2)
      Smallestk <- which(x)[1] + 1
      if (JackOver & (Smallestk < k+1))
        # Add 1 if overestimation is requested
        Smallestk <- Smallestk +1
      # Estimated value
      jackest <- gene[Smallestk, 1]
      # Estimated standard error
      sej <- gene[Smallestk, 2]
    }
    return(jackest) 
  }
}
