CheckentropartArguments <-
function() {

  # Get the list of arguments of the parent function
  ParentFunction <- sys.call(-1)[[1]]
  # If apply() or similar was used, the function name is not in ParentFunction: sys.call(-1)[[1]] returns "FUN"
  if (ParentFunction == "FUN") {
    warning("Function arguments cannot be checked, probably because you used apply(). Add CheckArguments=FALSE to suppress this warning.")
    return (TRUE)
  }
  
  ErrorFunction <- paste("Error in ", ParentFunction, ":")
  Args <- formals(match.fun(ParentFunction))

  ErrorMessage <- function(Message, Argument) {
    cat(paste(ErrorFunction, Message, deparse(substitute(Argument)), "cannot be:\n"))
    print(head(Argument))
    stop("Check the function arguments.", call. = FALSE)
  }

  # Correction 
  if (!is.na(names(Args["Correction"]))) {
    Correction <- eval(expression(Correction), parent.frame())
    if (!is.character(Correction)) {
      ErrorMessage("Correction must be a string.", Correction)
    }
  }
  
  # MC 
  if (!is.na(names(Args["MC"]))) {
    MC <- eval(expression(MC), parent.frame())
    if (!is.MetaCommunity(MC)) {
      ErrorMessage("MC must be a MetaCommunity object.", MC)
    }
  }
  
  # MClist
  if (!is.na(names(Args["MClist"]))) {
    MClist <- eval(expression(MClist), parent.frame())
    if (!is.list(MClist)) {
      ErrorMessage("MClist must be a list.", MClist)
    if (any(!unlist(lapply(MClist, function(x) is.MetaCommunity(x)))))
      ErrorMessage("All elements of MClist must be of class MetaCommunity.", MClist)
    }
  }

  # q 
  if (!is.na(names(Args["q"]))) {
    q <- eval(expression(q), parent.frame())
    if (!is.numeric(q) | length(q)!=1) {
      ErrorMessage("q must be a number.", q)
    }
  }
  # q.seq 
  if (!is.na(names(Args["q.seq"]))) {
    q.seq <- eval(expression(q.seq), parent.frame())
    if (!is.vector(q.seq)) {
      ErrorMessage("q.seq must be a numeric vector.", q.seq)
    }
  }
  
  # Normalize 
  if (!is.na(names(Args["Normalize"]))) {
    Normalize <- eval(expression(Normalize), parent.frame())
    if (!is.logical(Normalize)) {
      ErrorMessage("Normalize must be TRUE or FALSE.", Normalize)  
    }
  }

  # nPoints
  if (!is.na(names(Args["nPoints"]))) {
    nPoints <- eval(expression(nPoints), parent.frame())
    if (!is.numeric(nPoints) | length(nPoints)!=1)
      ErrorMessage("nPoints must be a number.", nPoints)
    if (any(nPoints < 1))
      ErrorMessage("nPoints must be at least 1.", nPoints)
  }
  
  # Ns 
  if (!is.na(names(Args["Ns"]))) {
    Ns <- eval(expression(Ns), parent.frame())
    if (!is.numeric(Ns))
      ErrorMessage("Ns must be numeric.", Ns)
    if (any(Ns < 0))
      ErrorMessage("All abundance values must be positive.", Ns)
  }
  # Nexp 
  if (!is.na(names(Args["Nexp"]))) {
    Nexp <- eval(expression(Nexp), parent.frame())
    if (!is.numeric(Nexp)) {
      ErrorMessage("Nexp must be numeric.", Nexp)
    } 
    if (any(Nexp < 0)){
      ErrorMessage("All abundance values must be positive.", Nexp)
    }
  } 
  
  # Ps 
  if (!is.na(names(Args["Ps"]))) {
    Ps <- eval(expression(Ps), parent.frame())
    if (!is.numeric(Ps))
      ErrorMessage("Ps must be numeric.", Ps)    
    # Probabilities must sum to 1
    if (!isTRUE(all.equal(sum(Ps), 1))) {
      ErrorMessage("Probabilities must sum to 1.", Ps)
    }
    if (any(Ps < 0)) {
      ErrorMessage("All probabilities must be positive.", Ps)
    }
  }
  # Pexp 
  if (!is.na(names(Args["Pexp"]))) {
    Pexp <- eval(expression(Pexp), parent.frame())
    if (!is.numeric(Pexp)) {
      ErrorMessage("Pexp must be numeric.", Pexp) 
    }
    # Probabilities must sum to 1
    if (!isTRUE(all.equal(sum(Pexp), 1))) {
      ErrorMessage("Probabilities must sum to 1.", Pexp) 
    }
    if (any(Pexp < 0)) {
      ErrorMessage("All probabilities must be positive.", Pexp)
    }
  }

  # NorP 
  if (!is.na(names(Args["NorP"]))) {
    NorP <- eval(expression(NorP), parent.frame())
    if (!is.numeric(NorP)) {
      ErrorMessage("NorP must be numeric.", NorP)
    }
    if (any(NorP < 0)) {
      ErrorMessage("All NorP values must be positive.", NorP)   
    }
    # NorP can be a true vector
    if (!is.vector(NorP)) {
      # or a "named vector" whose attributes are not "names". Then dim() return the vector's length.
      if (length(dim(NorP)) != 1) {
        # or a 2d numeric object
        if (length(dim(NorP)) == 2) {
          if (dim(NorP)[2] > 2) {
            # then it must have 1 or 2 columns
            ErrorMessage("NorP may be a vector or a two-column matrix.", NorP)                      
          }
        } else {
          ErrorMessage("NorP may be a vector or a two-column matrix.", NorP)          
        }
      } 
    }
    
    if (length(dim(NorP) != 1)) {
      if ((length(dim(NorP)) != 2) | (dim(NorP)[2] != 2)) {
        ErrorMessage("NorP may be a vector or a two-column matrix.", NorP)
      }  
    }
  }
  
  # Simulations 
  if (!is.na(names(Args["Simulations"]))) {
    Simulations <- eval(expression(Simulations), parent.frame())
    if (!is.numeric(Simulations))
      ErrorMessage("Simulations must be numeric.", Simulations)
    if (!is.vector(Simulations) | length(Simulations) > 1)
      ErrorMessage("Simulations must be a single number.", q.seq)
    if (any(Simulations < 1))
      ErrorMessage("Simulations must be at least 1.", Simulations)
  }

  # Tree 
  if (!is.na(names(Args["Tree"]))) {
    Tree <- eval(expression(Tree), parent.frame())
    if (!is.null(Tree)) {
      if (!inherits(Tree, "phylog") & !inherits(Tree, "hclust") & !inherits(Tree, "PPtree")) {
        ErrorMessage("Tree may be NULL or an object of class hclust or phylog or PPtree.", Tree)
      }
      if (inherits(Tree, "phylog")) {
        if (is.null(Tree$Wdist))
          ErrorMessage("phylog Tree must contain a distance matrix (use add.tools=TRUE when creating it).", Tree)
      }
      
    }  
  }
  # PhyloTree
  if (!is.na(names(Args["PhyloTree"]))) {
    PhyloTree <- eval(expression(PhyloTree), parent.frame())
    if (!is.null(PhyloTree)) {
      if (!inherits(PhyloTree, "phylog") & !inherits(PhyloTree, "hclust") & !inherits(PhyloTree, "PPtree")) {
        ErrorMessage("PhyloTree may be NULL or an object of class hclust or phylog or PPtree.", PhyloTree)
      }  
    }  
  }
  
  # RealValue
  if (!is.na(names(Args["RealValue"]))) {
    RealValue <- eval(expression(RealValue), parent.frame())
    if (!is.numeric(RealValue) | length(RealValue)!=1) {
      ErrorMessage("RealValue must be a number.", RealValue)
    }
  }
  # SimulatedValues
  if (!is.na(names(Args["SimulatedValues"]))) {
    SimulatedValues <- eval(expression(SimulatedValues), parent.frame())
    if (!is.vector(SimulatedValues)) {
      ErrorMessage("SimulatedValues must be a numeric vector.", SimulatedValues)
    }
  }

  # Weights
  if (!is.na(names(Args["Weights"]))) {
    Weights <- eval(expression(Weights), parent.frame())
    if (!is.vector(Weights)) {
      ErrorMessage("Weights must be a numeric vector.", Weights)
    if (any(Weights < 0))
      ErrorMessage("All weights must be positive.", Weights)
    }
  }

  # xy 
  if (!is.na(names(Args["xy"]))) {
    xy <- eval(expression(xy), parent.frame())
    if (!is.numeric(xy) | length(xy)!=2) {
      ErrorMessage("xy must be a numeric vector of length 2.", xy)
    }
  }

  # Z
  if (!is.na(names(Args["Z"]))) {
    Z <- eval(expression(Z), parent.frame())
    if (!is.null(Z)) {
      if (!is.matrix(Z)) {
        ErrorMessage("Z must be a square matrix.", Z)
      } else {
        if (dim(Z)[1] != dim(Z)[2]) {
          ErrorMessage("Z must be a square matrix.", Z)
        }
        if (!is.null(colnames(Z)) | !is.null(rownames(Z))) {
          # If the matrix is named, rows and columns must have the same names
          if (!identical(colnames(Z), rownames(Z)))
            ErrorMessage("Z row and column names must be identical.", Z)
        }
        # Must be a relatedness matrix
        if (any(Z<0))
          ErrorMessage("All terms of the relatedness matrix Z must be positive.", Z)
        if (any(diag(Z)<0))
          ErrorMessage("All terms of the relatedness matrix Z diagonal must be strictly positive.", Z)
      }
    }  
  }
  
  return (TRUE)
}
