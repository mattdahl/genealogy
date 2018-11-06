## UTILITY FUNCTIONS ##

logistic <- function (x) {
  return(1/(1+exp(-x)))
}

logit <- function (x) {
  return(log(x/(1-x)))
}

dpolya <- function (x, p, nu) {
  p <- p / sum(p)
  x <- as.integer(x)
  N <- sum(x)
  return(lgamma(N + 1) - sum(lgamma(x + 1)) + lgamma(sum(nu*p)) - lgamma(N + sum(nu*p)) + sum(lgamma(x + nu*p)) - sum(lgamma(nu*p)))
}

## GENERAL LIKELIHOOD FUNCTION ##
TreeLogLikelihood <- function (alpha, beta, nu, Xi) {
  
  T <- dim(Xi)[1]
  TreeDistance <- matrix(0, T, T)
  
  for (t in Descendents) {
    XipXit <- matrix(0,dim(Xi)[1],dim(Xi)[2])
    XipXit[1:t,1:t] <- Xi[1:t,1:t] + t(Xi)[1:t,1:t]
    
    tMat <- matrix(0,T,1)
    tMat[t,1] <- 1
    
    found <- c(rep(0,(t-1)),rep(1,(T-t+1)))
    foundtwice <- c(rep(0,(t)),rep(1,(T-t)))
    
    d <- 0
    while (sum(foundtwice) < sum(found)) {
      d <- d+1
      tMat <- XipXit %*% tMat
      finds <- (tMat > 0)
      newfinds <- which(finds & (found == 0))
      oldfinds <- which(finds & (found == 1))
      TreeDistance[t,newfinds] <- d
      found[newfinds] <- 1
      foundtwice[oldfinds] <- 1
    }
  }
  
  Unconnected <- TreeDistance == 0
  piMat <- exp(alpha * Unconnected  + beta * TreeDistance) * as.numeric(lower.tri(Unconnected, diag = FALSE))
  piMat <- piMat/rowSums(piMat,na.rm=TRUE)
  piMat[Founders,] <- 0
  LL <- 0
  
  for (t in Descendents) {	
    LL <- LL + dpolya(Y[t,1:(t-1)],piMat[t,1:(t-1)],nu)
  }
  
  return(LL)
}

## PREDICTED PROBABILITIES FUNCTION ##
TreePredictedProbabilities <- function (alpha, beta, Xi) {
  
  T <- dim(Xi)[1]
  XipXit <- Xi + t(Xi)
  
  TreeDistance <- matrix(0,dim(Xi)[1],dim(Xi)[2])
  
  for (t in Descendents) {
    tMat <- matrix(0,T,1)
    tMat[t,1] <- 1
    
    found <- c(rep(0,(t-1)),rep(1,(T-t+1)))
    foundtwice <- c(rep(0,(t)),rep(1,(T-t)))
    
    d <- 0
    while (sum(foundtwice) < sum(found)) {
      d <- d+1
      tMat <- XipXit %*% tMat
      finds <- (tMat > 0)
      newfinds <- which(finds & (found == 0))
      oldfinds <- which(finds & (found == 1))
      TreeDistance[t,newfinds] <- d
      found[newfinds] <- 1
      foundtwice[oldfinds] <- 1
    }
  }
  
  Unconnected <- TreeDistance == 0 
  piMat <- exp(alpha * Unconnected  + beta * TreeDistance) * as.numeric(lower.tri(Unconnected, diag = FALSE))
  piMat <- piMat/rowSums(piMat,na.rm=TRUE)
  piMat[Founders,] <- 0
  
  return(piMat)
}