############################################
#	Original script from Clark and Lauderdale's 2012 article "The Genealogy of Law"
# https://doi.org/10.1093/pan/mps019
#
#	Code tidied and modified by Matthew Dahl
############################################

#############################################
#	Tree generator script
#############################################


## SAMPLE POSTERIOR BY MCMC
# Start Values

Xi <- matrix(0, T, T)
for (t in Descendents) {
	Xi[t, 1] <- 1
}
Xi[2, 1] <- 1

XiMHProposal <- (Y + 1) / rowSums(Y + 1, na.rm = TRUE)
colnames(XiMHProposal) <- NULL
rownames(XiMHProposal) <- NULL

alpha <- -5
beta <- -1
nu <- 10

CurrentLL <- TreeLogLikelihood(alpha, beta, nu, Xi)
LL.chain <- rep(NA, sims)
alpha.chain <- rep(NA, sims)
beta.chain <- rep(NA, sims)
nu.chain <- rep(NA, sims)
Xi.chain <- array(NA, c(T, T, sims))

accept.alpha.chain <- rep(NA, sims)
accept.beta.chain <- rep(NA, sims)
accept.nu.chain <- rep(NA, sims)
accept.Xi.chain <- matrix(NA, T, sims)

accept.alpha <- 0
accept.beta <- 0
accept.nu <- 0
accept.Xi <- rep(0, T)

print(paste('MCMC started at:', date()))


# Sample Tree Connections by Gibbs #

for (sim in -(burn-1):sims) {
  for (t in Descendents) {
  	currentparent <- which(Xi[t,1:(t-1)] == 1)
  	Xi.proposed <- Xi
  	Xi.proposed[t,] <- 0
  	
  	currentRow <- Xi[t,]
  		
  	proposedparent <- which(rowSums(rmultinom(1,1,XiMHProposal[t,1:(t-1)])) > 0)
  	Xi.proposed[t, proposedparent] <- 1
  	proposedRow <- Xi.proposed[t,]
  		
  	ProposalLL <- TreeLogLikelihood(alpha,beta,nu,Xi.proposed)
  	likelihoodRatio <- exp(ProposalLL-CurrentLL)
  		
  	proposalRatio <- dmultinom(currentRow[1:(t-1)],1,XiMHProposal[t,1:(t-1)])/dmultinom(proposedRow[1:(t-1)],1,XiMHProposal[t,1:(t-1)])
  	
  	acceptRatio <- likelihoodRatio*proposalRatio
  	accept.Xi[t] <- rbinom(1,1,min(c(acceptRatio,1)))
  	Xi <- Xi*(1-accept.Xi[t]) + Xi.proposed*(accept.Xi[t])
  	CurrentLL <- CurrentLL*(1-accept.Xi[t]) + ProposalLL*(accept.Xi[t])
  }
  
  	
  # Sample Alpha by Independence Metropolis With Uniform Proposal#
  
  alpha.proposed <- runif(1,-10,-5)
  ProposalLL <- TreeLogLikelihood(alpha.proposed,beta,nu,Xi)
  acceptRatio <- exp(ProposalLL-CurrentLL)
  accept.alpha <- rbinom(1,1,min(c(acceptRatio,1)))
  alpha <- alpha*(1-accept.alpha) + alpha.proposed*(accept.alpha)
  CurrentLL <- CurrentLL*(1-accept.alpha) + ProposalLL*(accept.alpha)
  
  # Sample Beta by Independence Metropolis With Uniform Proposal#
  
  beta.proposed <- runif(1,-2,0)
  ProposalLL <- TreeLogLikelihood(alpha,beta.proposed,nu,Xi)
  acceptRatio <- exp(ProposalLL-CurrentLL)
  accept.beta <- rbinom(1,1,min(c(acceptRatio,1)))
  beta <- beta*(1-accept.beta) + beta.proposed*(accept.beta)
  CurrentLL <- CurrentLL*(1-accept.beta) + ProposalLL*(accept.beta)
  
  # Sample Gamma by Independence Metropolis With Uniform Proposal#
  
  nu.proposed <- runif(1,0,50)
  ProposalLL <- TreeLogLikelihood(alpha,beta,nu.proposed,Xi)
  acceptRatio <- exp(ProposalLL-CurrentLL)
  accept.nu <- rbinom(1,1,min(c(acceptRatio,1)))
  nu <- nu*(1-accept.nu) + nu.proposed*(accept.nu)
  CurrentLL <- CurrentLL*(1-accept.nu) + ProposalLL*(accept.nu)
  
  ## Store Draws to Chain ##
  
  if (sim > 0) {
    LL.chain[sim] <- CurrentLL
    
    alpha.chain[sim] <- alpha
    beta.chain[sim] <- beta
    nu.chain[sim] <- nu
    Xi.chain[,,sim] <- Xi
    
    accept.alpha.chain[sim] <- accept.alpha
    accept.beta.chain[sim] <- accept.beta
    accept.nu.chain[sim] <- accept.nu
    accept.Xi.chain[,sim] <- accept.Xi
  }
  
  if ((sim %% verbose) == 0) {
    print(paste(sim, 'iterations completed at:', date()))
  }
}

print(paste('MCMC finished at:', date()))

# Compute Posterior Statistics #

Xi.average <- rowMeans(Xi.chain, dims = 2)
alpha.est <- mean(alpha.chain)
alpha.err <- sd(alpha.chain)
beta.est <- mean(beta.chain)
beta.err <- sd(beta.chain)
nu.est <- mean(nu.chain)
nu.err <- sd(nu.chain)
LL.est <- mean(LL.chain)
LL.err <- sd(LL.chain)

BestSim <- which.max(LL.chain)

# Print Results #

print(paste("Alpha Acceptance Ratio:",mean(accept.alpha)))
print(paste("Beta Acceptance Ratio:",mean(accept.beta)))
print(paste("Nu Acceptance Ratio:",mean(accept.nu)))
print(paste("Mean Xi Acceptance Ratio:",round(mean(accept.Xi,na.rm=TRUE),3)))

print(paste("Best Posterior Xi (LL: ",round(max(LL.chain),2),"):",sep=""))
Xi.est <- Xi.chain[,,BestSim]
colnames(Xi.est) <- rownames(Xi.est) <- case_names
print(Xi.est)
print("Average Posterior Xi:")
Xi.mean <- rowMeans(Xi.chain,dims=2)
colnames(Xi.mean) <- rownames(Xi.mean) <- case_names
print(round(Xi.mean,2))

print(paste("Posterior Alpha: ",round(alpha.est,2)," (",round(alpha.err,2),")",sep=""))
print(paste("Posterior Beta: ",round(beta.est,2)," (",round(beta.err,2),")",sep=""))
print(paste("Posterior Nu: ",round(nu.est,2)," (",round(nu.err,2),")",sep=""))
print(paste("Posterior Deviance: ",round(LL.est,2)," (",round(LL.err,2),")",sep=""))


Ypredicted <- round(TreePredictedProbabilities(alpha,beta,Xi.est) * rowSums(Y, na.rm = TRUE), 1)
Yresiduals <- Y - Ypredicted

print("Data:")
colnames(Y) <- rownames(Y) <- case_names
print(Y)
print("Model Predicted Data:")
colnames(Ypredicted) <- rownames(Ypredicted) <- case_names
print(Ypredicted)
print("Residuals:")
colnames(Yresiduals) <- rownames(Yresiduals) <- case_names
print(Yresiduals)

