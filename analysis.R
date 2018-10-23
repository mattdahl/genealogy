############################################
#	Original script from Clark and Lauderdale's 2012 article "The Genealogy of Law"
# https://doi.org/10.1093/pan/mps019
#
#	Code tidied and modified by Matthew Dahl
############################################

## Set working directory
setwd('/Users/mattdahl/Documents/nd/research/projects/free_expression_doctrine')

## Libs
library(diagram)	# For plotXi function

## Settings
# Labels
caseType <- c('Abortion')
caseTypesLabel <- c('Abortion')
title <- paste(caseTypesLabel, 'Cases')

# Simulation
burn <- 200
sims <- 1000
verbose <- 100

# Bayesian priors
eta <- 0
Founders <- c(1)

#############################################
#	Simulations
#############################################

## Load data, create unique ids for citations
Citation.Data <- readRDS(file = paste('Data/', caseType, '/', caseType, 'CitationData.rds', sep = ''))
Case.Data <- readRDS(file = paste('Data/', caseType, '/', caseType, 'Cases.rds', sep = ''))

## Calculate unique case IDs from the set of cases and precendents
# (remove duplicate IDs that occur when multiple cases cite the same precedents)
UniqueCaseIDs <- sort(unique(c(Citation.Data$CaseID, Citation.Data$PrecedentID)))
UniqueCaseIDs.length <- length(UniqueCaseIDs)

## Create vectors of case names and dates for citing cases
Names <- rep('NA', UniqueCaseIDs.length)
Year <- rep(NA, UniqueCaseIDs.length)
for (i in 1:UniqueCaseIDs.length) {
	Names[i] <- Case.Data$caseName[Case.Data$CaseID == UniqueCaseIDs[i]]
	Year[i] <- Case.Data$date[Case.Data$CaseID == UniqueCaseIDs[i]]
}

## Create lower triangular matrix of citations
Y <- matrix(NA, nrow = UniqueCaseIDs.length, ncol = UniqueCaseIDs.length)
Y[lower.tri(Y)] <- 0
for (k in 1:dim(Citation.Data)[1]) {
  Y[
    match(c(Citation.Data$CaseID[k]), UniqueCaseIDs), # Rows are cases
    match(c(Citation.Data$PrecedentID[k]), UniqueCaseIDs) # Cols are precedents
  ] <- Citation.Data$citations2[k]
}

## Prepare data for the tree generator
T <- dim(Y)[1]
Descendents <- setdiff(1:T, Founders)

## Call the tree generator
source('tree_generator.R')

## Save results
save(Founders, Descendents, Y, Xi.chain, nu.chain, beta.chain, alpha.chain, BestSim, list = c('Xi.chain', 'nu.chain', 'beta.chain', 'alpha.chain'), file = paste('mcmc_samples/', caseType, 'SavedMCMCSample.Data' , sep = ''))
plotXi(Xi.est, file = paste('figures/TreePlots/', caseType, 'TreePlotBest2.pdf' , sep = ''), Names = Names, Years = Year, Title = title)
plotXi(round(Xi.mean), file = paste('figures/TreePlots/', caseType, 'TreePlotMean.pdf', sep=''), Names = Names, Years = Year, Title = title)
