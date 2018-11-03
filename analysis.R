############################################
#	Original script from Clark and Lauderdale's 2012 article "The Genealogy of Law"
# https://doi.org/10.1093/pan/mps019
#
#	Code tidied and modified by Matthew Dahl
############################################

## Set working directory
setwd('/Users/mattdahl/Documents/nd/research/projects/free_expression_doctrine')

## Settings
# Labels
caseType <- c('Abortion')
caseTypesLabel <- caseType
title <- paste(caseTypesLabel, 'Cases')

# Simulation
burn <- 200
sims <- 1000
verbose <- 100

# Number of parents per tree
Founders <- c(1)

#############################################
#	Analysis script
#############################################

## Load data
citation_data <- readRDS(file = paste('Data/', caseType, '/', caseType, 'CitationData.rds', sep = ''))
case_data <- readRDS(file = paste('Data/', caseType, '/', caseType, 'Cases.rds', sep = ''))

## Calculate unique case IDs from the set of cases and precendents
# (remove duplicate IDs that occur when multiple cases cite the same precedents)
case_ids <- sort(unique(c(citation_data$CaseID, citation_data$PrecedentID)))
case_ids.length <- length(case_ids)

## Create vectors of case names and dates for citing cases
case_names <- rep('NA', case_ids.length)
case_years <- rep(NA, case_ids.length)
for (i in 1:case_ids.length) {
  case_names[i] <- case_data$caseName[case_data$CaseID == case_ids[i]]
  case_years[i] <- case_data$date[case_data$CaseID == case_ids[i]]
}

## Create lower triangular matrix of citations
Y <- matrix(NA, nrow = case_ids.length, ncol = case_ids.length)
Y[lower.tri(Y)] <- 0
for (i in 1:dim(citation_data)[1]) {
  # Put the number of case-precedent citations in the Y matrix at position [case, precedent]
  Y[
    match(c(citation_data$CaseID[i]), case_ids), # Rows are cases
    match(c(citation_data$PrecedentID[i]), case_ids) # Cols are precedents
  ] <- citation_data$citations2[i]
}

## Prepare data for the tree generator
T <- dim(Y)[1]
Descendents <- setdiff(1:T, Founders)

## Estimate the tree
source('tree_generator.R')

## Save the resulting MCMC estimates
save(Founders, Descendents, Y, Xi.chain, nu.chain, beta.chain, alpha.chain, BestSim, list = c('Xi.chain', 'nu.chain', 'beta.chain', 'alpha.chain'), file = paste('mcmc_samples/', caseType, 'SavedMCMCSample.Data' , sep = ''))

## Plot the results
source('plotter.R')
plotXi(Xi.est, file = paste('figures/TreePlots/', caseType, 'TreePlotBest2.pdf' , sep = ''), case_names = case_names, case_years = case_years, Title = title)
plotXi(round(Xi.mean), file = paste('figures/TreePlots/', caseType, 'TreePlotMean.pdf', sep=''), case_names = case_names, case_years = case_years, Title = title)
