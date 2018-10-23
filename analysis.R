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
#	SIMULATIONS
#############################################

## Load data, create unique ids for citations
Citation.Data <- readRDS(file = paste('Data/', caseType, '/', caseType, 'CitationData.rds', sep = ''))
Case.Data <- readRDS(file = paste('Data/', caseType, '/', caseType, 'Cases.rds', sep = ''))

## Calculate unique case IDs from the set of cases and precendents
# (remove duplicate IDs that occur when multiple cases cite the same precedents)
UniqueCaseIDs <- sort(unique(c(Citation.Data$CaseID, Citation.Data$PrecedentID)))
UniqueCaseIDs.length <- length(UniqueCaseIDs)

IntExtIDMap <- UniqueCaseIDs
ExtIntIDMap <- rep(NA, max(UniqueCaseIDs))
for (i in 1:length(UniqueCaseIDs)) {
	ExtIntIDMap[UniqueCaseIDs[i]] <- i
}

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
save(Founders,Descendents,Y,Xi.chain,nu.chain,beta.chain,alpha.chain,BestSim,list=c("Xi.chain","nu.chain","beta.chain","alpha.chain"),file=paste("mcmc_samples/",caseType,"SavedMCMCSample.Data",sep=""))
plotXi(Xi.est,file=paste("figures/TreePlots/",caseType,"TreePlotBest2.pdf",sep=""),Names=Names,Years=Year,Title=title)
plotXi(round(Xi.mean),file=paste("figures/TreePlots/",caseType,"TreePlotMean.pdf",sep=""),Names=Names,Years=Year,Title=title)



#############################################
#	ANALYSIS OF SAMPLES
#############################################

## create vectors of variables for figures
allChildren <- NULL
AllAuthority <- NULL
AllCentrality <- NULL
AllYears <- NULL
AllNames <- NULL
AllHub <- NULL
AllOxford <- NULL
AllLii <- NULL
AllNyt <- NULL
##for(i in 1:length(caseTypes)){
	caseType <- caseTypes[1]

	load(file=paste('MCMCSamples/',caseType,'SavedMCMCSample.Data',sep=''))
	Xi <- rowMeans(Xi.chain,dims=2)
	ParentId <- apply(Xi,1,which.max)
	parentMatrix <- matrix(0,nrow=dim(Xi)[1],ncol=dim(Xi)[1])
	for(j in 1:dim(parentMatrix)[1]){
		parentMatrix[j, ParentId[j]] <- 1
	}
	childCases <- colSums(parentMatrix)
	allChildren <- c(allChildren,childCases)

	Citation.Data <- haven::read_dta(file=paste('Data/',caseType,'/',caseType,'CitationData.dta',sep=""))
	Case.Data <- haven::read_dta(file=paste('Data/',caseType,'/',caseType,'Cases.dta',sep=""))


	UniqueCaseIDs <- sort(unique(c(Citation.Data$CaseID,Citation.Data$PrecedentID)))
	IntExtIDMap <- UniqueCaseIDs
	ExtIntIDMap <- rep(NA,max(UniqueCaseIDs))
	for (j in 1:length(UniqueCaseIDs)){
		ExtIntIDMap[UniqueCaseIDs[j]] <- j
		}


	#create vector of dates for citing cases
	Year <- rep(NA,max(ExtIntIDMap, na.rm=T))
	Names <- rep("NA",max(ExtIntIDMap, na.rm=T))
	Authority <- rep(NA,max(ExtIntIDMap, na.rm=T))
	Centrality <- rep(NA,max(ExtIntIDMap, na.rm=T))
	Hub <- rep(NA,max(ExtIntIDMap, na.rm=T))
	Lii <- rep(NA,max(ExtIntIDMap, na.rm=T))
	Oxford <- rep(NA,max(ExtIntIDMap, na.rm=T))
	Nyt <- rep(NA,max(ExtIntIDMap, na.rm=T))
	for(j in 1:max(ExtIntIDMap, na.rm=T)){
	 	Year[j] <- Case.Data$date[Case.Data$CaseID==IntExtIDMap[j]]
	 	Authority[j] <- Case.Data$auth[Case.Data$CaseID==IntExtIDMap[j]]
	 	Centrality[j] <- Case.Data$incent[Case.Data$CaseID==IntExtIDMap[j]]
	 	Hub[j] <- Case.Data$hub[Case.Data$CaseID==IntExtIDMap[j]]
	 	Oxford[j] <- Case.Data$oxford[Case.Data$CaseID==IntExtIDMap[j]]
	 	Lii[j] <- Case.Data$liihc[Case.Data$CaseID==IntExtIDMap[j]]
	 	Nyt[j] <- Case.Data$nyt[Case.Data$CaseID==IntExtIDMap[j]]
		Names[j] <- Case.Data$caseName[Case.Data$CaseID==IntExtIDMap[j]]
	 }
	AllYears <- c(AllYears,Year)
	AllAuthority <- c(AllAuthority,Authority)
	AllCentrality <- c(AllCentrality,Centrality)
	AllOxford <- c(AllOxford,Oxford)
	AllLii <- c(AllLii,Lii)
	AllNyt <- c(AllNyt,Nyt)
	AllHub <- c(AllHub,Hub)	
	AllNames <- c(AllNames,Names)
##}


all.indicators <- cbind(allChildren,AllYears,AllAuthority,AllCentrality,AllOxford,AllLii,AllNyt,AllHub,AllNames)
colnames(all.indicators) <- c('allChildren','AllYears','AllAuthority','AllCentrality','AllOxford','AllLii','AllNyt','AllHub','AllNames')
clean.all.indicators <- all.indicators[is.na(allChildren)==FALSE,]
for(j in 2:9){
	clean.all.indicators <- clean.all.indicators[is.na(clean.all.indicators[,j])==FALSE,]
}



## model fit ##

##	Function to	Create Posterior Predictive Plot ##
TreePredictedProbabilities <- function(alpha,beta,Xi){
	
	T <- dim(Xi)[1]
	XipXit <- Xi + t(Xi)
	
TreeDistance <- matrix(0,dim(Xi)[1],dim(Xi)[2])

for (t in Descendents){
		tMat <- matrix(0,T,1)
		tMat[t,1] <- 1
		
		found <- c(rep(0,(t-1)),rep(1,(T-t+1)))
		foundtwice <- c(rep(0,(t)),rep(1,(T-t)))
		
		d <- 0
		while(sum(foundtwice) < sum(found)){
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


##	figure comparing actual and predicted citation counts	##
pdf('PredictedCitations.pdf',width=15,height=12)
par(mfrow=c(4,5))
##for(i in 1:length(caseTypes)){
	caseType <- caseTypes[1]
  load(file=paste('MCMCSamples/',caseType,'SavedMCMCSample.Data',sep=''))
  
  sims <- dim(Xi.chain)[3]
  T <- dim(Xi.chain)[1]
  Founders <- 1
  Descendents <- setdiff(1:dim(Xi.chain)[1], Founders)
  
  PredictedCitationProbabilities <- matrix(0,T,T)
  
  for (sim in 1:sims){
  PredictedCitationProbabilities <- PredictedCitationProbabilities + TreePredictedProbabilities(alpha.chain[sim],beta.chain[sim], Xi.chain[,,sim])/sims
  }
  
  PredictedVsActual <- data.frame(as.vector(matrix(rowSums(Y,na.rm=TRUE),T,T) * PredictedCitationProbabilities),as.vector(matrix(Y,T,T,dimnames=NULL)))
  PredictedVsActual <- as.matrix(PredictedVsActual)
  colnames(PredictedVsActual) <- NULL
  PredictedVsActual <- PredictedVsActual[rowSums(is.na(PredictedVsActual)) == 0,]
  PredictedVsActual <- PredictedVsActual[sort(PredictedVsActual[,1],index.return=TRUE)$ix,]
  rownames(PredictedVsActual) <- NULL
  colnames(PredictedVsActual) <- NULL
  
  SqrtPredictedVsActual <- sqrt(PredictedVsActual)
  plot(SqrtPredictedVsActual[,1],SqrtPredictedVsActual[,2],xlim=c(0,max(SqrtPredictedVsActual)),ylim=c(0,max(SqrtPredictedVsActual)),xlab="Posterior Predicted Citation Count",ylab="Actual Citation Count",main=caseTypesLabels[i],col="grey",axes=FALSE)
  axis(1,at=0:trunc(max(SqrtPredictedVsActual)),labels=(0:trunc(max(SqrtPredictedVsActual)))^2)
  axis(2,at=0:trunc(max(SqrtPredictedVsActual)),labels=(0:trunc(max(SqrtPredictedVsActual)))^2)
  box()
  abline(0,1)
  loess.out <- loess(PredictedVsActual[,2]~PredictedVsActual[,1])
  lines(sqrt(loess.out$x),sqrt(loess.out$fitted),lwd=2)


##}
dev.off()



## make figure with distribution of parent assignment probabilities
pdf('Drafts/NewParentProbs.pdf',15,12)
par(mfrow=c(4,5))
##for(i in 1:length(caseTypes)){
	caseType <- caseTypes[1]

	load(file=paste('MCMCSamples/',caseType,'SavedMCMCSample.Data',sep=''))
	Xi <- rowMeans(Xi.chain,dims=2)
	XiDescendentCountProbs <- apply(Xi,1,max)
	plot(density(XiDescendentCountProbs),main=caseTypesLabels[i], xlab='Posterior Parent Probability', ylab='Density'
			,xlim=c(0,1))
	abline(v=mean(XiDescendentCountProbs), col='grey')
##}
dev.off()


## plot distribution of total children
pdf('DistributionOfTotalChildren.pdf',6,6)
barplot(table(allChildren), main='Distribution of Total Number\nof Child Cases (Truncated)',
	xlab='Number of Estimated Child Cases',ylab='Density',
	xlim=c(0,14),col='grey')
dev.off()


## figure with distribution of parent age at childbirth
allAfterChildren <- NULL
allBeforeChildren <- NULL
pdf('DistributionOfParentAge.pdf',6,6)
plot(c(0,40),c(0,.15), main='Distribution of Parent\nAge at Childbirth',
	xlab='Age of Estimated Parent Case',ylab='Density',type='n')
colors <- c('black',rep('grey',10))
##for(i in 1:length(caseTypes)){
	caseType <- caseTypes[1]

	load(file=paste('MCMCSamples/',caseType,'SavedMCMCSample.Data',sep=''))
	Xi <- rowMeans(Xi.chain,dims=2)
	ParentId <- apply(Xi,1,which.max)
	
	Citation.Data <- haven::read_dta(file=paste('Data/',caseType,'/',caseType,'CitationData.dta',sep=""))
	Case.Data <- haven::read_dta(file=paste('Data/',caseType,'/',caseType,'Cases.dta',sep=""))

	UniqueCaseIDs <- sort(unique(c(Citation.Data$CaseID,Citation.Data$PrecedentID)))
	IntExtIDMap <- UniqueCaseIDs
	ExtIntIDMap <- rep(NA,max(UniqueCaseIDs))
	for (j in 1:length(UniqueCaseIDs)){
		ExtIntIDMap[UniqueCaseIDs[j]] <- j
		}

	#create vector of dates for citing cases
	Year <- rep(NA,max(ExtIntIDMap, na.rm=T))
	for(j in 1:max(ExtIntIDMap, na.rm=T)){
		Year[j] <- Case.Data$date[Case.Data$CaseID==IntExtIDMap[j]]
	}
	Year <- floor(Year)

	# create matrix telling me which years a parent had a child
	parentMatrix <- matrix(0,nrow=dim(Xi)[1],ncol=dim(Xi)[1])
	parentYears <- matrix(0,nrow=dim(Xi)[1],ncol=dim(Xi)[1])
	ageAtChildbirth <- matrix(0,nrow=dim(Xi)[1],ncol=dim(Xi)[1])
	
	for(j in 1:dim(parentMatrix)[1]){
		parentMatrix[j, ParentId[j]] <- 1
	}
	for(j in 1:dim(parentMatrix)[2]){
		parentYears[,j] <- parentMatrix[,j]*Year
		ageAtChildbirth[,j] <- parentYears[,j] - Year[j]
		
	}
	ageAtChildbirth <- ifelse(ageAtChildbirth<0,NA,ageAtChildbirth)
	lines(density(ageAtChildbirth, na.rm=TRUE), col=colors[i])
	
##}
dev.off()


## compare estimates of children with expert lists and NY Times list
# first, clean up data
children <- as.numeric(clean.all.indicators[,1])
lii <- as.numeric(clean.all.indicators[,6])
ox <- as.numeric(clean.all.indicators[,5])
nyt <- as.numeric(clean.all.indicators[,7])

AllNyt <- AllNyt-1
LiiByChildren <- NULL
OxfordByChildren <- NULL
NytByChildren <- NULL
for(j in 0:15){
	LiiByChildren[j+1] <- mean(AllLii[allChildren==j], na.rm=TRUE)
	OxfordByChildren[j+1] <- mean(AllOxford[allChildren==j], na.rm=TRUE)
	NytByChildren[j+1] <- mean(AllNyt[allChildren==j], na.rm=TRUE)
}

LiiXs <- seq(0:15)[which(!is.na(LiiByChildren))]
OxfordXs <- seq(0:15)[which(!is.na(OxfordByChildren))]
NytXs <- seq(0:15)[which(!is.na(NytByChildren))]

pdf('ComparisonWithLists.pdf',12,4)
par(mfrow=c(1,3))
plot(LiiXs, LiiByChildren[which(!is.na(LiiByChildren))], 
	ylab='Legal Information Institute List Inclusion',xlab='Number of Child Cases',
	main='Legal Information Institute List and Fertility',
	xlim=c(0,15), ylim=c(0,1), type='b')
points(jitter(allChildren),AllLii,pch="|", col='grey')
plot(OxfordXs, OxfordByChildren[which(!is.na(OxfordByChildren))],
	ylab='Oxford Important Case List Inclusion',xlab='Number of Child Cases',
	main='Oxford List and Fertility',
	xlim=c(0,15), ylim=c(0,1), type='b')
points(jitter(allChildren),AllOxford, pch="|", col='grey')
plot(NytXs, NytByChildren[which(!is.na(NytByChildren))],
	ylab='New York Times Coverage',xlab='Number of Child Cases',
	main='Salience and Fertility',
	xlim=c(0,15), ylim=c(0,1), type='b')
points(jitter(allChildren),AllNyt, pch="|", col='grey')
dev.off()


## compare estimates of children with Fowler measures
pdf('Drafts/ComparisonWithFowler.pdf',12,6)
par(mfrow=c(1,2))
plot(jitter(rescale(AllAuthority)),jitter(rescale(allChildren)), pch=21, cex=0.75, col='grey',bg='grey',
	xlab='Normalized Fowler and Jeon Authority Score',ylab='Normalized Number of Child Cases',
	main='Legal Authority and Fertility')
abline(lm(rescale(allChildren)~rescale(AllAuthority)))
plot(jitter(rescale(AllHub)),jitter(rescale(allChildren)), pch=21, cex=0.75, col='grey',bg='grey',
	xlab='Normalized Fowler and Jeon Hub Score',ylab='Normalized Number of Child Cases',
	main='Legal Grounding and Fertility')
abline(lm(rescale(allChildren)~rescale(AllHub)))
dev.off()