############################################
#	Original script from Clark and Lauderdale's 2012 article "The Genealogy of Law"
# https://doi.org/10.1093/pan/mps019
#
#	Code tidied and modified by Matthew Dahl
############################################

## Libs
library(scales)

#############################################
#	Robustness checking script
#############################################

## Define case types
caseTypes <- c('Abortion')
caseTypesLabels <- c('Abortion')
caseType <- 'Abortion'

## Create vectors of variables for figures
allChildren <- NULL
AllAuthority <- NULL
AllCentrality <- NULL
AllYears <- NULL
AllNames <- NULL
AllHub <- NULL
AllOxford <- NULL
AllLii <- NULL
AllNyt <- NULL

## Load case data
citation_data <- readRDS(file = paste('Data/', caseType, '/', caseType, 'CitationData.rds', sep = ''))
case_data <- readRDS(file = paste('Data/', caseType, '/', caseType, 'Cases.rds', sep = ''))

## Load saved MCMC data
load(file = paste('mcmc_samples/', caseType, 'SavedMCMCSample.Data', sep = ''))

## Calculate unique case IDs from the set of cases and precendents
# (remove duplicate IDs that occur when multiple cases cite the same precedents)
case_ids <- sort(unique(c(citation_data$CaseID, citation_data$PrecedentID)))
case_ids.length <- length(case_ids)

## Create vectors of various variables for citing cases
Year <- rep(NA, max(case_ids.length, na.rm = T))
Names <- rep('NA', max(case_ids.length, na.rm = T))
Authority <- rep(NA, max(case_ids.length, na.rm = T))
Centrality <- rep(NA, max(case_ids.length, na.rm = T))
Hub <- rep(NA, max(case_ids.length, na.rm = T))
Lii <- rep(NA, max(case_ids.length, na.rm = T))
Oxford <- rep(NA, max(case_ids.length, na.rm = T))
Nyt <- rep(NA, max(case_ids.length, na.rm = T))

for (i in 1:case_ids.length) {
  Year[i] <- case_data$date[case_data$CaseID == case_ids[i]]
  Authority[i] <- case_data$auth[case_data$CaseID == case_ids[i]]
  Centrality[i] <- case_data$incent[case_data$CaseID == case_ids[i]]
  Hub[i] <- case_data$hub[case_data$CaseID == case_ids[i]]
  Oxford[i] <- case_data$oxford[case_data$CaseID == case_ids[i]]
  Lii[i] <- case_data$liihc[case_data$CaseID == case_ids[i]]
  Nyt[i] <- case_data$nyt[case_data$CaseID == case_ids[i]]
  Names[i] <- case_data$caseName[case_data$CaseID == case_ids[i]]
}

AllYears <- c(AllYears,Year)
AllAuthority <- c(AllAuthority,Authority)
AllCentrality <- c(AllCentrality,Centrality)
AllOxford <- c(AllOxford,Oxford)
AllLii <- c(AllLii,Lii)
AllNyt <- c(AllNyt,Nyt)
AllHub <- c(AllHub,Hub)	
AllNames <- c(AllNames,Names)


all.indicators <- cbind(allChildren,AllYears,AllAuthority,AllCentrality,AllOxford,AllLii,AllNyt,AllHub,AllNames)
colnames(all.indicators) <- c('allChildren','AllYears','AllAuthority','AllCentrality','AllOxford','AllLii','AllNyt','AllHub','AllNames')
clean.all.indicators <- all.indicators[is.na(allChildren)==FALSE,]
for (j in 2:9) {
  clean.all.indicators <- clean.all.indicators[is.na(clean.all.indicators[,j])==FALSE,]
}

Xi <- rowMeans(Xi.chain, dims = 2)
ParentId <- apply(Xi, 1, which.max)
parentMatrix <- matrix(0, nrow = dim(Xi)[1], ncol = dim(Xi)[1])
for (i in 1:dim(parentMatrix)[1]) {
  parentMatrix[i, ParentId[i]] <- 1
}
childCases <- colSums(parentMatrix)
allChildren <- c(allChildren, childCases)


## Figure comparing actual and predicted citation counts
(function () {
  
  pdf('figures/Robustness/PredictedCitations.pdf',width=15,height=12)
  par(mfrow=c(4,5))

  sims <- dim(Xi.chain)[3]
  T <- dim(Xi.chain)[1]
  Founders <- 1
  Descendents <- setdiff(1:dim(Xi.chain)[1], Founders)
  
  PredictedCitationProbabilities <- matrix(0,T,T)
  
  for (sim in 1:sims) {
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
})()


## Figure with distribution of parent assignment probabilities
(function () {
  pdf('figures/Robustness/NewParentProbs.pdf',15,12)
  par(mfrow=c(4,5))

    Xi <- rowMeans(Xi.chain,dims=2)
  XiDescendentCountProbs <- apply(Xi,1,max)
  plot(density(XiDescendentCountProbs),main=caseTypesLabels[i], xlab='Posterior Parent Probability', ylab='Density'
       ,xlim=c(0,1))
  abline(v=mean(XiDescendentCountProbs), col='grey')
  dev.off()
})()

## Figure with distribution of total children
(function () {
  pdf('figures/Robustness/DistributionOfTotalChildren.pdf',6,6)
  barplot(table(allChildren), main='Distribution of Total Number\nof Child Cases (Truncated)',
          xlab='Number of Estimated Child Cases',ylab='Density',
          xlim=c(0,14),col='grey')
  dev.off()
})()

## Figure with distribution of parent age at childbirth
(function () {
  allAfterChildren <- NULL
  allBeforeChildren <- NULL
  pdf('figures/Robustness/DistributionOfParentAge.pdf',6,6)
  plot(c(0,40),c(0,.15), main='Distribution of Parent\nAge at Childbirth',
       xlab='Age of Estimated Parent Case',ylab='Density',type='n')
  colors <- c('black',rep('grey',10))

  Xi <- rowMeans(Xi.chain,dims=2)
  ParentId <- apply(Xi,1,which.max)

  # create matrix telling me which years a parent had a child
  parentMatrix <- matrix(0,nrow=dim(Xi)[1],ncol=dim(Xi)[1])
  parentYears <- matrix(0,nrow=dim(Xi)[1],ncol=dim(Xi)[1])
  ageAtChildbirth <- matrix(0,nrow=dim(Xi)[1],ncol=dim(Xi)[1])
  
  for (j in 1:dim(parentMatrix)[1]) {
    parentMatrix[j, ParentId[j]] <- 1
  }
  for (j in 1:dim(parentMatrix)[2]) {
    parentYears[,j] <- parentMatrix[,j]*Year
    ageAtChildbirth[,j] <- parentYears[,j] - Year[j]
    
  }
  ageAtChildbirth <- ifelse(ageAtChildbirth<0,NA,ageAtChildbirth)
  lines(density(ageAtChildbirth, na.rm=TRUE), col=colors[i])
  
  dev.off()
})()


## Figure that compares estimates of children with expert lists and NY Times list
(function () {
  # first, clean up data
  children <- as.numeric(clean.all.indicators[,1])
  lii <- as.numeric(clean.all.indicators[,6])
  ox <- as.numeric(clean.all.indicators[,5])
  nyt <- as.numeric(clean.all.indicators[,7])
  
  AllNyt <- AllNyt-1
  LiiByChildren <- NULL
  OxfordByChildren <- NULL
  NytByChildren <- NULL
  for (j in 0:15) {
    LiiByChildren[j+1] <- mean(AllLii[allChildren==j], na.rm=TRUE)
    OxfordByChildren[j+1] <- mean(AllOxford[allChildren==j], na.rm=TRUE)
    NytByChildren[j+1] <- mean(AllNyt[allChildren==j], na.rm=TRUE)
  }
  
  LiiXs <- seq(0:15)[which(!is.na(LiiByChildren))]
  OxfordXs <- seq(0:15)[which(!is.na(OxfordByChildren))]
  NytXs <- seq(0:15)[which(!is.na(NytByChildren))]
  
  pdf('figures/Robustness/ComparisonWithLists.pdf',12,4)
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
})()


## Figures that compares estimates of children with Fowler measures
(function () {
  pdf('figures/Robustness/ComparisonWithFowler.pdf', 12, 6)
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
})()