############################################
#	Original script from Clark and Lauderdale's 2012 article "The Genealogy of Law"
# https://doi.org/10.1093/pan/mps019
#
#	Code tidied and modified by Matthew Dahl
############################################

## Libs
library(diagram)

#############################################
#	Plotting script
#############################################

plotXi <- function(Xi, file = 'TreePlot.pdf', case_names, case_years, Title) {
  
  getchildYcoords <- function(posMat, parentID, childrenID, parentminY, parentmaxY) {
    minvec <- (cumsum(OffspringCount[childrenID]) - OffspringCount[childrenID]) / sum(OffspringCount[childrenID]) * (parentmaxY - parentminY) + parentminY
    maxvec <- (cumsum(OffspringCount[childrenID])) / sum(OffspringCount[childrenID]) * (parentmaxY - parentminY) + parentminY
    
    print(minvec)
    print(maxvec)
    
    posMat[childrenID, 2] <- (minvec + maxvec) / 2
    
    for (i in 1:length(childrenID)) {
      NextGen <- which((t(Xi) %*% diag(1,T,T))[childrenID[i],] == 1)
      if (length(NextGen) > 0) {
        posMat <- getchildYcoords(posMat,i,NextGen,minvec[i],maxvec[i])
      } 
    }
    
    return(posMat)
  }	
  
  T <- dim(Xi)[1]
  posMat <- matrix(0, T, 2)
  rowCount <- rep(0, T)
  posMat[,1] <- case_years
  
  for (f in Founders) {
    CurrentGeneration <- t(Xi) %*% diag(f,T,T)
    Generation <- rep(0,T)
    OffspringCount <- rep(0,T)
    GenerationCount <- 0
    
    while (sum(CurrentGeneration > 0)) {
      GenerationCount <- GenerationCount + 1
      Generation[which(CurrentGeneration[1,] == 1)] <- GenerationCount
      OffspringCount <- OffspringCount + rowSums(CurrentGeneration)
      CurrentGeneration <- t(Xi) %*% CurrentGeneration
    }
    
    OffspringCount <- OffspringCount + 1
  }
  
  fN <- length(Founders)
  founderminvec <- (1:fN - 1) / fN
  foundermaxvec <- (1:fN) / fN
  foundercentervec <- (founderminvec + foundermaxvec) / 2
  
  for (f in 1:fN) {	
    posMat[Founders[f],2] <- foundercentervec[f]
    NextGen <- which((t(Xi) %*% diag(1,T,T))[Founders[f],] == 1)
    NextGen <- intersect(NextGen, Descendents)
    posMat <- getchildYcoords(posMat, Founders[f], NextGen, founderminvec[f], foundermaxvec[f])	
  }
  
  posMat[,1:2] <- posMat[,2:1] 
  
  pdf(file = file, width = 13, height = 18)
  plot(
    posMat[,1],
    posMat[,2],
    col = 'grey',
    pch = 8,
    ylim = range(c(case_years - 5), max(case_years + 5)),
    xlim = c(-0.2, 1.2),
    ylab = 'Year',
    xlab = '',
    main = Title,
    axes = FALSE
  )
  axis(2)
  
  for (t in 1:T) {
    for (t2 in 1:t) {
      if (Xi[t, t2] == 1) {
        lines(
          c(posMat[t2, 1], posMat[t, 1]),
          c(posMat[t2, 2], posMat[t,2]),
          lwd = 1,
          col = 'grey'
        )
      }
    }
  }
  
  for (t in 1:T) {
    text(posMat[t, 1], posMat[t, 2], case_names[t], cex = 0.5)
  }	
  
  dev.off()
}