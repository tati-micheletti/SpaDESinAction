updateState <- function(sim){ #calcuate merchantable volume by strata
  sim$cutCells <- numeric(0)
  sim$volMap[] <- sim$volMap[]*0
  ytList <- sim$yieldTables[[1]]
  #age the previous cuts
  idx <- which(sim$harvestStateMap[]>0)
  ages <- sim$harvestStateMap[idx] - 1
  sim$harvestStateMap[idx] <- ages
  
  for (i in 1:length(ytList)){
    idx <- which(sim$strataMap[] == i)
    isOK <- which(sim$harvestStateMap[idx] == 0)
    idx <- idx[isOK]
    ageVec <- sim$ageMap[idx]
    isOK <- which(ageVec >= P(sim)$minAge)
    ageVec <- ageVec[isOK]
    idx <- idx[isOK]
    #right truncate ageVec to yieldTableMaxAge
    maxAge <- dim(ytList[[i]])[1]
    ageVec <- ifelse(ageVec > maxAge, maxAge, ageVec)
    volVec <- ytList[[i]][ageVec,3]
    isOK <- which(volVec >= P(sim)$minVolume)
    volVec <- volVec[isOK]
    idx <- idx[isOK]
    sim$volMap[idx] <- volVec
    
    ages <- sim$ageMap[idx]
    if (length(ages) > 0)
      x <- order(ages, decreasing=TRUE)
    else
      x <- integer(0)
    idx <-idx[x]     #put in oldest first
    sim$merchCells[[i]] <- idx
  }
  sim$cutCellIndex <- rep(1,length(ytList))
  return(invisible(sim))  
}