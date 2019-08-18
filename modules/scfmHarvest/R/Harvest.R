Harvest <- function(sim){
  
  sim <- updateState(sim)
  aac <- unlist(sim$annualCut) #level for region missing
  vol <- rep(0, length(aac))
  maxCells <- unlist(lapply(sim$merchCells,length))
  nStrata <- length(maxCells)
  while (any(vol < aac & sim$cutCellIndex <= maxCells)){
    for (i in 1:nStrata){
      idx <- sim$cutCellIndex[i]
      if (vol[i] < aac[i] && idx <= maxCells[i]){ #need more and there is potentially some
        while (sim$harvestStateMap[sim$merchCells[[i]][idx]] != 0){
          idx = idx + 1
          if (idx > maxCells[i])
            break
        }

        if (idx <= maxCells[i]){ #could we find a cell to cut? 
          cellX <- sim$merchCells[[i]][idx]
          #for a simple block simulator, get adjacent cells
          adjX  <- adjacent(sim$volMap, cellX, directions=8, pairs=FALSE, include=TRUE)
          #select those in same strata, that are merch and not blocked. 
          isOK <- which(sim$strataMap[adjX] == i & sim$volMap[adjX] > 0 & sim$harvestStateMap[adjX] == 0)
          adjX <- adjX[isOK]
          #if the block is big enough, cut them all and mark them
        
          if (length(adjX) * sim$landscapeAttr[[1]]$cellSize >= P(sim)$minBlockSize){
            vol[i] <- vol[i] +  sum(sim$volMap[adjX]) * sim$landscapeAttr[[1]]$cellSize #then get its Volume/ha
            sim$cutCells <- c(sim$cutCells, adjX)
            #should mark all cells adjacent also. 
            adjX2 <- adjacent(sim$volMap, adjX, directions=4, pairs=FALSE, include=TRUE)
            sim$harvestStateMap[adjX2] <- P(sim)$greenUpPeriod #mark it cut
          }
        }
        sim$cutCellIndex[i] <- idx + 1 # this now does not ensure that something is cut every step.
      }
    }
  }
  #sim<-sequentialCut(sim)
  #sim<-sequenceOldest(sim)
  #sim<-blockOldest(sim)
  sim$harvestStats[nrow(sim$harvestStats)+1,] <- vol #area, vol, vol/area; where to find these^
  
  return(invisible(sim))
  
}