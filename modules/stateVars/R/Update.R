Update <- function(sim){

  # UPDTAE DISTURBANCE MAP
  
  dtidx <- which(sim$dtMap[] > 0)
  dtval <- sim$dtMap[dtidx]
  dtval <- dtval - 1
  sim$dtMap[dtidx] <- dtval
  unMark <- dtidx[which(dtval == 0)]
  sim$disturbanceMap[unMark] <- 0
  
  
  if (is.data.table(sim$spreadState) && nrow(sim$spreadState) > 0){   #existant and non-empty?
    if (!is.null(sim$rstCurrentBurn)){
      valsCurrBurn <- data.table(pixelID = 1:length(sim$rstCurrentBurn), vals = getValues(sim$rstCurrentBurn))
      idx <- valsCurrBurn[vals == 1, pixelID]
    } else {
      idx <- sim$spreadState[state == "activeSource", pixels]         #scfmSpread will define "indices" for active burning pixels if the rstCurrBuen doesn't exist yet
    }
    raster::values(sim$disturbanceMap)[idx] <- 1
    raster::values(sim$dtMap)[idx] <- P(sim)$persistTimes[1]
  }      
  
  if (is.numeric(sim$cutCells) && length(sim$cutCells) > 0){
    adjx <- raster::adjacent(sim$disturbanceMap, sim$cutCells, pairs=FALSE)
    raster::values(sim$disturbanceMap)[adjx] <- 3
    raster::values(sim$dtMap)[adjx] <- P(sim)$persistTimes[3] 
    raster::values(sim$disturbanceMap)[sim$cutCells] <- 2 #do these after, because overlap when cutting blocks of cells.
    raster::values(sim$dtMap)[sim$cutCells] <- P(sim)$persistTimes[2]
    
    # UPDTAE AGE MAP
    raster::values(sim$ageMap)[sim$cutCells] <- 0
    
    # UPDTAE HABITAT MAP
    
        # Cells cut from deciduous become 12 in habitat 
    DEC <- which(sim$habitatMap[]==4|sim$habitatMap[]==5)
    cutDEC <- sim$cutCells[sim$cutCells %in% DEC]
    sim$habitatMap[cutDEC] <- 12
       
        # Cells between 31 -  90 AND class 12 become 4 
    yDEC <- which(sim$habitatMap[]==12&(sim$ageMap[]>30&sim$ageMap[]<91))
    sim$habitatMap[yDEC] <- 4
    
        # Cells older than 90 AND class 4 become class 5
    oDEC <- which(sim$habitatMap[]==4&sim$ageMap[]>90)
    sim$habitatMap[oDEC] <- 5

  }
  
  # ADD TO OUTPUTS IF I CREATE NEW STUFF IN SIM
  
  sim <- heightFromAge(sim)
  
  return(invisible(sim))
}
