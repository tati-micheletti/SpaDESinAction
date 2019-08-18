
getNeighborhoodCovars <- function(covarTable = sim$covar,
                                  disturbanceMap = sim$disturbanceMap,
                                  ageMap = sim$ageMap,
                                  habitatMap = sim$habitatMap,
                                  vegMap = sim$vegMap){
  
  distMatrix <- matrix(c(rep(1, 12), 0, rep(1,12)), ncol = 5)
  
  cellsN <- data.table::data.table(raster::adjacent(habitatMap, 
                                                    cells = c(1:ncell(habitatMap)), 
                                                    directions = distMatrix, 
                                                    include = FALSE, 
                                                    pairs = TRUE, 
                                                    id = FALSE))
  
  extrHabitat <- habitatMap[cellsN$to]
  extrAge <- ageMap[cellsN$to]
  extrDisturbance <- disturbanceMap[cellsN$to]
  
  cellsN$habitat <- extrHabitat
  cellsN$age <- extrAge
  cellsN$disturbance <- extrDisturbance
  
  countPerNeig <- cellsN[, .N, by = from] # Summary Table
  
  cropMap <- vegMap[]*0
  cropMap[] <- ifelse(cropMap[]==0, 1, NA)
  
  N_CUT <-  cellsN[, sum(disturbance==2 | age<31, na.rm = TRUE), by = from]
  N_CUT <- N_CUT$V1/countPerNeig$N
  N_CUT <- N_CUT*cropMap
  
  N_LATE <-  cellsN[, sum(disturbance==0 & age>90 & habitat==5, na.rm = TRUE), by = from]
  N_LATE <- N_LATE$V1/countPerNeig$N
  N_LATE <- N_LATE*cropMap
  
  N_DEC <-  cellsN[, sum(disturbance==0 & (habitat==4 | habitat==5), na.rm = TRUE), by = from]
  N_DEC <- N_DEC$V1/countPerNeig$N
  N_DEC <- N_DEC*cropMap
  
  N_MIX <-  cellsN[, sum(disturbance==0 & habitat==9, na.rm = TRUE), by = from]
  N_MIX <- N_MIX$V1/countPerNeig$N
  N_MIX <- N_MIX*cropMap
  
  N_SB <-  cellsN[, sum(disturbance==0 & habitat==7, na.rm = TRUE), by = from]
  N_SB <- ifelse(N_SB$V1>0,1,0)
  N_SB <- N_SB*cropMap
  
  N_RICH <-  cellsN[, length(unique(habitat)), by = from]
  N_RICH <- N_RICH$V1
  N_RICH <- N_RICH*cropMap
  
  covarTable$N_CUT <- N_CUT
  covarTable$N_DEC <- N_DEC
  covarTable$N_MIX <- N_MIX
  covarTable$N_LATE <- N_LATE
  covarTable$N_SB <- N_SB
  covarTable$N_RICH <- N_RICH

  return(invisible(covarTable))
}