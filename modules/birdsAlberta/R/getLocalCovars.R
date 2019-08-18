
getLocalCovars <- function(covarTable = sim$covar,
                           disturbanceMap = sim$disturbanceMap,
                           ageMap = sim$ageMap,
                           habitatMap = sim$habitatMap,
                           vegMap = sim$vegMap){
  
invisible({
  
  if(!"ID" %in% names(covarTable)){
    covarTable$ID <- seq(1:nrow(covarTable))
  }
  naCells <- which(is.na(habitatMap[]))
  covarTable <- covarTable[ID %in% naCells, names(covarTable)[1:12] := NA_real_]
  
  # L_CUT
  if(!length(which(disturbanceMap[]==2 | ageMap[]<31))==0){
    L_CUT.rows <- which(disturbanceMap[]==2 | ageMap[]<31)
    covarTable[ID %in% L_CUT.rows, L_CUT := 1L]
  }
  
  # L_ODEC
  if(!length(which(disturbanceMap[]==0 & ageMap[]>90 & habitatMap[]==5))==0){
    L_ODEC.rows <- which(disturbanceMap[]==0 & ageMap[]>90 & habitatMap[]==5) 
    covarTable[ID %in% L_ODEC.rows, L_ODEC := 1L]  
  }
  
  # L_MIX (old)
  if(!length(which(disturbanceMap[]==0 & ageMap[]>90 & habitatMap[]==9))==0){
    L_MIX.rows <- which(disturbanceMap[]==0 & ageMap[]>90 & habitatMap[]==9)
    covarTable[ID %in% L_MIX.rows, L_MIX := 1L]
  }
  
  # L_YDEC
  if(!length(which(disturbanceMap[]==0 & (ageMap[]>30 & ageMap[]<91) & habitatMap[]==4))==0){
    L_YDEC.rows <- which(disturbanceMap[]==0 & (ageMap[]>30 & ageMap[]<91) & habitatMap[]==4)
    covarTable[ID %in% L_YDEC.rows, L_YDEC := 1L]  
  }
  
  # L_MIX (young)
  if(!length(which(disturbanceMap[]==0 & (ageMap[]>30 & ageMap[]<91) & habitatMap[]==9))==0){
    L_MIX.rows <- which(disturbanceMap[]==0 & (ageMap[]>30 & ageMap[]<91) & habitatMap[]==9)
    covarTable[ID %in% L_MIX.rows, L_MIX := 1L]
  }
  
  # L_CC
  closed <- which(vegMap[] %in% c(1:5))
  medium <- which(vegMap[] %in% c(6:7, 12))
  low <- which(vegMap[] %in% c(8:11, 13:15, 25))
  covarTable[ID %in% closed, L_CC := 0.8]
  covarTable[ID %in% medium, L_CC := 0.6]
  covarTable[ID %in% low, L_CC := 0.3]
  
  #L_WDIS: already 0
  
})
  
  return(invisible(covarTable))
}