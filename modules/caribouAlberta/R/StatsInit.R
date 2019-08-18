StatsInit <- function(sim){
  
    sim$maxCaribPop <- (sim$landscapeAttr[[1]]$burnyArea/100)*6
  
  sim$Nt <- P(sim)$N0
  
  if(is.null(sim$Nt)|is.na(sim$Nt)){
    
    sim$Nt <- sim$maxCaribPop
  }
  
  sim$SorensenStats <- data.frame(Lambda = numeric(0), Nt = numeric(0))
  
  return(invisible(sim))
}