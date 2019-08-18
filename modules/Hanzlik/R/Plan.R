Plan <- function(sim) {
  region <- sim$hanzlikPars[[1]]   # we are going to use the first region, if we can
  
  for (i in 1:length(region)){
    
    aac <- 0
    idx <- which(sim$strataMap[] == i)
    
    if (length(idx)>0) {
      hVec <-region[[i]]$hVec  #in future versions there will be list of these,
      #one for each collumn in the YT
      nh <- length(hVec)
      aVec<-rep(0,nh)
      
      x <- tabulate(sim$ageMap[idx])
      
      nx<-length(x)
      if (nx <= nh){
        aVec[1:nx] <- x
      }
      else {
        aVec  <- x[1:nh]
        x[nh] <- x[nh] + sum(x[nh+1:nx]) #accumulate any "missing ones"
      }
      aac<-sum(aVec*hVec) * sim$landscapeAttr[[1]]$cellSize #this needs to be changed.
    }
    sim$annualCut[[i]] <- aac    # 
  }
  return(invisible(sim))
}