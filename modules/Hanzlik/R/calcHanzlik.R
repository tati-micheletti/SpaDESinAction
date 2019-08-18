calcHanzlik <- function(ytList,rpm){
  cHelper <- function(yt,rpm){
    #browser()
    if (!is.matrix(yt))
      stop("invalid yield table")
    nr <- dim(yt)[1]
    nc <- dim(yt)[2] #we will assume last collumn has the total volume, and manage on that for now
    vt<- yt[,nc]/1:nr
    #plot(vt)
    R <- order(vt,decreasing=TRUE)[1] #trouve le age de la culmination
    vt <- yt[,nc]
    inc <- c(0, diff(vt))
    #tmp <- yt[R:nr]/(P(sim)$rationPeriodMultiplier*R)  
    tmp <- vt[R:nr]/(rpm*R)    #contribution of each age class to Vm/R in m^3/ha
    tmp <- c(inc[1:R-1], tmp)
    return(list(R=R,I=inc,hVec=tmp))
  }
  # assumes yt is just a vector of volumes in 1 yr age classes starting at 1
  # Vm/R + I
  res <- lapply(ytList, cHelper,rpm=rpm)
  return(res)
}