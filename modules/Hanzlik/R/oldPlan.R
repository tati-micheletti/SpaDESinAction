oldPlan <- function(sim) {

  # browser()
  # 
  # region <- sim$hanzlikPars[[1]]   # we are going to use the first region, if we can
  # 
  # hVec<-sim$hanzlikPars[[1]][["SwAw"]]$hVec       #1 == "AB": for now, assume hanzlikPars has at least one object, and take the 1st
  # #Comment peut-on le generaliser pour plusieurs courbes de rendemment?
  # nh <- length(hVec)
  # res<-rep(0,nh)
  # x<-tabulate(sim$ageMap[])         # yieldTableMazAge = nh
  # # tabulate fills gaps with 0s
  # nx<-length(x)
  # if (nx <= nh){
  #   res[1:nx] <- x
  # }
  # else {
  #   res <- x[1:nh]
  #   x[nh] <- x[nh] + sum(x[nh+1:nx]) #accumulate any "missing ones"
  # }
  # res<-sum(res*hVec) * sim$cellSize #this needs to be changed.
  # print(sprintf("%d AAC = %5.1f\n",time(sim),res/1e3))
  # sim$annualCut <- res
  return(invisible(sim))
}

