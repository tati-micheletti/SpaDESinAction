birdModelVernier <- function(covarTable = sim$covar,
                             covarParams = sim$covarParams){
  
  species <- names(covarParams)
  pixels <- 1:nrow(covarTable)
  
lp.species <-  lapply(X = species, FUN = function(species) {
  lp.pixel <- sapply(X = pixels, FUN = function(pixels){
  
linearPredictor <-  covarParams[[species]]$Constant + # Constant
                    covarParams[[species]]$L_CUT  * covarTable$L_CUT[pixels]  +  # Station in recently harvested cutblock (≤30 years) proportion ==> HarvestStateMap > 0
                    covarParams[[species]]$L_YDEC * covarTable$L_YDEC[pixels] + # Station in young deciduous forest (31–90 years) proportion ==> Age map * LCC05
                    covarParams[[species]]$L_ODEC * covarTable$L_ODEC[pixels] +  # Station in old deciduous forest (>90 years) proportion ==> Age map = LCC05
                    covarParams[[species]]$L_MIX  * covarTable$L_MIX[pixels]  +  # Station in mixed deciduous/conifer forest (>30 years) proportion ==> Age map * LCC05
                    covarParams[[species]]$L_CC   * covarTable$L_CC[pixels]   +  # Area weighted mean forest crown closure ==> 0.3 low, 0.6 medium, 0.8 closed
                    covarParams[[species]]$L_WDIS * covarTable$L_WDIS[pixels] +  # Distance from station centre to nearest river or lake ==> ??? just 0 for all?
                    covarParams[[species]]$N_CUT  * covarTable$N_CUT[pixels]  +  # Station in recently harvested cutblock (≤30 years) proportion ==> HarvestStateMap > 0
                    covarParams[[species]]$N_LATE * covarTable$N_LATE[pixels] +  # Proportion in old forest (90* year) proportion ==> Age map (independently of class?)
                    covarParams[[species]]$N_DEC  * covarTable$N_DEC[pixels]  +  # Proportion in deciduous forest ==> LCC05 (independently of age?)
                    covarParams[[species]]$N_MIX  * covarTable$N_MIX[pixels]  +  # Station in mixed deciduous/conifer forest (>30 years) proportion ==> Age map
                    covarParams[[species]]$N_SB   * covarTable$N_SB[pixels]   +  # Presence of black spruce forest ==> Y/N
                    covarParams[[species]]$N_RICH * covarTable$N_RICH[pixels]  # Number of habitat classes in neighbourhood ==> LCC05

p <- exp(linearPredictor)/1 + exp(linearPredictor)

    return(p)
  })
  return(lp.pixel)
})
names(lp.species) <- species

  return(invisible(lp.species))
}
