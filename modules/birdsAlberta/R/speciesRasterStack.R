speciesRasterStack <- function(birdModelVernier = sim$birdModelVernier,
                               rasterTemplate = sim$vegMap){
 
  species <- names(birdModelVernier)

  rasterStack <- lapply(species, FUN = function(species){
    
    rasterTemplate[] <- birdModelVernier[[species]]

    return(rasterTemplate)
  })
  
  names(rasterStack) <- species
  # rasterStack <- raster::stack(rasterStack)

  return(rasterStack)
}