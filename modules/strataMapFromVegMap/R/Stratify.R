Stratify <- function(sim) {

convertTable <- data.frame(read.csv(file.path(dataPath(sim),"ConvertTable.csv")))

reclassToVernier <- matrix(c(convertTable[,2],convertTable[,2]+1,convertTable[,4]), ncol = 3)
sim$habitatMap <- raster::reclassify(x = sim$vegMap, rcl = reclassToVernier, include.lowest = TRUE, right = FALSE)

reclassToStrata <- matrix(c(convertTable[,4],convertTable[,4]+1,convertTable[,6]), ncol = 3)
sim$strataMap <- raster::reclassify(x = sim$habitatMap, rcl = reclassToStrata, include.lowest = TRUE, right = FALSE)

return(invisible(sim))
}