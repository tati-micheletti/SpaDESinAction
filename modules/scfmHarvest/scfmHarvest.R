
defineModule(sim, list(
  name = "scfmHarvest",
  description = NA, #"insert module description here",
  keywords = NA, # c("insert key words here"),
  authors = person("Steve", "Cumming", email = "first.last@example.com", role = c("aut", "cre")),
  childModules = character(0),
  version = list(SpaDES.core = "0.1.1.9011", scfmHarvest = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "scfmHarvest.Rmd"),
  reqdPkgs = list(),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("minVolume", "numeric", 50, 0, NA, "Minimum volume in m^3/ha for eligibility to harvest"),
    defineParameter("minAge", "numeric", 40, 0, NA, "Minimum age to harvest"),
    defineParameter("startTime", "numeric", 1, 0, NA, "Start time"),  # IF BREAKS, BACK TO 0 and change birds
    defineParameter("returnInterval", "numeric", 1, 0, 0, "waddya think?"),
    defineParameter("minBlockSize", "numeric", 20, 1, NA, "minimum block size (ha)"),
    defineParameter("greenUpPeriod", "numeric", 30, 0, 30, "For how many years do cut blocks affect adjacent areas?"),
    defineParameter(".plotInitialTime", "numeric", 0, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Caching the module")
  ),
  inputObjects = bind_rows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = "ageMap", objectClass = "RasterLayer", desc = "forest Age structure" ),
    expectsInput(objectName = "yieldTables", objectClass = "list", desc = "what it is"),
    expectsInput(objectName = "strataMap", objectClass = "Rasterayer", desc = "map assigning cells to yield tables"),
    expectsInput(objectName = "landscapeAttr", objectClass = "list", desc = "info about landscape we are running on; e.g. cellSize"),
    expectsInput(objectName = "annualCut", objectClass = "list", desc = "m^3/yr per yield Class")
    
  ),
  outputObjects = bind_rows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = "harvestStateMap", objectClass = "RasterLayer", desc = "areas locked up by spatial constraints"),
    createsOutput(objectName = "volMap", objectClass = "RasterLayer", desc = "volume denisty per cell"),
    createsOutput(objectName = "harvestStats", objectClass = "data.frame" , desc = "record of harvest activity by strata"),
    createsOutput(objectName = "cutCells", objectClass = "numeric", desc = "vector of cells harvested in current period")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.scfmHarvest = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {

      sim <- Init(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$startTime, "scfmHarvest", "harvest")
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "scfmHarvest", "plot")
    },
    harvest = {
      
      sim <- Harvest(sim) 
      
      sim <- scheduleEvent(sim, time(sim) + P(sim)$returnInterval,"scfmHarvest","harvest")
    },
    plot = {
      
      # Plot(sim$harvestStateMap, legendRange=0:1, title = "Harvest State Map", new = TRUE)
      
      sim <- scheduleEvent(sim, time(sim)+P(sim)$.plotInterval, "scfmHarvest", "plot")
    },
  
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}


Init <- function(sim) {

  sim$harvestStats <- data.frame(rep(list(numeric(0)),length(sim$yieldTables[[1]]))) #again, use the 1st and only region
  
  #sim$disturbanceMap<-raster(sim$strataMap) * 0 
  #setColors(sim$disturbanceMap, n=3) <- c("white", "blue","red")
  
  if(!suppliedElsewhere("strataMap", sim)){stop("No strata map provided. This would normally come from module 'vegMapToStrataMap")}
  sim$harvestStateMap <- raster::raster(sim$strataMap)
  sim$harvestStateMap[] <- sim$strataMap[]*0
  setColors(sim$harvestStateMap, n=10) <- c("white","black")
  
  sim$volMap <- raster::raster(sim$strataMap) 
  sim$volMap[] <- sim$strataMap[]*0
  sim$merchCells <- list()    #of cells currently in the list for each stratum
  sim$cutCells <- numeric()
  
  
  return(invisible(sim))
}

.inputObjects <- function(sim) {

  if(!suppliedElsewhere("strataMap", sim)){
    stop("No strata map provided. This would normally come from module 'vegMapToStrataMap") # [ IMPROVE ] add a strataMap for the module to work alone
  }
  return(invisible(sim))
}