
defineModule(sim, list(
  name = "caribouAlberta",
  description = "Caribou population model. Calculates changes in boreal caribou herd size based on Sorensen model (2008)",
  keywords = c("caribou", "population dynamics"),
  authors = person("Antoine", "Adde", email = "antoine.adde.1@ulaval.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = list(SpaDES.core = "0.1.1.9011", caribouAlberta = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "caribouAlberta.Rmd"),
  reqdPkgs = list("ggplot2"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("N0", "numeric", NA_real_, 0, 5000, "Initial herd size"),
    defineParameter(".statsInitialTime", "numeric", 0, NA_real_, NA_real_, "This describes the simulation time at which the first stat event should occur"),
    defineParameter(".statsInterval", "numeric", 1, NA_real_, NA_real_, "This describes the simulation time interval between stat events"),
    defineParameter(".plotInterval", "numeric", 1, NA_real_, NA_real_, "This describes the simulation time interval between stat events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "maxCaribPop", objectClass ="numeric", desc = "Max Caribou population. Defaults is 200"),
    expectsInput(objectName="disturbanceMap", objectClass = "RasterLayer", desc="state of burned or harvestd cells"),
    expectsInput(objectName="landscapeAttr", objectClass = "list", desc="landscape parameters from scfm {Regime|Driver}")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "SorensenStats", objectClass ="data.frame", desc = "lambda and herd size table"),
    createsOutput(objectName = "Nt", objectClass ="numeric", desc = "Caribou population size in time t"),
#    createsOutput(objectName = "caribouPopGrowth", objectClass ="data.frame", desc = "Caribou population growth graph"),
createsOutput(objectName = "plotCaribou", objectClass ="ggplot", desc = "Caribou population/lambda over time")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.caribouAlberta = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      
      sim <- StatsInit(sim)
      
      sim <- scheduleEvent(sim, P(sim)$.statsInitialTime, "caribouAlberta", "stats")  #get the initial values
      sim <- scheduleEvent(sim, P(sim)$.statsInitialTime, "caribouAlberta", "plot", eventPriority = 10)  #get the initial values
      
    },
    stats = {
      sim$SorensenStats <- Stats(disturbanceMap = sim$disturbanceMap, 
                   landscapeAttr = sim$landscapeAttr,
                   maxCaribPop = sim$maxCaribPop,
                   SorensenStats = sim$SorensenStats)
      
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.statsInterval, "caribouAlberta", "stats")
      
    },
    plot = {
      
      sim$plotCaribou <- plotCaribou(sim = sim, SorensenStats = sim$SorensenStats)
      
 #     Plot(sim$plotCaribou$populationCaribou, title = "Caribou population size", new = TRUE)
      Plot(sim$plotCaribou, title = "Caribou growth rate dynamics", new = TRUE)
      
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "caribouAlberta", "plot")
      
    },

    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  
  return(invisible(sim))
}


.inputObjects <- function(sim) {
  if (!suppliedElsewhere("maxCaribPop", sim)){
    sim$maxCaribPop <- 200
  }
  
  return(invisible(sim))
}