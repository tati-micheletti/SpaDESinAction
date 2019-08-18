defineModule(sim, list(
  name = "stateVars",
  description = "keep track of stat transitions affecting multiple moduls; also classification",
  keywords = c("insert key words here"),
  authors = c(person(c("First", "Middle"), "Last", email="email@example.com", role=c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("1.1.1.9006"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit =  "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "stateVars.Rmd"),
  reqdPkgs = list("raster", "data.table", "RColorBrewer"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
    defineParameter("startTime", "numeric", 0, 0, NA, "Start time"),
    defineParameter("returnInterval", "numeric", 1, 0, NA, "waddya think?"),
    defineParameter(".plotInitialTime", "numeric", 0, NA, NA, desc="This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", 1, NA, NA, desc="This describes the simulation time at which the first plot event should occur"),
    defineParameter("persistTimes", "numeric", c(20,10,10), c(0,0,0), c(100,100,100), desc=  "For how many years do disturbances effect indicators?")
    
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "disturbanceMap", objectClass = "RasterLayer", desc = "state of burned or harvestd cells"),
    expectsInput(objectName = "spreadState", objectClass = "data.table", desc = "table of active initial spread cells"),
    expectsInput(objectName = "ageMap", objectClass = "RasterLayer", desc = "time since last disturbance"),
    expectsInput(objectName = "cutCells", objectClass = "numeric", desc = "vector of recently cut cells")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "heightMap", objectClass = "RasterLayer", desc = "crude height from age model"),
    createsOutput(objectName = "disturbanceMap", objectClass = "RasterLayer", desc = "update for caribou adjacency"),
    createsOutput(objectName = "dtMap", objectClass = "RasterLayer", desc = "timer for disturbances"),
    createsOutput(objectName = "habitatMap", objectClass = "RasterLayer", desc = "habitat map"),
    createsOutput(objectName = "harvestStateMap", objectClass = "RasterLayer", desc = "harvest")
  )
))


doEvent.stateVars = function(sim, eventTime, eventType) {
  switch (eventType,
  init = {
    sim <- Init(sim)
    
    sim <- scheduleEvent(sim, P(sim)$startTime,"stateVars","update")
    sim <- scheduleEvent(sim, P(sim)$.plotInitialTime,"stateVars","plot")
  },
  update = {
    
    sim <- Update(sim)
    sim <- scheduleEvent(sim, time(sim) + P(sim)$returnInterval, "stateVars", "update")
  },
  plot = {
    
    setColors(sim$disturbanceMap,n=4) <- c("grey80", "red", "blue", "yellow")
    Plot(sim$disturbanceMap, title = "Disturbance Map", new = TRUE)

    sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "stateVars", "plot")
  },
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

heightFromAge <-function(sim){
  x <- sim$ageMap[]/ (sim$ageMap[] + 100)
  sim$heightMap[] <- 80*x
  
  return(invisible(sim))
}

Init <- function(sim) {
  
  #use disturbanceMap as template to copy.
  #browser()
  sim$disturbanceMap <- sim$ageMap
  sim$disturbanceMap[] <- sim$ageMap[] * 0
  sim$dtMap <-  sim$ageMap
  sim$dtMap[] <-  sim$ageMap[] * 0
  
  #raster::setValues(sim$harvestStateMap, values=0) The reason CS people fucking HATE R 
  #is because of random non-orthogonality like this, and the documentation that tries to be 
  #like UNIX, but fails. We shoulda gone with Python.
  #0 = none
  #1 = burn
  #2 = cut
  #3 = adjacent to cut
  
  sim$heightMap <- raster::raster(sim$harvestStateMap)
  sim$heightMap[] <- sim$harvestStateMap[] * 0
  setColors(sim$heightMap, n=10) <- colorRampPalette(c("white","green4"))(10)
  sim<-heightFromAge(sim)
  
  return(invisible(sim))
}


### add additional events as needed by copy/pasting from above

