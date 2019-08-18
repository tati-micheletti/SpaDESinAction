
defineModule(sim, list(
  name = "birdsAlberta",
  description = "This bird module is an adaptation of the work of Vernier et al. 2008 for 2 bird species: COWA and RBNU",
  keywords = NA, # c("insert key words here"),
  authors = person("Tati", "Micheletti", email = "tati.micheletti@gmail.com", role = c("aut", "cre")),
  childModules = character(0),
  version = list(SpaDES.core = "0.1.1.9010", birdsAlberta = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "birdsAlberta.Rmd"),
  reqdPkgs = list("data.table", "reproducible"),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter("modelTime", "numeric", 0, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInitialTime", "numeric", 0, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".plotInterval", "numeric", end(sim)-start(sim) , NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "covarParams", objectClass = "list", desc = "Table with model parameters", sourceURL = NA),
    expectsInput(objectName = "vegMap", objectClass = "RasterLayer", desc = "vegetation map", sourceURL = sim$url.vegMap),
    expectsInput(objectName = "covar", objectClass = "data.table", desc = "Table with covariate values", sourceURL = NA)
    
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "birdModelVernier", objectClass = "list", desc = "Lists with table containing cell ID and abundance value per species"),
    createsOutput(objectName = "covarTable", objectClass = "list", desc = "List of covariates for all years"),
    createsOutput(objectName = "birdAbundance", objectClass = "list", desc = "List of rasters presenting bird probability of presence"),
    createsOutput(objectName = "covar", objectClass = "list", desc = "Table with covariate values")
    )
))

doEvent.birdsAlberta = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {

      sim <- Init(sim)

      # schedule future event(s)
      
      sim <- scheduleEvent(sim, P(sim)$modelTime, "birdsAlberta", "model")
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "birdsAlberta", "plot")
    },
    
    model = {

      sim$covar <- getLocalCovars(covarTable = sim$covar,
                                  disturbanceMap = sim$disturbanceMap,
                                  ageMap = sim$ageMap,
                                  habitatMap = sim$habitatMap,
                                  vegMap = sim$vegMap)
      sim$covar <- getNeighborhoodCovars(covarTable = sim$covar, 
                                         disturbanceMap = sim$disturbanceMap,
                                         ageMap = sim$ageMap,
                                         habitatMap = sim$habitatMap,
                                         vegMap = sim$vegMap)
      
      sim$covarTable[[paste0("YEAR",time(sim))]] <- sim$covar      
      
      sim$birdModelVernier <- Cache(birdModelVernier, covarTable = sim$covar,
                                               covarParams = sim$covarParams)
      # schedule future event(s)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "birdsAlberta", "model")
    },
    
    plot = {
      
      sim$birdAbundance[[paste0("Year", time(sim))]] <- speciesRasterStack(birdModelVernier = sim$birdModelVernier,
                                              rasterTemplate = sim$vegMap)
      
        Plot(sim$birdAbundance[[paste0("Year", time(sim))]][["COWA"]], title = "COWA presence probability")
        Plot(sim$birdAbundance[[paste0("Year", time(sim))]][["RBNU"]], title = "RBNU presence probability")

      # schedule future event(s)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "birdsAlberta", "plot")
      
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

Init <- function(sim) {


  
  return(invisible(sim))
}

.inputObjects <- function(sim) {
 
  if(!suppliedElsewhere("covarParams", sim)){
    
     sim$covarParams <- invisible(createParams())
  }
  
  if(suppliedElsewhere("habitatMap",sim)){
    
    sim$covar <- data.table::data.table(matrix(0,
                                               ncol = length(colnames(sim$covarParams[[1]])[2:length(colnames(sim$covarParams[[1]]))]), 
                                               nrow = ncell(sim$habitatMap)))
  } else {
    
        if(suppliedElsewhere("vegMap",sim)){
        sim$covar <- data.table::data.table(matrix(0, 
                                                   ncol = length(colnames(sim$covarParams[[1]])[2:length(colnames(sim$covarParams[[1]]))]), 
                                                   nrow = ncell(sim$vegMap))) 
        } else stop("Please, provide a vegetation or habitat map.") }
  
  names(sim$covar) <- colnames(sim$covarParams[[1]])[2:length(colnames(sim$covarParams[[1]]))]

  return(invisible(sim))
}

