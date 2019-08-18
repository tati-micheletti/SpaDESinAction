
defineModule(sim, list(
  name = "strataMapFromVegMap",
  description =  "Assign yield tables (strata) to vegetation classes", #"insert module description here",
  keywords = NA, # c("insert key words here"),
  authors = person("Tati", "Micheletti", email = "tati.micheletti@gmail.com", role = c("aut", "cre")),
  childModules = character(0),
  version = list(SpaDES.core = "0.1.1.9011", strataMapFromVegMap = "0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "strataMapFromVegMap.Rmd"),
  reqdPkgs = list(),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description"),
    defineParameter(".plotInitialTime", "numeric", 0, NA, NA, "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".initialTime", "numeric", 0, NA, NA, "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA, "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA, "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "logical", FALSE, NA, NA, "Should this entire module be run with caching activated? This is generally intended for data-type modules, where stochasticity and time are not relevant")
  ),
  inputObjects = bind_rows(
    #expectsInput("objectName", "objectClass", "input object description", sourceURL, ...),
    expectsInput(objectName = "vegMap", objectClass = "RasterLayer", desc = "Vegetation Map", sourceURL = NA)
  ),
  outputObjects = bind_rows(
    #createsOutput("objectName", "objectClass", "output object description", ...),
    createsOutput(objectName = "strataMap", objectClass = "RasterLayer", desc = "Harvest strata map"),
    createsOutput(objectName = "habitatMap", objectClass = "RasterLayer", desc = "Habitat map based on Vernier et. al 2008 [Table 1]")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.strataMapFromVegMap = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      
      sim <- Init(sim)
      sim <- Stratify(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "strataMapFromVegMap", "plot")
    },
    
    plot = {
      
      Plot(sim$strataMap, title = "Strata for Harvesting", new = TRUE)
      Plot(sim$habitatMap, title = "Habitat Classification", new = TRUE)
      
    },
    
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}


Init <- function(sim) {
  
  if (!suppliedElsewhere("vegMap", sim)){
    stop("Can't create strata maps. Provide vegetation map.")
  }
  
  return(invisible(sim))
}

.inputObjects <- function(sim){
  
  if (!suppliedElsewhere("ecodistrict", sim)){
    sim$ecodistrict <- 644
  }
  
  if (!suppliedElsewhere("studyArea", sim)){
    if (sim$ecodistrict == 644)
      message(crayon::yellow("Study area was not provided. Using ecodistrict number 644 (Primrose Lake, AB/SK)")) else
        message(crayon::yellow("Study area was not provided but ecodistrict was. Using ecodistrict number", sim$ecodistrict))
    sim$studyArea <- prepInputs(url = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip", 
                                destinationPath = inputPath(sim)) %>%
      .[.[["ECODISTRIC"]] %in% sim$ecodistrict, ]
  }
  
  if (!suppliedElsewhere("ageMap",sim)){
    message(crayon::yellow("Age map was not provided. Using Canada's age map from 2004 1km resolution"))
    sim$ageMap <- Cache(prepInputs, url = "https://drive.google.com/open?id=1WkTFn6RssKHgUBKSfAhDYcKS5Kb1Wsef", 
                        targetFile = "can_age04_1km.tif",
                        studyArea = sim$studyArea,
                        destinationPath = dataPath(sim), userTags = "objectName:ageMap")
  }
  
  if (!suppliedElsewhere("vegMap",sim)){
    message(crayon::yellow("Veg map was not provided. Using Canada's LCC2005"))
    sim$vegMap <- Cache(prepInputs, url = "ftp://ftp.ccrs.nrcan.gc.ca/ad/NLCCLandCover/LandcoverCanada2005_250m/LandCoverOfCanada2005_V1_4.zip", 
                        studyArea = sim$studyArea,
                        rasterToMatch = sim$ageMap,
                        destinationPath = dataPath(sim), userTags = "objectName:vegMap")
  }
  return(invisible(sim))
}
