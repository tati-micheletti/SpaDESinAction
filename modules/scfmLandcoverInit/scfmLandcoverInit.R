stopifnot(packageVersion("SpaDES") >= "0.99.0")

defineModule(sim,list(
    name = "scfmLandcoverInit",
    description = "Takes the LCC05 classification of 39 land cover classes, and reclassifies it to flammable and inflammable [1,0]",
    keywords = c("fire", "LCC05", "land cover classification 2005", "BEACONs"),
    childModules = character(),
    authors = c(
      person(c("Eliot", "J", "B"),"McIntire", email = "Eliot.McIntire@canada.ca", role = c("aut", "cre")),
      person("Steve", "Cumming", email = "stevec@sbf.ulaval.ca", role = c("aut")),
      person("Ian", "Eddy", email = "ian.eddy@canada.ca", role = c("aut"))),
    version = numeric_version("0.1.0"),
    spatialExtent = raster::extent(rep(NA_real_, 4)),
    timeframe = as.POSIXlt(c("2005-01-01", NA)),
    documentation = list("README.txt", "scfmLandcoverInit.Rmd"),
    timeunit = "year",
    citation = list(),
    reqdPkgs = list("fasterize", "purrr", "raster", "sf", 'rgeos',
                    "PredictiveEcology/LandR@development",
                    "PredictiveEcology/reproducible@development"),
    parameters = rbind(
      defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA, desc = "Initial time for plotting"),
      defineParameter(".plotInterval", "numeric", NA_real_, NA, NA, desc = "Interval between plotting"),
      defineParameter(".saveInitialTime", "numeric", NA_real_, NA, NA, desc = "Initial time for saving"),
      defineParameter(".saveInterval", "numeric", NA_real_, NA, NA, desc = "Interval between save events"),
      defineParameter("useCache", "logical", TRUE, NA, NA, desc = "Use cache"),
      defineParameter("neighbours", "numeric", 8, NA, NA, desc = "Number of immediate cell neighbours"),
      defineParameter("sliverThreshold", "numeric", 10000, NA, NA,
                      desc = "fire regime polygons with area less than this number will be merged")
    ),
    inputObjects = bind_rows(
      expectsInput(objectName = "ecodistrict", objectClass = "numeric", desc = "Ecodistrict number to use as studyArea"), # added
      expectsInput(objectName="nNbrs", objectClass = "numeric", desc = "raster cell neighborhood defaults to 8"), # added
      expectsInput(objectName = "studyArea", objectClass = "SpatialPolygonsDataFrame", desc = "",
                   sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"),
      expectsInput(objectName = "vegMap", objectClass = "RasterLayer", desc = "Landcover to build flammability map",
                   sourceURL = "ftp://ftp.ccrs.nrcan.gc.ca/ad/NLCCLandCover/LandcoverCanada2005_250m/LandCoverOfCanada2005_V1_4.zip"),
      expectsInput(objectName = "rasterToMatch", objectClass = "RasterLayer",
                   desc = "template raster for raster GIS operations.", 
                   sourceURL = "https://drive.google.com/open?id=1WkTFn6RssKHgUBKSfAhDYcKS5Kb1Wsef"),
      expectsInput(objectName = "fireRegimePolys", objectClass = "SpatialPolygonsDataFrame",
                   desc = "Areas to calibrate individual fire regime parameters. Defaults to ecoregions",
                   sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip"),
      expectsInput(objectName = "ageMap", objectClass = "RasterLayer",
                   desc = "Default map to define flammable map",
                   sourceURL = NA),
      expectsInput(objectName = "flammableMap", objectClass = "RasterLayer", desc = "map of landscape flammability")
    ),
    outputObjects = bind_rows(
      createsOutput(objectName = "cellsByZone", objectClass = "data.frame",
                    desc = "explains which raster cells are in which polygon"),
      createsOutput(objectName = "flammableMap", objectClass = "RasterLayer", desc = "map of landscape flammability"),
      createsOutput(objectName = "landscapeAttr", objectClass = "list", desc = "list of polygon attributes inc. area"),
      createsOutput(objectName = "landscapeAttrHARVEST", objectClass = "list", desc = "list of polygon attributes inc. area especifically for harvest"),
      createsOutput(objectName = "fireRegimeRas", objectClass = "RasterLayer",
                    desc = "Rasterized version of fireRegimePolys with values representing polygon ID")
    )
)
)

doEvent.scfmLandcoverInit = function(sim, eventTime, eventType, debug = FALSE) {
  switch(eventType,
         init = {

           sim <- Init(sim)
           sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "scfmLandcoverInit", "plot")
         },
         plot =  {
           Plot(sim$vegMap, new = TRUE, title = "Vegetation Classes")
           
           # schedule future event(s)
           sim <- scheduleEvent(sim, time(sim) + P(sim)$.plotInterval, "scfmLandcoverInit", "plot")

         },
         warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                       "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = "")))

  return(invisible(sim))
}
Init <- function(sim) {
  message("checking sim$fireRegimePolys for sliver polygons...")
  if (sf::st_is_longlat(sim$fireRegimePolys) & is.na(P(sim)$sliverThreshold)) {
    sim$fireRegimePolys <- projectInputs(x = sim$fireRegimePolys, targetCRS = crs(sim$vegMap))
    stop("You must supply P(sim)$sliverThreshold for fireRegimePolys with an unprojected CRS")
  }
  sim$fireRegimePolys$trueArea <- gArea(sim$fireRegimePolys, byid = TRUE)
  if (is.na(P(sim)$sliverThreshold)) {
    sim@params[[currentModule(sim)]]$sliverThreshold <- 1e4

  }
  if (any(sim$fireRegimePolys$trueArea < P(sim)$sliverThreshold)) {
    message("sliver polygon(s) detected. Merging to their nearest valid neighbour")
  sim$fireRegimePolys <- Cache(deSliver, sim$fireRegimePolys, threshold = P(sim)$sliverThreshold,
                               userTags = c("deSliver", currentModule(sim)))
  }
  if (is.null(sim$fireRegimePolys$PolyID)) {
    if (is.null(sim$fireRegimePolys$REGION_)) {
      sim$fireRegimePolys$PolyID <- row.names(sim$fireRegimePolys)
    } else {
      sim$fireRegimePolys$PolyID <- as.numeric(sim$fireRegimePolys$REGION_)
    }
  }
  temp <- sf::st_as_sf(sim$fireRegimePolys)
  temp$PolyID <- as.numeric(temp$PolyID) #fasterize needs numeric; row names must stay char
  sim$fireRegimeRas <- fasterize::fasterize(sf = temp, raster = sim$vegMap, field = "PolyID")
  sim$flammableMap <- LandR::defineFlammable(sim$vegMap, filename2 = NULL)
  sim$flammableMap <- setValues(raster(sim$flammableMap), sim$flammableMap[])
  sim$flammableMap <-  mask(sim$flammableMap, mask = sim$fireRegimeRas)
  
  # This makes sim$landscapeAttr & sim$cellsByZone
  outs <- Cache(genFireMapAttr,
                sim$flammableMap,
                sim$fireRegimePolys,
                P(sim)$neighbours)

  sim$landscapeAttr <- outs$landscapeAttr
  
  sim$landscapeAttrHARVEST <- genFireMapAttrHARVEST(flammableMap = sim$flammableMap, 
                          nNbrs = sim$nNbrs) # TO COMPARE WITH THE OTHER OUT
  
  sim$cellsByZone <- outs$cellsByZone
  return(invisible(sim))
}

genFireMapAttr <- function(flammableMap, fireRegimePolys, neighbours) {
  #calculate the cell size, total area, and number of flammable cells, etc.
  #All areas in ha
  cellSize <- prod(res(flammableMap)) / 1e4 # in ha
  

  if (neighbours == 8)
    w <- matrix(c(1, 1, 1, 1, 0, 1, 1, 1, 1), nrow = 3, ncol = 3)
  else if (neighbours == 4)
    w <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3, ncol = 3)
  else
    stop("illegal neighbours specification")

  makeLandscapeAttr <- function(flammableMap, weight, fireRegimePolys) {
    neighMap <- focal(x = flammableMap, w = w, na.rm = TRUE) #default function is sum(...,na.rm)

    # extract table for each polygon
    valsByPoly <- raster::extract(neighMap, fireRegimePolys, cellnumbers = TRUE)
    valsByPoly <- lapply(valsByPoly, na.omit)
    names(valsByPoly) <- fireRegimePolys$PolyID
    uniqueZoneNames <- fireRegimePolys$PolyID #get unique zones.
    valsByZone <- lapply(uniqueZoneNames, function(ecoName) {
      aa <- valsByPoly[names(valsByPoly) == ecoName]
      if (is.list(aa))
        aa <- do.call(rbind, aa)
      return(aa)
    })
    names(valsByZone) <- uniqueZoneNames

    # Derive frequency tables of number of flammable cells, per polygon type, currently ECOREGION
    nNbrs <- lapply(valsByZone, function(x) {
      nNbrs <- tabulate(x[, 2] + 1, 9)#depends on sfcmLandCoverInit
      names(nNbrs) <- 0:8
      return(nNbrs)
    })
    nFlammable <- lapply(valsByZone, function(x) {
      sum(getValues(flammableMap)[x[, 1]], na.rm = TRUE) #sums flammable pixels in FRI polygons
    })

    landscapeAttr <- purrr::transpose(list(cellSize = rep(list(cellSize), length(nFlammable)),
                                           nFlammable = nFlammable,
                                           nNbrs = nNbrs,
                                           cellsByZone = lapply(valsByZone, function(x) x[, 1])))

      landscapeAttr <- lapply(landscapeAttr, function(x) {
        append(x, list(burnyArea = x$cellSize * x$nFlammable))
      })
      names(landscapeAttr) <- names(valsByZone)

    return(landscapeAttr)
  }

  landscapeAttr <- makeLandscapeAttr(flammableMap, w, fireRegimePolys)

  cellsByZoneFn <- function(flammableMap, landscapeAttr) {

    cellsByZone <- data.frame(cell = 1:ncell(flammableMap), zone = NA_character_, stringsAsFactors = FALSE)

    for (x in names(landscapeAttr)) {
      cellsByZone[landscapeAttr[[x]]$cellsByZone, "zone"] <- x
      }
    return(cellsByZone)
  }

  cellsByZone <- cellsByZoneFn(flammableMap, landscapeAttr)

  return(invisible(list(landscapeAttr = landscapeAttr, cellsByZone = cellsByZone)))
}

### template initilization

.inputObjects <- function(sim) {
  invisible(tryCatch(dev.off(), error = function(e){}))
  dev.new()
  dev.new()
  quickPlot::clearPlot()
  
  dPath <- dataPath(sim) #where files will be downloaded
  cacheTags = c(currentModule(sim), "function:.inputObjects")
  
  if (!suppliedElsewhere("ecodistrict", sim)) {
  sim$ecodistrict <- 644
  }
  
  if (!suppliedElsewhere("studyArea", sim)) {
    message("study area not supplied. Using Ecodistrict 644")
    sim$studyArea <- Cache(prepInputs, url = extractURL(objectName = "studyArea", sim = sim), 
                       destinationPath = getPaths()$inputPath, targetCRS = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84", filename2 = "studyArea", overwrite = TRUE) %>%
      .[.[["ECODISTRIC"]] %in% sim$ecodistrict, ]
  }

  if (!suppliedElsewhere("rasterToMatch", sim)) {
    message("rasterToMatch not supplied. generating from Canada Forest's Age Map using cropped to studyArea")
    sim$rasterToMatch <- Cache(prepInputs, url = extractURL(objectName = "rasterToMatch", sim = sim), 
                           targetFile = "can_age04_1km.tif",
                           studyArea = sim$studyArea,
                           destinationPath = getPaths()$inputPath, userTags = "objectName:RTM", 
                           filename2 = NULL)
  }

  if (!suppliedElsewhere("vegMap", sim)) {
    message("vegMap not supplied. Using default LandCover of Canada 2005 V1_4a")

    sim$vegMap <- LandR::prepInputsLCC(year = 2005,
                                destinationPath = dPath,
                                studyArea = sim$studyArea,
                                rasterToMatch = sim$rasterToMatch,
                                filename2 = TRUE,
                                overwrite = TRUE,
                                userTags = c("cacheTags", "vegMap"))
  }

  if (!suppliedElsewhere("fireRegimePolys", sim)) {
    message("fireRegimePolys not supplied. Using default ecoregions of Canada cropped to study Area")
    sim$fireRegimePolys <- Cache(prepInputs, url = extractURL("fireRegimePolys", sim),
                                      destinationPath = dPath,
                                      studyArea = sim$studyArea,
                                      rasterToMatch = sim$rasterToMatch,
                                      filename2 = TRUE,
                                      overwrite = TRUE,
                                      userTags = c("cacheTags", "fireRegimePolys"))

  }
  
  if (!suppliedElsewhere("ageMap",sim)){
    message(crayon::yellow("Age map was not provided. Using Canada's age map from 2004 1km resolution"))
    sim$ageMap <- Cache(prepInputs, url = "https://drive.google.com/open?id=1WkTFn6RssKHgUBKSfAhDYcKS5Kb1Wsef", 
                        targetFile = "can_age04_1km.tif",
                        studyArea = sim$studyArea,
                        destinationPath = dataPath(sim), userTags = "objectName:ageMap", 
                        filename2 = NULL)
  }
  
  if (!suppliedElsewhere("flammableMap", sim)){
    sim$flammableMap <- sim$ageMap
    names(sim$flammableMap) <- "flammableMap"
  }
  
  if (!suppliedElsewhere("sliverThreshold", sim)){
    sim$sliverThreshold <- 10000
  }
  if (!suppliedElsewhere("numTypes", sim)){
    sim$numTypes <- 8
  }
  if (!suppliedElsewhere("landscapeAttr", sim)){
    sim$landscapeAttr = list(cellSize = 6.25)
  }
  if (!suppliedElsewhere("nNbrs", sim)){
    sim$nNbrs <- 8
  }
  if (!suppliedElsewhere("areaInHa", sim)){
    sim$areaInHa = pi
  }

  return(invisible(sim))
}
