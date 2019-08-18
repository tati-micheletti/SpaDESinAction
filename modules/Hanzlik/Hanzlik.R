defineModule(sim, list(
  name = "Hanzlik",
  description = "insert module description here",
  keywords = c("insert key words here"),
  authors = c(person(c("First", "Middle"), "Last", email="email@example.com", role=c("aut", "cre"))),
  childModules = character(),
  version = numeric_version("1.1.1.9001"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit =  "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Hanzlik.Rmd"),
  reqdPkgs = list(),
  parameters = rbind(
    #defineParameter("paramName", "paramClass", value, min, max, "parameter description")),
    defineParameter("rationPeriodMultiplier", "numeric", 2, NA, NA, "Vm/R*x"),
    defineParameter("replanInterval", "numeric", 10, NA, NA, "How frequently do we replan?")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = "yieldTables", objectClass = "list", desc = "a list of named yield table objects"),
    expectsInput(objectName = "ageMap", objectClass = "RasterLayer", desc = "map from whence comes theforest age-class structure"),
    expectsInput(objectName = "strataMap", objectClass = "Rasterayer", desc = "map assigning cells to yield tables"),
    expectsInput(objectName = "landscapeAttr", objectClass = "list", desc = "info about landscape we are running on; e.g. cellSize")
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = "hanzlikPars", objectClass = "list", desc = "hanzlik structures per yield class"),
    createsOutput(objectName = "annualCut", objectClass = "list", desc = "m^3/yr per yield Class")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.Hanzlik = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      sim <- Init(sim)
      sim <- scheduleEvent(sim, time(sim), "Hanzlik", "plan")
      },
    plan = {
      sim <- Plan(sim)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$replanInterval, "Hanzlik", "plan")
    },
    warning(paste("Undefined event type: '", events(sim)[1, "eventType", with = FALSE],
                  "' in module '", events(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

Init <- function(sim) {
  rpm <- P(sim)$rationPeriodMultiplier
  sim$hanzlikPars <- lapply(sim$yieldTables, calcHanzlik, rpm = rpm)
  return(invisible(sim))
}
