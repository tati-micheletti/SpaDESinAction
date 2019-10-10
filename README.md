---
title: "SpaDESinAction"
author: "Tati Micheletti"
date: "October 10, 2019"
output:
  pdf_document: default
  html_document: default
---

## Integrated simulation models: fire, harvesting, birds and caribou

The present toy model is an exercise to understand how SpaDES works when integrating different modules. 
It is composed of 4 main models (with one to several submodules within each):  

1. SCFM (Steve Cumming's fire model): this model is composed of 3 main modules (ignition, escape and spread), and 3 auxiliary modules to establish the landscape, fire drivers, and fire regime.
  See more information on this model in  [modules/harvestREADME.Rmd](modules/harvestREADME.Rmd)    

2. SCSHM (Steve Cumming's Simple Harvest model): this model is composed of 5 main modules:
  `loadYieldTables` - loads yield tables necessary for identifying cut blocks;
  `Hanzlik` - determines AAC;
  `strataMapFromVegMap` - creates a strata map based on the vegetation;
  `scfmHarvest` - the harvest component itself;
  `stateVars` - a module that updates the state variables (what was cut and what burned).
  It also has an auxiliary module to establish the landscape. See more information on this model in [modules/scfmREADME.Rmd](modules/scfmREADME.Rmd)  

3. birdsAlberta (bird model based on Vernier's et al., 2008): this model uses coefficients from Vernier et al., 2008 to predict abundance of 2 species of bird (COWA and RBNU) in function of disturbance (fire and harvest).  

4. caribouAlberta (caribou model based on Sorensen's et al., 2008): this model uses coefficients from Sorensen et al., 2008 to predict caribou lambda (growth rate) in function of disturbance (fire and harvest).  

### Running the simulations

#### 1. Setting up all necessary packages, paths and loading all necessary libraries  

When setting up the paths using the `SpaDES.core` function `setPaths()`, note that the modules folder should point to where all your modules' *folders* are.
`SpaDES` will find all of these by their names.

```{r setup, echo = TRUE, message = FALSE}
devtools::install_github("PredictiveEcology/SpaDES.core", ref = "development")
devtools::install_github("PredictiveEcology/SpaDES.tools", ref = "development")
devtools::install_github("PredictiveEcology/reproducible", ref = "development")
devtools::install_github("PredictiveEcology/reproducible", ref = "development")
devtools::install_github("tati-micheletti/usefun")

# Set a few options
options(spades.lowMemory = TRUE)
options(spades.moduleCodeChecks = FALSE)
 
# Load libraries  
  library("SpaDES")
  library("magrittr")
  library("raster")
  
# Set the directories
  workDirectory <- getwd()

    SpaDES.core::setPaths(
           modulePath = file.path(workDirectory, "modules"),
           inputPath = file.path(workDirectory, "inputs"), 
           outputPath = file.path(workDirectory, "outputs"),
           cachePath = file.path(workDirectory, "cache"))
```

#### 2. Setting up all inputs for the simulation runs

To setup the simulation, we will set 3 'inputs' for the model:  

1. **TIME**: start by setting the time for how long you want to run the simulations for.
  This needs to be a list, in the format of `list(start = 0, end = 10)`.  

2. **PARAMETERS**: Next, pass the values of any parameters you would like to change.
  These values should be a `list()` of the name of the module that contains these.
  For example, if you want to modify how often the bird model should be run (by default only happening in the beginning and end of the simulation), you can pass inside the parameters list the following: `birdsAlberta = list(.plotInterval = 5)`, for it to run every 5 years.
  You can find all parameters that can be changed in the metadata of each module, together with a default value and the description of what that parameter does.  

3. **OBJECTS**: At last, pass any objects you might want to override from the defaults.
  One example is the study area.
  You can load a shapefile in R, and pass it in the form `list(studyArea = shapefile)`.
  The names of all objects passed need to match the names these have in the modules.
  These can also be found in the metadata of each module.  

It is interesting to note that while **objects are shared** accross modules, **parameters are not**.
The modules are in fact integrated by the objects these use as *inputs* and objects these produce as *outputs*.
These also have to be defined in the metadata.  

Module developers should provide default values for all parameters and objetcs, even if the objects are expected as inputs.
This guarantees that anyone running your module, even without any data, is capable of seeing it work.

```{r setupSimulation, message = FALSE}
times <- list(start = 0, end = 10)
parameters <- list()
objects <- list()
```

#### 3. Setting up the modules

**MODULES**: The last step to run the simulation is to identify which modules you want to run together. This should also be passed as a list of the modules names i.e. list("module1", "module2", "module3").

Here, we will mix and match a set of different modules to demonstrate the power of SpaDES.
First, we will run a fire model, followed by a harvesting model, and at last, an integrated model that has fire, harvesting, birds and caribou.

```{r modules, message = FALSE}
modulesFire <- list("scfmLandcoverInit", "scfmRegime", 
                    "scfmDriver", "scfmIgnition", 
                    "scfmEscape", "scfmSpread")

modulesHarvest <- list("scfmLandcoverInit", "loadYieldTables", 
                       "Hanzlik", "strataMapFromVegMap", 
                       "scfmHarvest", "stateVars") 

allModules <- list("scfmLandcoverInit", "scfmRegime", 
                   "scfmDriver", "scfmIgnition", 
                   "scfmEscape", "scfmSpread", # Fire model
                   "loadYieldTables", "Hanzlik", 
                   "strataMapFromVegMap", "scfmHarvest", # Harvesting
                   "stateVars", # Updating maps and variables
                   "birdsAlberta", # Bird model
                   "caribouAlberta") # Caribou population model
```

#### 4. `SpaDES` in ACTION!

At last, we will run the simulation calling the `simInitAndSpades()` function.
This function needs as inputs the paths, times, parameters, objects, and modules, apart from the *loadOrder* of the modules.
If the load order is not provided, the modules might be loaded in a different order than these are supposed to, which might return an error due to missing objects.

##### **FIRE**

First, we will run the fire model:

```{r fireModel, message = FALSE, warning = FALSE}
FIREinAction <- simInitAndSpades(times = times,
                                   params = parameters,
                                   modules = modulesFire,
                                   objects = objects,
                                   paths = getPaths(),
                                   loadOrder = unlist(modulesFire))
```

```{r, include = FALSE}
dev.off() # This closes the device opened, as we open a new one when we run the next call
```

Note that you can see all events happening, as well as any messages printed during the simulation.
Really handy is also to see all plots while the simulation is running.
This is really helpful to visually identify potential problems even before the simulation ends.

Now we will check some of the results from this simulation.
First, lets check the age map, flammable map (where fire can burn), the last years' fires and the cummulative fire map:

```{r fireResultsHid, include = FALSE}
# Fire plots
clearPlot()
Plot(FIREinAction$ageMap, title = "Age Map")
Plot(FIREinAction$flammableMap, title = "Flammable Map", col = c("skyblue", "red"))
Plot(FIREinAction$rstCurrentBurn, title = "Last year's burns", col = c("lightgrey", "red"))
Plot(FIREinAction$burnMap, title = "Cummulative burn map")
```

We can also see some of the fire statistics such as:

1. Summary of burns

```{r fireSummary1}
# Fire Summary
knitr::kable(FIREinAction$burnSummary)
```

2. Probability of ignition of a pixel

```{r fireSummary2}
# Fire Summary
print(FIREinAction$pIg)
```

3. Fire driver parameters

```{r fireSummary3}
# Fire Summary
FIREinAction$scfmDriverPars$`108`
```

##### HARVEST

Now, we will run the harvest model:

```{r harvestModel, message = FALSE, warning = FALSE}
HARVESTinAction <- simInitAndSpades(times = times, 
                                   params = parameters, 
                                   modules = modulesHarvest, 
                                   objects = objects, 
                                   paths = getPaths(),
                                   loadOrder = unlist(modulesHarvest))
```

```{r, include = FALSE}
dev.off()
```

Now we will check some of the results from this simulation.
First, lets check the age map, disturbance map (where all the harvest was made), and the strata for harvesting:

```{r harvestResults, include = FALSE}
# Harvest plots
clearPlot()
Plot(HARVESTinAction$ageMap, title = "Age Map")
Plot(HARVESTinAction$disturbanceMap, title = "Disturbance Map")
Plot(HARVESTinAction$strataMap, title = "Strata for Harvesting")
Plot(HARVESTinAction$volMap, title = "Volume for Harvesting")
```

We can also see some of the fire statistics such as harvest statistics (total harvested per strata)

```{r harvestSummary1, echo = FALSE}
# Harvest Summary
names(HARVESTinAction$harvestStats) <- paste0("Strata", 1:ncol(HARVESTinAction$harvestStats))
knitr::kable(HARVESTinAction$harvestStats)
```

##### INTEGRATED MODEL

Now, we will integrate the fire, harvesting, birds and caribou:

```{r integratedModel, message = FALSE, warning = FALSE}
SpaDESinAction <- simInitAndSpades(times = times, 
                                   params = parameters, 
                                   modules = allModules, 
                                   objects = objects, 
                                   paths = getPaths(), 
                                   loadOrder = unlist(allModules))
```

```{r, include = FALSE}
dev.off()
```

Now we will check some of the results from this simulation. First, lets check the birds:

```{r birdsResults, include = FALSE}
# Birds plots
clearPlot()
Plot(SpaDESinAction$birdAbundance[[paste0("Year", times$start)]][["COWA"]], 
     title = paste0("COWA presence probability year ", times$start))
Plot(SpaDESinAction$birdAbundance[[paste0("Year", times$end)]][["COWA"]], 
     title = paste0("COWA presence probability year ", times$end))
Plot(SpaDESinAction$birdAbundance[[paste0("Year", times$start)]][["RBNU"]], 
     title = paste0("RBNU presence probability year ", times$start))
Plot(SpaDESinAction$birdAbundance[[paste0("Year", times$end)]][["RBNU"]], 
     title = paste0("RBNU presence probability year ", times$end))
```

We can even make some cool analysis with these results, such as answering how many birds were lost in the 10 years of the simulation:

```{r birdsStats, include = FALSE}
# bird delta plots
clearPlot()
library("RColorBrewer")
cols <- brewer.pal(5, "PRGn")
pal <- colorRampPalette(cols)
COWAChange <- (SpaDESinAction$birdAbundance[[paste0("Year", times$end)]][["COWA"]] - 
                 SpaDESinAction$birdAbundance[[paste0("Year", times$start)]][["COWA"]]/
                 SpaDESinAction$birdAbundance[[paste0("Year", times$start)]][["COWA"]])
Plot(COWAChange, title = paste0("COWA Delta presence probability from year ", 
                                times$start, " to ", times$end), cols = pal(3))

RBNUChange <- (SpaDESinAction$birdAbundance[[paste0("Year", times$end)]][["RBNU"]] -
                 SpaDESinAction$birdAbundance[[paste0("Year", times$start)]][["RBNU"]]/
                 SpaDESinAction$birdAbundance[[paste0("Year", times$start)]][["RBNU"]])
```

```{r, echo = FALSE}
Plot(RBNUChange, title = paste0("RBNU Delta presence probability from year ", 
                                times$start, " to ", times$end), cols = pal(3))
```

We can also see some of the caribou outputs:

1. Lambda through time

```{r caribouPlot, echo = FALSE}
# caribou plot
quickPlot::clearPlot()
Plot(SpaDESinAction$plotCaribou, title = "Caribou lambda through time")
```

2. Lambda and population size through time

```{r caribouStats, echo = FALSE}
# caribou Summary
knitr::kable(SpaDESinAction$SorensenStats)
```

### Hands-on

Are you ready to start using `SpaDES`?
You can either play around with these modules (change parameters or change the study the area) or you can start you very own module.
For that, you can simply type:

```{r newMod, eval = FALSE, echo = TRUE}
newModule("ModuleName", path = getPaths()$modulePath)
```

**Happy `SpaDES`-ing!!!**
