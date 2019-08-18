
VernierTable1 <- data.frame(Class = seq(1:11), Description = c( 
"Water (lakes and rivers)",
"Agriculture, non-forested openings, and wetland",
"Fire-origin stands ≤30 years and >70% deciduous",
"Stands 31–90 years and >70% deciduous",
"Stands >90 years and >70% deciduous",
"Stands with >70% white spruce",
"Black spruce stands",
"Pine stands",
"Stands with <70% deciduous and <70% spruce",
"Clearcuts and burns that were salvage-logged",
"Well sites, large cutlines, roads"))

write.csv(VernierTable1, file.path(getwd(), "modules", "vegMapToStrataMap","data", "VernierTable1.csv"))

LCCTable <- data.frame(
    Class = c(
  seq(1:39)
  ),Description = c(
  "Temperate or subpolar needle-leaved evergreen closed tree canopy",
  "Cold deciduous closed tree canopy",
  "Mixed needle-leaved evergreen – cold deciduous closed tree canopy",
  "Mixed needle-leaved evergreen – cold deciduous closed young tree canopy",
  "Mixed cold deciduous – needle-leaved evergreen closed tree canopy",
  "Temperate or subpolar needle-leaved evergreen medium density, moss-shrub understory",
  "Temperate or subpolar needle-leaved evergreen medium density, lichen-shrub understory",
  "Temperate or subpolar needle-leaved evergreen low density, shrub-moss understory",
  "Temperate or subpolar needle-leaved evergreen low density, lichen (rock) understory",
  "Temperate or subpolar needle-leaved evergreen low density, poorly drained",
  "Cold deciduous broad-leaved, low to medium density",
  "Cold deciduous broad-leaved, medium density, young regenerating",
  "Mixed needle-leaved evergreen – cold deciduous, low to medium density",
  "Mixed cold deciduous - needle-leaved evergreen, low to medium density",
  "Low regenerating young mixed cover",
  "High-low shrub dominated",
  "Grassland",
  "Herb-shrub-bare cover",
  "Wetlands",
  "Sparse needle-leaved evergreen, herb-shrub cover",
  "Polar grassland, herb-shrub",
  "Shrub-herb-lichen-bare",
  "Herb-shrub poorly drained",
  "Lichen-shrub-herb-bare soil",
  "Low vegetation cover",
  "Cropland-woodland",
  "High biomass cropland",
  "Medium biomass cropland",
  "Low biomass cropland",
  "Lichen barren",
  "Lichen-sedge-moss-low shrub wetland",
  "Lichen-spruce bog",
  "Rock outcrops",
  "Recent burns",
  "Old burns (<10 years)",
  "Urban and Built-up",
  "Water bodies",
  "Mixes of water and land",
  "Snow/ ice"
  )
)

write.csv(LCCTable, file.path(getwd(), "modules", "vegMapToStrataMap","data", "LCCTable.csv"))
