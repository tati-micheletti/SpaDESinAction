plotCaribou <- function(sim = sim, SorensenStats = sim$SorensenStats){

  
populationCaribou <- data.frame(Time = 0:time(sim), lambda = SorensenStats[,1])

plotsCarib <- ggplot2::ggplot(data = populationCaribou, aes(x = Time, y = lambda)) +
                      geom_line()

# growthCaribou <- data.frame(Time = 1:time(sim), lambda = SorensenStats[,1])
#   
#   plotsCarib <- list(populationCaribou = populationCaribou, growthCaribou = growthCaribou)
  
  return(plotsCarib)

}
