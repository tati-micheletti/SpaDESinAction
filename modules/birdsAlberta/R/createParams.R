# Covariates Table
# from Vernier et al., 2008

createParams <- function(){ # [ IMPROVE ] BY ADDING THE WHOLE TABLE OF COVARIATES FROM VERNIER AND PASSING SPECIES ARGUMENT

covars <- list(COWA = data.table::data.table(Constant = -3.012,  
                                               L_CUT = -0.526,
                                               L_MIX = 0.638,
                                               L_ODEC = 0.866,
                                               L_YDEC = 1.376, 
                                               L_CC = -0.002,
                                               L_WDIS = 0.101,
                                               N_DEC =  0.934,
                                               N_MIX = 0.472,
                                               N_CUT = 0.645,
                                               N_LATE = -1.119,
                                               N_SB = 0.369,
                                               N_RICH = 0.058),
               RBNU = data.table::data.table(Constant = -2.549,
                                               L_CUT = -2.679,
                                               L_MIX = 0.499,
                                               L_ODEC = -0.135,
                                               L_YDEC = -0.444,
                                               L_CC = 0.004,
                                               L_WDIS = 0.164,
                                               N_DEC = -0.486,
                                               N_MIX = 0.130,
                                               N_CUT = 0.224,  
                                               N_LATE = 1.140,  
                                               N_SB = 0.210, 
                                               N_RICH = 0.164)
               )
               
               
               
      # covar <- data.table::data.table(species = "COWA", "MOWA", "RBNU", "TEWA", "YWAR", 
      #                            Constant = c(-3.012, -3.271, -2.549, 1.138, -4.190),
      #                            L_CUT = c(-0.526, 0.521, -2.679, 0.069, 0.357),
      #                            L_MIX = c(0.638, 0.416, 0.499, 0.497, -0.453),
      #                            L_ODEC = c(0.866, 1.752, -0.135, 0.973, 1.372),
      #                            L_YDEC = c(1.376, 1.518, -0.444, 0.218, 0.508),
      #                            L_CC = c(-0.002, -0.014, 0.004, -0.009, -0.015),
      #                            L_WDIS = c(0.101, -0.124, 0.164, -0.013, -0.163),
      #                            N_DEC = c(0.934, 2.088, -0.486, 1.138, 2.606),
      #                            N_MIX = c(0.472, -0.196, 0.130, -0.066, -0.022),
      #                            N_CUT = c(0.645, 1.187, 0.224, 0.518, 2.869),
      #                            N_LATE = c(-1.119, 1.085, 1.140, -1.330, 1.379),
      #                            N_SB = c(0.369, -0.351, 0.210, 0.140, -0.674),
      #                            N_RICH = c(0.058, 0.163, 0.164, -0.131, 0.067))

message("Loading covariate table for COWA and RBNU. [see Vernier et al. 2008, Appendix 1]")

return(covars)

}